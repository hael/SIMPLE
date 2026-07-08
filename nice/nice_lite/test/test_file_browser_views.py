import tempfile
from types import SimpleNamespace
from unittest.mock import patch

from django.http import HttpResponse
from django.template.loader import render_to_string
from django.test import RequestFactory
from django.test import SimpleTestCase

from ..views import file_browser_views


class _AuthUser:
    is_authenticated = True
    username = "tester"


class _ProjectIdQuery:
    def __init__(self, ids):
        self._ids = ids

    def distinct(self):
        return self

    def values_list(self, *_args, **_kwargs):
        return self._ids


class FileBrowserViewTests(SimpleTestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def _request(self):
        request = self.factory.get("/filebrowser/file")
        request.user = _AuthUser()
        return request

    def test_invalid_project_selection_sets_error(self):
        request = self._request()

        with patch.object(file_browser_views, "get_project_id", return_value=None), patch.object(file_browser_views.ProjectModel.objects, "filter", return_value=_ProjectIdQuery([1])), patch.object(file_browser_views, "render", return_value=HttpResponse("ok")) as mock_render:
            response = file_browser_views.view_file_browser(request, "file")

        self.assertEqual(response.status_code, 200)
        context = mock_render.call_args[0][2]
        self.assertTrue(context["error"])
        self.assertEqual(context["errortext"], "invalid project selection")

    def test_file_browser_selected_entry_style_overrides_hover(self):
        html = render_to_string(
            "filebrowser.html",
            {
                "type": "dir",
                "purpose": "external_input",
                "path": "/tmp",
                "parentdir": "/",
                "error": False,
                "errortext": "",
                "files": [],
                "dirs": ["movies"],
            },
        )

        self.assertIn(".browser-entry-selected:hover", html)
        self.assertIn('classList.add("browser-entry-selected")', html)

    def test_file_selector_directory_and_file_entries_share_hover_highlight(self):
        html = render_to_string(
            "filebrowser.html",
            {
                "type": "file",
                "purpose": "",
                "path": "/tmp",
                "parentdir": "/",
                "error": False,
                "errortext": "",
                "files": ["movie.mrc"],
                "dirs": ["movies"],
            },
        )

        self.assertRegex(
            html,
            r'class="[^"]*browser-entry[^"]*hover:bg-streambg[^"]*cursor-pointer[^"]*"\s+onclick="selectDir',
        )
        self.assertRegex(
            html,
            r'class="[^"]*browser-entry[^"]*hover:bg-streambg[^"]*cursor-pointer[^"]*"\s+onclick="selectFile',
        )
        self.assertIn('selectDir(this,', html)
        self.assertIn('selectBrowserEntry(element);\n            watchDouble += 1;', html)

    def test_project_root_picker_allows_no_project_selection(self):
        with tempfile.TemporaryDirectory() as base_dir:
            request = self.factory.get("/filebrowser/dir", {"purpose": "project_root", "selectedpath": base_dir})
            request.user = _AuthUser()

            with patch.object(file_browser_views, "get_project_id", return_value=None), patch.object(file_browser_views.ProjectModel.objects, "filter", return_value=_ProjectIdQuery([])), patch.object(file_browser_views, "render", return_value=HttpResponse("ok")) as mock_render:
                response = file_browser_views.view_file_browser(request, "dir")

        self.assertEqual(response.status_code, 200)
        context = mock_render.call_args[0][2]
        self.assertFalse(context["error"])
        self.assertEqual(context["purpose"], "project_root")
        self.assertEqual(context["path"], base_dir)

    def test_external_input_directory_picker_allows_outside_project(self):
        with tempfile.TemporaryDirectory() as project_dir, tempfile.TemporaryDirectory() as source_dir:
            request = self.factory.get("/filebrowser/dir", {"purpose": "external_input", "selectedpath": source_dir})
            request.user = _AuthUser()

            with patch.object(file_browser_views, "get_project_id", return_value=1), patch.object(file_browser_views.ProjectModel.objects, "filter", return_value=_ProjectIdQuery([1])), patch.object(file_browser_views, "Project", return_value=SimpleNamespace(absdir=project_dir)), patch.object(file_browser_views, "render", return_value=HttpResponse("ok")) as mock_render:
                response = file_browser_views.view_file_browser(request, "dir")

        self.assertEqual(response.status_code, 200)
        context = mock_render.call_args[0][2]
        self.assertFalse(context["error"])
        self.assertEqual(context["purpose"], "external_input")
        self.assertEqual(context["path"], source_dir)

    def test_path_outside_project_is_rejected(self):
        request = self._request()
        with tempfile.TemporaryDirectory() as base_dir:
            with patch.object(file_browser_views, "get_project_id", return_value=1), patch.object(file_browser_views.ProjectModel.objects, "filter", return_value=_ProjectIdQuery([1])), patch.object(file_browser_views, "Project", return_value=SimpleNamespace(absdir=base_dir)), patch.object(file_browser_views, "render", return_value=HttpResponse("ok")) as mock_render:
                response = file_browser_views.view_file_browser(request, "file", path="../../etc")

        self.assertEqual(response.status_code, 200)
        context = mock_render.call_args[0][2]
        self.assertTrue(context["error"])
        self.assertEqual(context["errortext"], "path outside project")

    def test_directory_listing_hides_dotfiles_and_sorts(self):
        request = self._request()
        with tempfile.TemporaryDirectory() as base_dir:
            # Visible entries.
            open(base_dir + "/b.txt", "w", encoding="utf-8").close()
            open(base_dir + "/a.dat", "w", encoding="utf-8").close()
            os_dir = base_dir + "/z_folder"
            hidden_file = base_dir + "/.hidden.txt"
            hidden_dir = base_dir + "/.hidden_dir"

            import os
            os.makedirs(os_dir, exist_ok=True)
            os.makedirs(hidden_dir, exist_ok=True)
            open(hidden_file, "w", encoding="utf-8").close()

            with patch.object(file_browser_views, "get_project_id", return_value=1), patch.object(file_browser_views.ProjectModel.objects, "filter", return_value=_ProjectIdQuery([1])), patch.object(file_browser_views, "Project", return_value=SimpleNamespace(absdir=base_dir)), patch.object(file_browser_views, "render", return_value=HttpResponse("ok")) as mock_render:
                response = file_browser_views.view_file_browser(request, "file", path=base_dir)

        self.assertEqual(response.status_code, 200)
        context = mock_render.call_args[0][2]
        self.assertFalse(context["error"])
        self.assertEqual(context["files"], ["a.dat", "b.txt"])
        self.assertEqual(context["dirs"], ["z_folder"])
