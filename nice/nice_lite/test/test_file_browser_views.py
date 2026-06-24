import tempfile
from types import SimpleNamespace
from unittest.mock import patch

from django.http import HttpResponse
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
