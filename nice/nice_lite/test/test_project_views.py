from unittest.mock import Mock
from unittest.mock import patch

from django.test import RequestFactory
from django.test import SimpleTestCase

from ..views import project_views


class _AuthUser:
    is_authenticated = True
    username = "tester"


class ProjectViewTests(SimpleTestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def test_create_project_redirects_to_index_shell(self):
        request = self.factory.post(
            "/createproject",
            {
                "new_project_name": "test project",
                "new_project_dirc": "/tmp",
            },
        )
        request.user = _AuthUser()

        project = Mock()
        project.id = 11
        project.new.return_value = True

        workspace = Mock()
        workspace.id = 22
        workspace.new.return_value = True

        with patch.object(project_views, "Project", return_value=project), patch.object(project_views, "Workspace", return_value=workspace), patch.object(project_views.os.path, "isdir", return_value=True), patch.object(project_views, "clear_checksum_cookies"):
            response = project_views.view_create_project(request)

        self.assertEqual(response.status_code, 302)
        self.assertEqual(response["Location"], "/")
        self.assertEqual(response.cookies["selected_project_id"].value, "11")
        self.assertEqual(response.cookies["selected_workspace_id"].value, "22")

    def test_close_new_project_preserves_selection_and_returns_to_index_shell(self):
        request = self.factory.get("/closeprojectform")
        request.user = _AuthUser()
        request.COOKIES["selected_project_id"] = "7"
        request.COOKIES["selected_workspace_id"] = "9"

        with patch.object(project_views, "clear_checksum_cookies"):
            response = project_views.view_close_new_project(request)

        self.assertEqual(response.status_code, 302)
        self.assertEqual(response["Location"], "/")
        self.assertNotIn("selected_project_id", response.cookies)
        self.assertNotIn("selected_workspace_id", response.cookies)
