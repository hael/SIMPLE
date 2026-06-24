from types import SimpleNamespace
from unittest.mock import Mock
from unittest.mock import patch

from django.http import HttpResponse
from django.test import RequestFactory
from django.test import SimpleTestCase

from ..views import index_views
from ..views import job_builder_views
from ..views import workspace_views


class _AuthUser:
    is_authenticated = True
    username = "tester"


class _FakeProjectQuery:
    def __init__(self, ids):
        self._ids = ids

    def distinct(self):
        return self

    def values_list(self, *_args, **_kwargs):
        return self._ids

    def __len__(self):
        return len(self._ids)


def _render_with_context(_request, _template, context):
    response = HttpResponse("ok")
    response._ctx = context
    return response


def _reverse_with_query(name, *args, **kwargs):
    result = f"rev:{name}"
    query = kwargs.get("query")
    if isinstance(query, dict) and len(query) > 0:
        query_items = [f"{key}={value}" for key, value in query.items()]
        result += "?" + "&".join(query_items)
    return result


class IndexViewBranchTests(SimpleTestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def test_project_sentinel_routes_to_new_project(self):
        request = self.factory.get("/")
        request.user = _AuthUser()

        with patch.object(index_views, "get_project_id", return_value=-1), patch.object(index_views, "get_workspace_id", return_value=None), patch.object(index_views.ProjectModel.objects, "filter", return_value=_FakeProjectQuery([1])), patch.object(index_views.WorkspaceModel.objects, "filter", return_value=[]), patch.object(index_views, "reverse", side_effect=_reverse_with_query), patch.object(index_views, "render", side_effect=_render_with_context), patch.object(index_views, "clear_checksum_cookies"), patch.object(index_views.messages, "add_message"):
            response = index_views.view_index(request)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response._ctx["iframeurl"], "rev:nice_lite:new_project")

    def test_workspace_sentinel_creates_workspace_when_project_accessible(self):
        request = self.factory.get("/")
        request.user = _AuthUser()

        workspace_obj = Mock()
        workspace_obj.new.return_value = True
        workspace_obj.get_id.return_value = 42

        with patch.object(index_views, "get_project_id", return_value=1), patch.object(index_views, "get_workspace_id", return_value=-1), patch.object(index_views.ProjectModel.objects, "filter", return_value=_FakeProjectQuery([1])), patch.object(index_views.WorkspaceModel.objects, "filter", return_value=["w1"]), patch.object(index_views, "Project", return_value=object()), patch.object(index_views, "Workspace", return_value=workspace_obj), patch.object(index_views, "reverse", side_effect=_reverse_with_query), patch.object(index_views, "render", side_effect=_render_with_context), patch.object(index_views, "clear_checksum_cookies"), patch.object(index_views.messages, "add_message"):
            response = index_views.view_index(request)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response._ctx["current_workspace_id"], 42)
        self.assertEqual(response._ctx["iframeurl"], "rev:nice_lite:workspace?selected_workspace_id=42")


class JobBuilderBranchTests(SimpleTestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def test_inaccessible_selected_job_clears_cookie(self):
        request = self.factory.get("/jobbuilder")
        request.user = _AuthUser()

        stream_job = Mock()
        stream_job.get_jobmodel.return_value = object()

        simple_stream = Mock()
        simple_stream.loadUIJSON.return_value = True
        simple_stream.get_ui.return_value = {"user_inputs": []}

        simple_batch = Mock()
        simple_batch.loadUIJSON.return_value = True
        simple_batch.get_ui.return_value = {}

        with patch.object(job_builder_views, "get_job_id", return_value=77), patch.object(job_builder_views, "StreamJob", return_value=stream_job), patch.object(job_builder_views, "_is_job_accessible", return_value=False), patch.object(job_builder_views, "SIMPLEStream", return_value=simple_stream), patch.object(job_builder_views, "SIMPLEBatch", return_value=simple_batch), patch.object(job_builder_views, "render", return_value=HttpResponse("ok")), patch.object(job_builder_views, "clear_checksum_cookies"), patch.object(job_builder_views.messages, "add_message"):
            response = job_builder_views.view_job_builder(request)

        self.assertEqual(response.status_code, 200)
        self.assertIn("selected_job_id", response.cookies)

    def test_accessible_job_prefills_stream_inputs_from_args(self):
        request = self.factory.get("/jobbuilder")
        request.user = _AuthUser()

        stream_job = Mock()
        stream_job.get_jobmodel.return_value = SimpleNamespace(args={"gain": "2.5"})

        simple_stream = Mock()
        simple_stream.loadUIJSON.return_value = True
        simple_stream.get_ui.return_value = {"user_inputs": [{"key": "gain", "keytype": "float"}]}

        simple_batch = Mock()
        simple_batch.loadUIJSON.return_value = True
        simple_batch.get_ui.return_value = {}

        with patch.object(job_builder_views, "get_job_id", return_value=88), patch.object(job_builder_views, "StreamJob", return_value=stream_job), patch.object(job_builder_views, "_is_job_accessible", return_value=True), patch.object(job_builder_views, "SIMPLEStream", return_value=simple_stream), patch.object(job_builder_views, "SIMPLEBatch", return_value=simple_batch), patch.object(job_builder_views, "clear_checksum_cookies"), patch.object(job_builder_views, "render", side_effect=_render_with_context):
            response = job_builder_views.view_job_builder(request)

        self.assertEqual(response.status_code, 200)
        stream_inputs = response._ctx["stream_user_inputs"]
        self.assertEqual(stream_inputs[0]["value"], "2.5")


class WorkspaceAccessBranchTests(SimpleTestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def test_view_workspace_inaccessible_returns_204(self):
        request = self.factory.get("/workspace")
        request.user = _AuthUser()

        with patch.object(workspace_views, "get_workspace_id", return_value=1), patch.object(workspace_views, "get_project_id", return_value=2), patch.object(workspace_views, "Workspace", return_value=SimpleNamespace()), patch.object(workspace_views, "_is_workspace_accessible", return_value=False), patch.object(workspace_views.messages, "add_message"), patch.object(workspace_views, "print_error"):
            response = workspace_views.view_workspace(request)

        self.assertEqual(response.status_code, 204)
