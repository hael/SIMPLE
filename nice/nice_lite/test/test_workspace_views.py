import hashlib
import json
from types import SimpleNamespace
from unittest.mock import patch

from django.http import HttpResponse
from django.test import RequestFactory
from django.test import SimpleTestCase

from ..views import workspace_views


class _AuthUser:
    is_authenticated = True
    username = "tester"


class _FakeJob:
    classification_2D_stats = {}


class _FakeQueryset:
    def __init__(self, values_payload):
        self._values_payload = values_payload
        self._jobs = [_FakeJob()]

    def order_by(self, *_args, **_kwargs):
        return self

    def values(self):
        return self._values_payload

    def __iter__(self):
        return iter(self._jobs)


class WorkspaceJobsViewTests(SimpleTestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def test_workspace_jobs_returns_empty_jobs_when_inaccessible(self):
        request = self.factory.get("/workspacejobs")
        request.user = _AuthUser()

        with patch.object(workspace_views, "get_workspace_id", return_value=1), patch.object(workspace_views, "get_project_id", return_value=2), patch.object(workspace_views, "Workspace", return_value=SimpleNamespace()), patch.object(workspace_views, "_is_workspace_accessible", return_value=False), patch.object(workspace_views, "render", return_value=HttpResponse("jobs")) as mock_render:
            response = workspace_views.view_workspace_jobs(request)

        self.assertEqual(response.status_code, 200)
        mock_render.assert_called_once_with(request, "jobs.html", {"jobs": []})

    def test_workspace_jobs_renders_and_sets_checksum_cookie_on_first_load(self):
        request = self.factory.get("/workspacejobs")
        request.user = _AuthUser()

        fake_workspace = SimpleNamespace(id=1)
        fake_queryset = _FakeQueryset([{"id": 1, "status": "running"}])

        with patch.object(workspace_views, "get_workspace_id", return_value=1), patch.object(workspace_views, "get_project_id", return_value=2), patch.object(workspace_views, "Workspace", return_value=fake_workspace), patch.object(workspace_views, "_is_workspace_accessible", return_value=True), patch.object(workspace_views.JobModel.objects, "filter", return_value=fake_queryset), patch.object(workspace_views, "render", return_value=HttpResponse("jobs")) as mock_render, patch.object(workspace_views, "_normalize_latest_cls2d") as mock_normalize:
            response = workspace_views.view_workspace_jobs(request)

        self.assertEqual(response.status_code, 200)
        self.assertIn("workspace_jobs_checksum", response.cookies)
        mock_render.assert_called_once()
        mock_normalize.assert_called_once_with(fake_queryset)

    def test_workspace_jobs_returns_204_when_checksum_matches(self):
        payload = [{"id": 1, "status": "running"}]
        checksum = hashlib.md5(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()

        request = self.factory.get("/workspacejobs")
        request.user = _AuthUser()
        request.COOKIES["workspace_jobs_checksum"] = checksum

        fake_workspace = SimpleNamespace(id=1)
        fake_queryset = _FakeQueryset(payload)

        with patch.object(workspace_views, "get_workspace_id", return_value=1), patch.object(workspace_views, "get_project_id", return_value=2), patch.object(workspace_views, "Workspace", return_value=fake_workspace), patch.object(workspace_views, "_is_workspace_accessible", return_value=True), patch.object(workspace_views.JobModel.objects, "filter", return_value=fake_queryset), patch.object(workspace_views, "render") as mock_render, patch.object(workspace_views, "_normalize_latest_cls2d") as mock_normalize:
            response = workspace_views.view_workspace_jobs(request)

        self.assertEqual(response.status_code, 204)
        mock_render.assert_not_called()
        mock_normalize.assert_not_called()
