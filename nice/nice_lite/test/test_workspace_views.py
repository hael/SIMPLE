import hashlib
import json
import os
import tempfile
from types import SimpleNamespace
from unittest.mock import Mock, patch

from django.http import HttpResponse
from django.test import RequestFactory
from django.test import SimpleTestCase

from ..views import workspace_views


class _AuthUser:
    is_authenticated = True
    username = "tester"


class _FakeJob:
    def __init__(self, status="running", master_stats=None, pckg=""):
        self.classification_2D_stats = {}
        self.particle_sieving_stats = {}
        self.master_stats = master_stats or {}
        self.status = status
        self.pckg = pckg


class _FakeQueryset:
    def __init__(self, values_payload, jobs=None):
        self._values_payload = values_payload
        self._jobs = jobs if jobs is not None else [_FakeJob()]

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
        mock_render.assert_called_once_with(request, "jobs_cards.html", {"jobs": []})

    def test_workspace_jobs_renders_and_sets_checksum_cookie_on_first_load(self):
        request = self.factory.get("/workspacejobs")
        request.user = _AuthUser()

        fake_workspace = SimpleNamespace(id=1, get_absdir=lambda: "/workspace")
        fake_queryset = _FakeQueryset([{"id": 1, "status": "running"}])

        with patch.object(workspace_views, "get_workspace_id", return_value=1), patch.object(workspace_views, "get_project_id", return_value=2), patch.object(workspace_views, "Workspace", return_value=fake_workspace), patch.object(workspace_views, "_is_workspace_accessible", return_value=True), patch.object(workspace_views.JobModel.objects, "filter", return_value=fake_queryset), patch.object(workspace_views, "render", return_value=HttpResponse("jobs")) as mock_render, patch.object(workspace_views, "_normalize_latest_cls2d") as mock_normalize:
            response = workspace_views.view_workspace_jobs(request)

        self.assertEqual(response.status_code, 200)
        self.assertIn("workspace_jobs_checksum", response.cookies)
        mock_render.assert_called_once_with(
            request,
            "jobs_cards.html",
            {"jobs": fake_queryset},
        )
        mock_normalize.assert_called_once_with(fake_queryset)

    def test_workspace_jobs_returns_204_when_checksum_matches(self):
        payload = [{"id": 1, "status": "running"}]
        checksum_payload = {
            "jobs": payload,
            "template": "jobs_cards.html",
            "batch_project_sources": [],
        }
        checksum = hashlib.md5(json.dumps(checksum_payload, sort_keys=True, default=str).encode()).hexdigest()

        request = self.factory.get("/workspacejobs")
        request.user = _AuthUser()
        request.COOKIES["workspace_jobs_checksum"] = checksum

        fake_workspace = SimpleNamespace(id=1, get_absdir=lambda: "/workspace")
        fake_queryset = _FakeQueryset(payload)

        with patch.object(workspace_views, "get_workspace_id", return_value=1), patch.object(workspace_views, "get_project_id", return_value=2), patch.object(workspace_views, "Workspace", return_value=fake_workspace), patch.object(workspace_views, "_is_workspace_accessible", return_value=True), patch.object(workspace_views.JobModel.objects, "filter", return_value=fake_queryset), patch.object(workspace_views, "render") as mock_render, patch.object(workspace_views, "_normalize_latest_cls2d") as mock_normalize:
            response = workspace_views.view_workspace_jobs(request)

        self.assertEqual(response.status_code, 204)
        mock_render.assert_not_called()
        mock_normalize.assert_not_called()

    def test_workspace_jobs_force_renders_when_checksum_matches(self):
        payload = [{"id": 1, "status": "running"}]
        checksum_payload = {
            "jobs": payload,
            "template": "jobs_cards.html",
            "batch_project_sources": [],
        }
        checksum = hashlib.md5(json.dumps(checksum_payload, sort_keys=True, default=str).encode()).hexdigest()

        request = self.factory.get("/workspacejobs", {"force": "1"})
        request.user = _AuthUser()
        request.COOKIES["workspace_jobs_checksum"] = checksum

        fake_workspace = SimpleNamespace(id=1, get_absdir=lambda: "/workspace")
        fake_queryset = _FakeQueryset(payload)

        with patch.object(workspace_views, "get_workspace_id", return_value=1), patch.object(workspace_views, "get_project_id", return_value=2), patch.object(workspace_views, "Workspace", return_value=fake_workspace), patch.object(workspace_views, "_is_workspace_accessible", return_value=True), patch.object(workspace_views.JobModel.objects, "filter", return_value=fake_queryset), patch.object(workspace_views, "render", return_value=HttpResponse("jobs")) as mock_render, patch.object(workspace_views, "_normalize_latest_cls2d") as mock_normalize:
            response = workspace_views.view_workspace_jobs(request)

        self.assertEqual(response.status_code, 200)
        self.assertIn("workspace_jobs_checksum", response.cookies)
        mock_render.assert_called_once_with(
            request,
            "jobs_cards.html",
            {"jobs": fake_queryset},
        )
        mock_normalize.assert_called_once_with(fake_queryset)

    def test_workspace_jobs_annotates_only_finished_batch_projects(self):
        def job(job_id, **overrides):
            values = {
                "id": job_id,
                "disp": job_id,
                "name": f"job {job_id}",
                "pckg": "simple",
                "prog": "demo",
                "status": "finished",
                "dirc": f"{job_id}_job",
            }
            values.update(overrides)
            return SimpleNamespace(**values)

        with tempfile.TemporaryDirectory() as workspace_dir:
            finished_dir = os.path.join(workspace_dir, "1_import_movies")
            running_dir = os.path.join(workspace_dir, "2_preprocess")
            os.mkdir(finished_dir)
            os.mkdir(running_dir)
            finished_project = os.path.join(finished_dir, "workspace.simple")
            with open(finished_project, "w", encoding="utf-8"):
                pass
            with open(os.path.join(running_dir, "workspace.simple"), "w", encoding="utf-8"):
                pass

            finished = job(1, dirc="1_import_movies")
            running = job(2, status="running", dirc="2_preprocess")
            stream = job(3, pckg="stream")
            missing = job(4, pckg="single", dirc="4_missing")
            unsafe = job(5, dirc="../outside")
            jobs = [finished, running, stream, missing, unsafe]

            sources = workspace_views._annotate_batch_project_drag_paths(
                jobs,
                workspace_dir,
            )

        self.assertEqual(sources, [(1, os.path.realpath(finished_project))])
        self.assertEqual(finished.batch_project_drag_path, os.path.realpath(finished_project))
        for unavailable in jobs[1:]:
            self.assertEqual(unavailable.batch_project_drag_path, "")


class WorkspaceJobRefreshTests(SimpleTestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def test_remove_missing_job_records_deletes_only_missing_safe_paths(self):
        with tempfile.TemporaryDirectory() as workspace_dir:
            os.mkdir(os.path.join(workspace_dir, "1_present"))
            reconcile_job_counter = Mock()
            workspace = SimpleNamespace(
                id=4,
                get_absdir=lambda: workspace_dir,
                reconcile_job_counter=reconcile_job_counter,
            )
            jobs = [
                SimpleNamespace(id=1, dirc="1_present"),
                SimpleNamespace(id=2, dirc="2_missing"),
                SimpleNamespace(id=3, dirc="../outside"),
            ]
            delete_queryset = Mock()

            with patch.object(
                workspace_views.JobModel.objects,
                "filter",
                side_effect=[jobs, delete_queryset],
            ) as job_filter:
                removed = workspace_views._remove_missing_job_records(workspace)

        self.assertEqual(removed, 1)
        job_filter.assert_any_call(dset=4, id__in=[2])
        delete_queryset.delete.assert_called_once_with()
        reconcile_job_counter.assert_called_once_with()

    def test_refresh_removes_missing_records_for_owned_workspace(self):
        request = self.factory.post("/refreshworkspacejobs")
        request.user = _AuthUser()
        workspace = SimpleNamespace(id=4)

        with patch.object(workspace_views, "get_workspace_id", return_value=4), patch.object(workspace_views, "get_project_id", return_value=3), patch.object(workspace_views, "Workspace", return_value=workspace), patch.object(workspace_views, "_is_workspace_accessible", return_value=True), patch.object(workspace_views, "_remove_missing_job_records", return_value=2) as remove_missing, patch.object(workspace_views, "clear_checksum_cookies") as clear_checksums, patch.object(workspace_views.messages, "add_message"):
            response = workspace_views.view_refresh_workspace_jobs(request)

        self.assertEqual(response.status_code, 302)
        remove_missing.assert_called_once_with(workspace)
        clear_checksums.assert_called_once_with(request, response)

    def test_refresh_rejects_inaccessible_workspace(self):
        request = self.factory.post("/refreshworkspacejobs")
        request.user = _AuthUser()

        with patch.object(workspace_views, "get_workspace_id", return_value=4), patch.object(workspace_views, "get_project_id", return_value=3), patch.object(workspace_views, "Workspace", return_value=SimpleNamespace(id=4)), patch.object(workspace_views, "_is_workspace_accessible", return_value=False), patch.object(workspace_views, "_remove_missing_job_records") as remove_missing, patch.object(workspace_views.messages, "add_message"):
            response = workspace_views.view_refresh_workspace_jobs(request)

        self.assertEqual(response.status_code, 302)
        remove_missing.assert_not_called()

    def test_refresh_requires_post(self):
        request = self.factory.get("/refreshworkspacejobs")
        request.user = _AuthUser()

        response = workspace_views.view_refresh_workspace_jobs(request)

        self.assertEqual(response.status_code, 405)
