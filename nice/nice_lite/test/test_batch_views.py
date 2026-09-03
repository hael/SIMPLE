from types import SimpleNamespace
from unittest.mock import Mock, patch

from django.core import signing
from django.http import HttpResponse, HttpResponseRedirect
from django.test import RequestFactory, SimpleTestCase

from ..views import batch_views


class _AuthUser:
    is_authenticated = True
    username = "tester"


class BatchViewTests(SimpleTestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def _request(self, path):
        request = self.factory.post(path, {"selected_job_id": "7"})
        request.user = _AuthUser()
        return request

    def _get_request(self, path):
        request = self.factory.get(path)
        request.user = _AuthUser()
        return request

    def test_batch_access_helper_rejects_job_owned_by_another_user(self):
        batch_job = Mock()
        batch_job.get_jobmodel.return_value = SimpleNamespace(
            dset=SimpleNamespace(user="another-user"),
        )
        with patch.object(batch_views, "BatchJob", return_value=batch_job):
            resolved_job, resolved_model = batch_views._get_accessible_batch_job(
                self._get_request("/viewbatch/7"),
                "view_batch",
                job_id=7,
            )

        self.assertIsNone(resolved_job)
        self.assertIsNone(resolved_model)

    def test_batch_access_helper_rejects_non_batch_job(self):
        batch_job = Mock()
        batch_job.get_jobmodel.return_value = None
        with patch.object(batch_views, "BatchJob", return_value=batch_job):
            resolved_job, resolved_model = batch_views._get_accessible_batch_job(
                self._get_request("/viewbatch/7"),
                "view_batch",
                job_id=7,
            )

        self.assertIsNone(resolved_job)
        self.assertIsNone(resolved_model)

    def test_batch_detail_renders_owned_job_and_sets_selection_cookies(self):
        project = SimpleNamespace(id=3, name="project")
        workspace = SimpleNamespace(
            id=4,
            name="workspace",
            user="tester",
            proj=project,
            proj_id=project.id,
        )
        jobmodel = SimpleNamespace(
            id=7,
            disp=2,
            name="Import Movie Data",
            desc="movies",
            status="finished",
            cdat="created",
            args={"nthr": "4"},
            master_stats={
                "job_type": "batch",
                "package": "simple",
                "program": "import_movies",
            },
            dset=workspace,
            dset_id=workspace.id,
        )
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = "/project/workspace/2_import_movies/workspace.simple"
        batch_job.get_artifact_summary.return_value = {"counts": [], "images": []}
        batch_job.get_safe_job_dir.return_value = "/project/workspace/2_import_movies"
        batch_job.get_log_tails.return_value = []
        launcher = Mock()
        launcher.get_ui.return_value = {
            "import_movies": {
                "program": {"executable": "simple_exec"},
                "compute": [
                    {"key": "nthr", "label": "Number of threads", "has_default": True, "default": 8},
                    {"key": "scale", "label": "Scale", "has_default": True, "default": 1.5},
                    {"key": "mskdiam", "label": "Mask diameter"},
                    {"key": "projfile", "label": "Project file"},
                ],
            },
        }
        rendered_response = HttpResponse("batch")
        request = self._get_request("/viewbatch/7")
        request.COOKIES["workspace_checksum"] = "stale"

        with (
            patch.object(batch_views, "_get_accessible_batch_job", return_value=(batch_job, jobmodel)),
            patch.object(batch_views, "SIMPLEBatch", return_value=launcher),
            patch.object(batch_views, "SIMPLEProjFile") as projfile,
            patch.object(batch_views, "render", return_value=rendered_response) as render,
        ):
            projfile.return_value.getGlobalStats.return_value = {
                "mic": {"n": 24, "info": "micrographs"},
            }
            response = batch_views.view_batch(request, 7)

        self.assertEqual(response.status_code, 200)
        context = render.call_args.args[2]
        self.assertEqual(context["name"], "Import Movie Data")
        self.assertEqual(context["project_sections"][0]["count"], 24)
        self.assertEqual(context["arguments"], [
            {
                "key": "nthr",
                "label": "Number of threads",
                "value": "4",
                "origin": "submitted",
                "submitted": True,
            },
            {
                "key": "scale",
                "label": "Scale",
                "value": 1.5,
                "origin": "default",
                "submitted": False,
            },
            {
                "key": "mskdiam",
                "label": "Mask diameter",
                "value": None,
                "origin": "unset",
                "submitted": False,
            },
        ])
        self.assertEqual(context["submitted_argument_count"], 1)
        self.assertEqual(response.cookies["selected_project_id"].value, "3")
        self.assertEqual(response.cookies["selected_workspace_id"].value, "4")
        self.assertEqual(response.cookies["workspace_checksum"]["max-age"], 0)

    def test_argument_rows_fall_back_to_saved_keys(self):
        launcher = Mock()
        launcher.get_ui.return_value = None
        jobmodel = SimpleNamespace(args={"nthr": "4"})

        with patch.object(batch_views, "SIMPLEBatch", return_value=launcher):
            arguments = batch_views._argument_rows(
                jobmodel,
                {"package": "simple", "program": "import_movies"},
            )

        self.assertEqual(arguments, [{
            "key": "nthr",
            "label": "nthr",
            "value": "4",
            "origin": "submitted",
            "submitted": True,
        }])

    def test_extract_batch_context_paginates_particle_stack_headers(self):
        project = SimpleNamespace(name="project")
        workspace = SimpleNamespace(name="workspace", proj=project)
        jobmodel = SimpleNamespace(
            id=7,
            disp=1,
            name="Extract Particle Images",
            desc="",
            status="finished",
            cdat="created",
            args={},
            master_stats={"package": "unknown", "program": "extract"},
            dset=workspace,
        )
        particle_stack_page = {
            "stacks": [{"name": "particles.mrcs", "count": 80}],
            "particles": [{"number": 41, "stack_name": "particles.mrcs", "stack_index": 41}],
            "total": 80,
        }
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = None
        batch_job.get_artifact_summary.return_value = {"counts": [], "images": []}
        batch_job.get_particle_stack_page.return_value = particle_stack_page
        batch_job.get_safe_job_dir.return_value = "/project/workspace/1_extract"
        batch_job.get_log_tails.return_value = []

        context = batch_views._batch_detail_context(batch_job, jobmodel, particle_page=2)

        batch_job.get_particle_stack_page.assert_called_once_with(page=2, page_size=40)
        self.assertIs(context["particle_stack_page"], particle_stack_page)
        self.assertEqual(context["artifact_counts"], [{"extension": "MRCS", "count": 1}])

    def test_particle_thumbnail_endpoint_renders_only_requested_particle(self):
        batch_job = Mock()
        batch_job.get_particle_thumbnail.return_value = b"png"
        request = self._get_request("/batchparticle/7/particles.mrcs/9")

        with patch.object(
            batch_views,
            "_get_accessible_batch_job",
            return_value=(batch_job, object()),
        ):
            response = batch_views.view_batch_particle_thumbnail(
                request,
                7,
                "particles.mrcs",
                9,
            )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "image/png")
        self.assertEqual(response.content, b"png")
        batch_job.get_particle_thumbnail.assert_called_once_with("particles.mrcs", 9)

    def test_import_movie_page_paginates_submitted_movie_source(self):
        batch_job = Mock()
        batch_job.get_import_movie_paths.return_value = [
            f"/data/movies/movie_{number:03d}.mrcs"
            for number in range(1, 482)
        ]

        movie_page = batch_views._import_movie_page(
            7,
            batch_job,
            page=2,
            page_size=40,
        )

        batch_job.get_import_movie_paths.assert_called_once_with()
        self.assertEqual(movie_page["first_movie"], 41)
        self.assertEqual(movie_page["last_movie"], 80)
        self.assertEqual(movie_page["movies"][0]["number"], 41)
        payload = signing.Signer(
            salt=batch_views._BATCH_MOVIE_THUMBNAIL_SALT,
        ).unsign_object(movie_page["movies"][0]["token"])
        self.assertEqual(payload, {
            "job_id": 7,
            "path": "/data/movies/movie_041.mrcs",
            "version": batch_views.BatchJob.MOVIE_THUMBNAIL_CACHE_VERSION,
        })
        self.assertEqual(movie_page["page_numbers"][-1], 13)

    def test_movie_thumbnail_endpoint_renders_only_signed_import_path(self):
        movie_path = "/data/movies/movie_001.mrcs"
        token = batch_views._movie_thumbnail_token(7, movie_path)
        batch_job = Mock()
        batch_job.get_import_movie_thumbnail.return_value = b"webp"
        request = self._get_request(f"/batchmovie/7/{token}")

        with patch.object(
            batch_views,
            "_get_accessible_batch_job",
            return_value=(batch_job, object()),
        ):
            response = batch_views.view_batch_movie_thumbnail(request, 7, token)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "image/webp")
        self.assertEqual(response.content, b"webp")
        batch_job.get_import_movie_thumbnail.assert_called_once_with(movie_path)

    def test_movie_thumbnail_endpoint_rejects_token_for_another_job(self):
        token = batch_views._movie_thumbnail_token(8, "/data/movie.mrcs")
        batch_job = Mock()
        request = self._get_request(f"/batchmovie/7/{token}")

        with patch.object(
            batch_views,
            "_get_accessible_batch_job",
            return_value=(batch_job, object()),
        ):
            response = batch_views.view_batch_movie_thumbnail(request, 7, token)

        self.assertEqual(response.status_code, 404)
        batch_job.get_import_movie_thumbnail.assert_not_called()

    def test_batch_detail_exposes_default_hidden_ctf_micrograph_toggle(self):
        jobmodel = SimpleNamespace(
            id=7,
            disp=2,
            name="Estimate CTF",
            desc="",
            status="running",
            cdat="created",
            args={},
            master_stats={
                "job_type": "batch",
                "package": "unknown",
                "program": "ctf_estimate",
            },
            dset=SimpleNamespace(
                name="workspace",
                proj=SimpleNamespace(name="project"),
            ),
        )
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = None
        artifact_summary = {
            "counts": [],
            "images": [{
                "path": "ctf.jpg",
                "previews": [
                    {"kind": "image"},
                    {"kind": "micrograph", "hidden_by_default": True},
                ],
            }],
        }
        batch_job.get_artifact_summary.return_value = artifact_summary
        batch_job.get_safe_job_dir.return_value = "/workspace/2_ctf_estimate"
        batch_job.get_log_tails.return_value = []

        context = batch_views._batch_detail_context(batch_job, jobmodel)

        batch_job.get_artifact_summary.assert_called_once_with()
        self.assertEqual(context["artifact_images"], artifact_summary["images"])
        self.assertTrue(context["ctf_artifact_micrograph_toggle_available"])

    def test_batch_detail_prefers_pick_previews_from_motion_artifacts(self):
        jobmodel = SimpleNamespace(
            id=7,
            disp=3,
            name="Pick Particles from Micrographs",
            desc="",
            status="finished",
            cdat="created",
            args={},
            master_stats={
                "job_type": "batch",
                "package": "unknown",
                "program": "pick",
            },
            dset=SimpleNamespace(
                name="workspace",
                proj=SimpleNamespace(name="project"),
            ),
        )
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = "/workspace/3_pick/workspace.simple"
        batch_job.get_artifact_summary.return_value = {"counts": [], "images": []}
        batch_job.get_safe_job_dir.return_value = "/workspace/3_pick"
        batch_job.get_log_tails.return_value = []
        batch_job.get_pick_micrograph_previews.return_value = [{
            "path": "/workspace/1_motion/movie_thumb.jpg",
            "number": 11,
            "xdim": 4096,
            "ydim": 3072,
            "boxes": [{"x": 101, "y": 202, "width": 180, "height": 180}],
        }]

        with patch.object(batch_views, "SIMPLEProjFile") as projfile:
            project_reader = projfile.return_value
            project_reader.getGlobalStats.return_value = {"mic": {"n": 30}}

            context = batch_views._batch_detail_context(batch_job, jobmodel)

        batch_job.get_pick_micrograph_previews.assert_called_once_with(
            max_previews=20,
            max_coordinates=1500,
        )
        project_reader.getFieldStats.assert_not_called()
        self.assertEqual(context["pick_micrographs"], [{
            "path": "/workspace/1_motion/movie_thumb.jpg",
            "number": 11,
            "xdim": 4096,
            "ydim": 3072,
            "boxes": [{"x": 101, "y": 202, "width": 180, "height": 180}],
        }])
        self.assertTrue(context["pick_box_overlay_available"])

    def test_batch_detail_rejects_inaccessible_job(self):
        with (
            patch.object(batch_views, "_get_accessible_batch_job", return_value=(None, None)),
            patch.object(batch_views, "render") as render,
            patch.object(batch_views.messages, "add_message"),
        ):
            response = batch_views.view_batch(self._get_request("/viewbatch/7"), 7)

        self.assertEqual(response.status_code, 302)
        render.assert_not_called()

    def test_stop_calls_batch_job_for_owned_running_job(self):
        job = Mock()
        job.stop.return_value = True
        jobmodel = SimpleNamespace(status="running")
        with patch.object(batch_views, "_get_accessible_batch_job", return_value=(job, jobmodel)), patch.object(batch_views.messages, "add_message"), patch.object(batch_views, "redirect", return_value=HttpResponseRedirect("/workspace")):
            response = batch_views.view_batch_stop(self._request("/stopbatch"))

        self.assertEqual(response.status_code, 302)
        job.stop.assert_called_once_with()

    def test_delete_calls_permanent_delete_for_owned_finished_job(self):
        job = Mock()
        job.delete.return_value = True
        jobmodel = SimpleNamespace(
            status="finished",
            dset_id=4,
            dset=SimpleNamespace(proj_id=3),
        )
        with patch.object(batch_views, "_get_accessible_batch_job", return_value=(job, jobmodel)), patch.object(batch_views, "Project") as project_cls, patch.object(batch_views, "Workspace") as workspace_cls, patch.object(batch_views.messages, "add_message"), patch.object(batch_views, "redirect", return_value=HttpResponseRedirect("/workspace")):
            response = batch_views.view_batch_delete(self._request("/deletebatch"))

        self.assertEqual(response.status_code, 302)
        job.delete.assert_called_once_with(project_cls.return_value, workspace_cls.return_value)

    def test_delete_calls_permanent_delete_for_stale_queued_job(self):
        job = Mock()
        job.queued_job_can_delete.return_value = True
        job.delete.return_value = True
        jobmodel = SimpleNamespace(
            status="queued",
            dset_id=4,
            dset=SimpleNamespace(proj_id=3),
        )
        with patch.object(batch_views, "_get_accessible_batch_job", return_value=(job, jobmodel)), patch.object(batch_views, "Project") as project_cls, patch.object(batch_views, "Workspace") as workspace_cls, patch.object(batch_views.messages, "add_message"), patch.object(batch_views, "redirect", return_value=HttpResponseRedirect("/workspace")):
            response = batch_views.view_batch_delete(self._request("/deletebatch"))

        self.assertEqual(response.status_code, 302)
        job.queued_job_can_delete.assert_called_once_with()
        job.delete.assert_called_once_with(project_cls.return_value, workspace_cls.return_value)

    def test_delete_rejects_active_queued_job(self):
        job = Mock()
        job.queued_job_can_delete.return_value = False
        jobmodel = SimpleNamespace(status="queued")
        with patch.object(batch_views, "_get_accessible_batch_job", return_value=(job, jobmodel)), patch.object(batch_views.messages, "add_message"), patch.object(batch_views, "redirect", return_value=HttpResponseRedirect("/workspace")):
            response = batch_views.view_batch_delete(self._request("/deletebatch"))

        self.assertEqual(response.status_code, 302)
        job.queued_job_can_delete.assert_called_once_with()
        job.delete.assert_not_called()

    def test_delete_rejects_running_job(self):
        job = Mock()
        jobmodel = SimpleNamespace(status="running")
        with patch.object(batch_views, "_get_accessible_batch_job", return_value=(job, jobmodel)), patch.object(batch_views.messages, "add_message"), patch.object(batch_views, "redirect", return_value=HttpResponseRedirect("/workspace")):
            response = batch_views.view_batch_delete(self._request("/deletebatch"))

        self.assertEqual(response.status_code, 302)
        job.delete.assert_not_called()
