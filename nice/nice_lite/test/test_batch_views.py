import gzip
import os
import struct
import tempfile
from types import SimpleNamespace
from unittest.mock import Mock, patch

from django.http import HttpResponse, HttpResponseRedirect
from django.test import RequestFactory, SimpleTestCase
from django.urls import resolve, reverse
import mrcfile
import numpy as np

from ..views import batch_views
from ..compress_volume import volume_cache_version


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

    def test_viewbatch_route_renders_active_template_and_sets_selection_cookies(self):
        self.assertEqual(
            reverse("nice_lite:view_batch", kwargs={"jobid": 7}),
            "/viewbatch/7",
        )
        match = resolve("/viewbatch/7")
        self.assertIs(match.func, batch_views.view_batch)

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
            pckg="simple",
            prog="import_movies",
            master_stats={
                "project_metadata": {"nmics": 24},
                "source": {
                    "type": "project_file",
                    "filename": "inputs/particles.simple",
                },
            },
            dset=workspace,
            dset_id=workspace.id,
        )
        batch_job = Mock()
        batch_job.get_absdir.return_value = "/project/workspace/2_import_movies"
        batch_job.get_log_tails.return_value = []
        launcher = Mock()
        launcher.get_ui.return_value = {
            "import_movies": {
                "program": {"executable": "simple_exec"},
                "compute": [
                    {
                        "key": "nthr",
                        "label": "Number of threads",
                        "has_default": True,
                        "default": 8,
                    },
                    {
                        "key": "scale",
                        "label": "Scale",
                        "has_default": True,
                        "default": 1.5,
                    },
                    {"key": "projfile", "label": "Project file"},
                ],
            },
        }
        rendered_response = HttpResponse("batch")
        request = self._get_request("/viewbatch/7")
        request.COOKIES["workspace_checksum"] = "stale"

        with (
            patch.object(
                batch_views,
                "_get_accessible_batch_job",
                return_value=(batch_job, jobmodel),
            ),
            patch.object(batch_views, "SIMPLEBatch", return_value=launcher),
            patch.object(
                batch_views,
                "render",
                return_value=rendered_response,
            ) as render,
        ):
            response = match.func(request, **match.kwargs)

        self.assertEqual(response.status_code, 200)
        self.assertEqual(render.call_args.args[1], "nice_batch/batchview.html")
        context = render.call_args.args[2]
        self.assertEqual(context["name"], "Import Movie Data")
        self.assertEqual(context["jobstats"], {"nmics": 24})
        self.assertEqual(context["folder"], "/project/workspace/2_import_movies")
        self.assertEqual(context["arguments"], [
            {
                "key": "projfile",
                "label": "Project file",
                "value": "inputs/particles.simple",
                "origin": "submitted",
                "submitted": True,
                "visibility": "standard",
            },
            {
                "key": "nthr",
                "label": "Number of threads",
                "value": "4",
                "origin": "submitted",
                "submitted": True,
                "visibility": "standard",
            },
            {
                "key": "scale",
                "label": "Scale",
                "value": 1.5,
                "origin": "default",
                "submitted": False,
                "visibility": "standard",
            },
        ])
        self.assertEqual(context["submitted_argument_count"], 2)
        self.assertEqual(response.cookies["selected_project_id"].value, "3")
        self.assertEqual(response.cookies["selected_workspace_id"].value, "4")
        self.assertEqual(response.cookies["workspace_checksum"]["max-age"], 0)

    def test_argument_rows_fall_back_to_saved_keys(self):
        launcher = Mock()
        launcher.get_ui.return_value = None
        jobmodel = SimpleNamespace(args={"nthr": "4"}, pckg="simple", prog="import_movies")

        with patch.object(batch_views, "SIMPLEBatch", return_value=launcher):
            arguments = batch_views._argument_rows(jobmodel)

        self.assertEqual(arguments, [{
            "key": "nthr",
            "label": "nthr",
            "value": "4",
            "origin": "submitted",
            "submitted": True,
            "visibility": "standard",
        }])

    def test_argument_rows_omit_missing_or_projectless_recorded_sources(self):
        launcher = Mock()
        launcher.get_ui.return_value = None

        for source in (
            None,
            {"type": "none"},
            {"type": "project_file", "filename": ""},
        ):
            with self.subTest(source=source):
                metadata = {} if source is None else {"source": source}
                jobmodel = SimpleNamespace(
                    args={},
                    pckg="simple",
                    prog="import_movies",
                    master_stats=metadata,
                )

                with patch.object(batch_views, "SIMPLEBatch", return_value=launcher):
                    self.assertEqual(batch_views._argument_rows(jobmodel), [])

    def test_argument_rows_resolve_recorded_batch_and_snapshot_projects(self):
        launcher = Mock()
        launcher.get_ui.return_value = None
        cases = (
            (
                {"type": "batch_job", "batch_job_id": 8},
                SimpleNamespace(dirc="7_solve2D", pckg="simple"),
                "7_solve2D/workspace.simple",
            ),
            (
                {
                    "type": "stream_snapshot",
                    "stream_job_id": 5,
                    "particle_set_id": 2,
                    "filename": "snapshot_2.simple",
                },
                SimpleNamespace(
                    dirc="3_simple_stream",
                    pckg="simple_stream",
                    particle_sets_stats={"particle_sets": [{"id": 2, "type": "snapshot2D"}]},
                ),
                "3_simple_stream/classification_2D/snapshots/"
                "snapshot_2/snapshot_2.simple",
            ),
            (
                {
                    "type": "stream_snapshot",
                    "stream_job_id": 5,
                    "particle_set_id": 4,
                    "filename": "snapshot_4.simple",
                },
                SimpleNamespace(
                    dirc="3_simple_stream",
                    pckg="simple_stream",
                    particle_sets_stats={"particle_sets": [{"id": 4, "type": "snapshot3D"}]},
                ),
                "3_simple_stream/solve3D_multistate/snapshots/"
                "snapshot_4/snapshot_4.simple",
            ),
        )

        for source, source_job, expected in cases:
            with self.subTest(source=source):
                jobmodel = SimpleNamespace(
                    args={},
                    pckg="simple",
                    prog="refine2D",
                    dset_id=4,
                    master_stats={"source": source},
                )
                source_jobs = Mock()
                source_jobs.first.return_value = source_job

                with (
                    patch.object(batch_views, "SIMPLEBatch", return_value=launcher),
                    patch.object(
                        batch_views.JobModel.objects,
                        "filter",
                        return_value=source_jobs,
                    ) as filter_jobs,
                ):
                    arguments = batch_views._argument_rows(jobmodel)

                self.assertEqual(arguments, [{
                    "key": "projfile",
                    "label": "Project file",
                    "value": expected,
                    "origin": "submitted",
                    "submitted": True,
                    "visibility": "standard",
                }])
                expected_source_id = source.get(
                    "batch_job_id",
                    source.get("stream_job_id"),
                )
                filter_jobs.assert_called_once_with(id=expected_source_id, dset_id=4)

    def test_argument_rows_omit_project_ignored_by_simple_dispatch(self):
        launcher = Mock()
        launcher.get_ui.return_value = None

        for program in ("new_project", "reproject"):
            with self.subTest(program=program):
                jobmodel = SimpleNamespace(
                    args={},
                    pckg="simple",
                    prog=program,
                    master_stats={
                        "source": {
                            "type": "project_file",
                            "filename": "input.simple",
                        },
                    },
                )

                with patch.object(batch_views, "SIMPLEBatch", return_value=launcher):
                    self.assertEqual(batch_views._argument_rows(jobmodel), [])

    def test_active_batch_volume_viewer_is_explicit_and_program_gated(self):
        jobmodel = SimpleNamespace(
            id=7,
            disp=9,
            name="Initial 3D Reconstruction",
            desc="",
            status="finished",
            cdat="created",
            args={},
            pckg="simple",
            prog="solve3D",
            master_stats={},
            dset=SimpleNamespace(
                name="workspace",
                proj=SimpleNamespace(name="project"),
            ),
        )
        batch_job = Mock()
        batch_job.get_log_tails.return_value = []
        batch_job.get_absdir.return_value = "/workspace/9_solve3D"
        batch_job.get_stage_volume_outputs.return_value = [{
            "path": "/workspace/9_solve3D/recvol_state01.mrc",
            "name": "recvol_state01.mrc",
            "stage": "1",
            "state": 1,
            "kind": "volpath",
            "width": 256,
            "height": 256,
            "depth": 256,
            "voxel_size": (1.3, 1.3, 1.3),
            "minimum": -2.0,
            "maximum": 8.0,
        }]

        with patch.object(batch_views, "_argument_rows", return_value=[]):
            closed_context = batch_views._batch_overview_context(
                batch_job,
                jobmodel,
            )
            open_context = batch_views._batch_overview_context(
                batch_job,
                jobmodel,
                volume_viewer_requested=True,
            )
            jobmodel.prog = "refine3D_states"
            jobmodel.status = "running"
            running_states_context = batch_views._batch_overview_context(
                batch_job,
                jobmodel,
                volume_viewer_requested=True,
            )
            jobmodel.prog = "solve2D"
            wrong_program_context = batch_views._batch_overview_context(
                batch_job,
                jobmodel,
                volume_viewer_requested=True,
            )

        self.assertFalse(closed_context["volume_viewer_requested"])
        self.assertEqual(closed_context["volume_outputs"], [])
        self.assertTrue(open_context["volume_viewer_requested"])
        self.assertEqual(
            open_context["volume_outputs"][0]["name"],
            "recvol_state01.mrc",
        )
        self.assertNotIn("path", open_context["volume_outputs"][0])
        self.assertTrue(running_states_context["volume_viewer_requested"])
        self.assertEqual(
            running_states_context["volume_outputs"][0]["name"],
            "recvol_state01.mrc",
        )
        self.assertFalse(wrong_program_context["volume_viewer_requested"])
        self.assertEqual(wrong_program_context["volume_outputs"], [])

    def test_batch_volume_data_serves_compressed_owned_mrc_file(self):
        temporary_job_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_job_dir.cleanup)
        volume_path = os.path.join(
            temporary_job_dir.name,
            "recvol_state01.mrc",
        )
        with mrcfile.new(volume_path) as volume_file:
            volume_file.set_data(np.arange(64, dtype=np.float32).reshape(4, 4, 4))
            volume_file.voxel_size = 1.5

        jobmodel = SimpleNamespace(
            status="finished",
            prog="solve3D",
            master_stats={},
        )
        batch_job = Mock()
        batch_job.get_volume_outputs.return_value = [{
            "name": "recvol_state01.mrc",
            "path": volume_path,
        }]
        batch_job.get_stage_volume_outputs.return_value = []

        with patch.object(
            batch_views,
            "_get_accessible_batch_job",
            return_value=(batch_job, jobmodel),
        ):
            response = batch_views.view_batch_volume_data(
                self._get_request("/batchvolume/7/recvol_state01.mrc"),
                7,
                "recvol_state01.mrc",
            )
            version = volume_cache_version(volume_path)
            cached_response = batch_views.view_batch_volume_data(
                self._get_request(f"/batchvolume/7/recvol_state01.mrc?v={version}"),
                7,
                "recvol_state01.mrc",
            )
            source_stat = os.stat(volume_path)
            os.utime(volume_path, ns=(source_stat.st_atime_ns,
                                      source_stat.st_mtime_ns + 1_000_000_000))
            stale_response = batch_views.view_batch_volume_data(
                self._get_request(f"/batchvolume/7/recvol_state01.mrc?v={version}"),
                7,
                "recvol_state01.mrc",
            )

        self.assertEqual(response.status_code, 200)
        self.assertFalse(response.streaming)
        self.assertEqual(response["Content-Type"], "application/octet-stream")
        self.assertEqual(response["Content-Encoding"], "gzip")
        self.assertIn("no-store", response["Cache-Control"])
        self.assertIn("private", cached_response["Cache-Control"])
        self.assertIn("max-age=86400", cached_response["Cache-Control"])
        self.assertNotIn("no-store", cached_response["Cache-Control"])
        self.assertNotEqual(version, volume_cache_version(volume_path))
        self.assertIn("no-store", stale_response["Cache-Control"])
        wire_data = response.content
        decoded = gzip.decompress(wire_data)
        self.assertEqual(struct.unpack_from("<i", decoded, 12)[0], 0)
        converted = np.frombuffer(decoded, dtype=np.int8, offset=1024).reshape(4, 4, 4)
        self.assertEqual(converted[0, 0, 0], -128)
        # The retained source value 21 maps through the header range 0..63.
        self.assertEqual(converted[1, 1, 1], -43)
        self.assertEqual(struct.unpack_from("<f", decoded, 40)[0], 6.0)
        self.assertEqual(os.listdir(temporary_job_dir.name), ["recvol_state01.mrc"])

    def test_batch_volume_data_serves_running_refine3d_states_stage(self):
        temporary_job_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_job_dir.cleanup)
        volume_path = os.path.join(
            temporary_job_dir.name,
            "recvol_state01_stage01_lp.mrc",
        )
        with mrcfile.new(volume_path) as volume_file:
            volume_file.set_data(np.arange(64, dtype=np.float32).reshape(4, 4, 4))

        jobmodel = SimpleNamespace(
            status="running",
            prog="refine3D_states",
            master_stats={"project_metadata": {"cls3D": {"stage1": []}}},
        )
        batch_job = Mock()
        batch_job.get_volume_outputs.return_value = []
        batch_job.get_stage_volume_outputs.return_value = [{
            "name": "recvol_state01_stage01_lp.mrc",
            "path": volume_path,
        }]

        with patch.object(
            batch_views,
            "_get_accessible_batch_job",
            return_value=(batch_job, jobmodel),
        ):
            response = batch_views.view_batch_volume_data(
                self._get_request("/batchvolume/7/recvol_state01_stage01_lp.mrc"),
                7,
                "recvol_state01_stage01_lp.mrc",
            )
            jobmodel.status = "queued"
            queued_response = batch_views.view_batch_volume_data(
                self._get_request("/batchvolume/7/recvol_state01_stage01_lp.mrc"),
                7,
                "recvol_state01_stage01_lp.mrc",
            )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Encoding"], "gzip")
        self.assertEqual(struct.unpack_from("<i", gzip.decompress(
            response.content), 12)[0], 0)

        self.assertEqual(queued_response.status_code, 404)

    def test_batch_volume_data_rejects_an_undeclared_filename(self):
        jobmodel = SimpleNamespace(
            status="finished",
            prog="solve3D",
            master_stats={},
        )
        batch_job = Mock()
        batch_job.get_volume_outputs.return_value = [{
            "name": "recvol_state01.mrc",
            "path": "/workspace/9_solve3D/recvol_state01.mrc",
        }]
        batch_job.get_stage_volume_outputs.return_value = []

        with patch.object(
            batch_views,
            "_get_accessible_batch_job",
            return_value=(batch_job, jobmodel),
        ):
            response = batch_views.view_batch_volume_data(
                self._get_request("/batchvolume/7/other.mrc"),
                7,
                "other.mrc",
            )

        self.assertEqual(response.status_code, 404)

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

    def test_mark_finished_calls_batch_job_for_owned_queued_job(self):
        job = Mock()
        job.markComplete.return_value = True
        jobmodel = SimpleNamespace(status="queued")
        request = self.factory.post(
            "/finishbatch",
            {"selected_job_id": "7", "mark_finished": "1"},
        )
        request.user = _AuthUser()
        with (
            patch.object(
                batch_views,
                "_get_accessible_batch_job",
                return_value=(job, jobmodel),
            ),
            patch.object(batch_views.messages, "add_message"),
            patch.object(
                batch_views,
                "redirect",
                return_value=HttpResponseRedirect("/workspace"),
            ),
        ):
            response = batch_views.view_batch_mark_finished(
                request
            )

        self.assertEqual(response.status_code, 302)
        job.markComplete.assert_called_once_with(None, None)

    def test_mark_finished_rejects_other_job_statuses(self):
        for status in ("running", "finished", "stopped"):
            with self.subTest(status=status):
                job = Mock()
                jobmodel = SimpleNamespace(status=status)
                request = self.factory.post(
                    "/finishbatch",
                    {"selected_job_id": "7", "mark_finished": "1"},
                )
                request.user = _AuthUser()
                with (
                    patch.object(
                        batch_views,
                        "_get_accessible_batch_job",
                        return_value=(job, jobmodel),
                    ),
                    patch.object(batch_views.messages, "add_message"),
                    patch.object(
                        batch_views,
                        "redirect",
                        return_value=HttpResponseRedirect("/workspace"),
                    ),
                ):
                    response = batch_views.view_batch_mark_finished(
                        request
                    )

                self.assertEqual(response.status_code, 302)
                job.markComplete.assert_not_called()

    def test_mark_finished_requires_explicit_confirmation_key(self):
        job = Mock()
        jobmodel = SimpleNamespace(status="queued")
        with (
            patch.object(
                batch_views,
                "_get_accessible_batch_job",
                return_value=(job, jobmodel),
            ),
            patch.object(batch_views.messages, "add_message"),
            patch.object(
                batch_views,
                "redirect",
                return_value=HttpResponseRedirect("/workspace"),
            ),
        ):
            response = batch_views.view_batch_mark_finished(
                self._request("/finishbatch")
            )

        self.assertEqual(response.status_code, 302)
        job.markComplete.assert_not_called()

    def test_mark_finished_requires_post(self):
        response = batch_views.view_batch_mark_finished(
            self._get_request("/finishbatch")
        )

        self.assertEqual(response.status_code, 405)

    def test_rerun_requires_post(self):
        response = batch_views.view_batch_rerun(
            self._get_request("/rerunbatch")
        )

        self.assertEqual(response.status_code, 405)

    def test_rerun_opens_job_builder_for_owned_terminal_batch_job(self):
        job = Mock()
        metadata = {
            "source": {"type": "project_file", "filename": "input.simple"},
        }
        jobmodel = SimpleNamespace(
            id=7,
            status="finished",
            pckg="simple",
            prog="refine2D",
            master_stats=metadata,
            dset_id=3,
        )

        with patch.object(
            batch_views,
            "_get_accessible_batch_job",
            return_value=(job, jobmodel),
        ):
            response = batch_views.view_batch_rerun(
                self._request("/rerunbatch")
            )

        self.assertEqual(response.status_code, 302)
        self.assertEqual(response.url, "/newstream?selected_job_id=7")
        job.rerun.assert_not_called()

    def test_rerun_rejects_active_or_unsupported_batch_jobs(self):
        for status, package in (
            ("running", "simple"),
            ("finished", "simple_stream"),
        ):
            with self.subTest(status=status, package=package):
                job = Mock()
                jobmodel = SimpleNamespace(
                    status=status,
                    pckg=package,
                    prog="demo",
                    master_stats={
                    },
                )
                with (
                    patch.object(
                        batch_views,
                        "_get_accessible_batch_job",
                        return_value=(job, jobmodel),
                    ),
                    patch.object(batch_views.messages, "add_message"),
                    patch.object(
                        batch_views,
                        "redirect",
                        return_value=HttpResponseRedirect("/workspace"),
                    ),
                ):
                    response = batch_views.view_batch_rerun(
                        self._request("/rerunbatch")
                    )

                self.assertEqual(response.status_code, 302)
                job.rerun.assert_not_called()

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
