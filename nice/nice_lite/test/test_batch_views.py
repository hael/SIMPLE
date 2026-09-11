import os
import tempfile
from types import SimpleNamespace
from unittest.mock import Mock, patch
from urllib.parse import parse_qs, urlparse

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
            pckg="simple",
            prog="import_movies",
            master_stats={
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
            patch.object(
                batch_views,
                "_project_sampling_distance",
                return_value=0.885,
            ),
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
        self.assertEqual(context["project_sampling_distance"], 0.885)
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
        }])

    def test_submitted_mask_diameter_accepts_only_positive_finite_numbers(self):
        for value, expected in (
            ("190", 190),
            (120.6, 120.6),
            (0.4, 0.4),
            (0, None),
            (-10, None),
            (True, None),
            ("nan", None),
            ("invalid", None),
            (None, None),
        ):
            with self.subTest(value=value):
                jobmodel = SimpleNamespace(args={"mskdiam": value})
                self.assertEqual(
                    batch_views._submitted_mask_diameter(jobmodel),
                    expected,
                )

    def test_project_sampling_distance_prefers_output_then_stack_and_micrograph(self):
        reader = Mock()
        reader.read_records.side_effect = lambda oritype: {
            "out": [{"imgkind": "cavg", "smpd": "1.5"}],
            "stk": [{"smpd": 1.3}],
            "mic": [{"smpd": 0.885}],
        }[oritype]

        with patch.object(
            batch_views,
            "SIMPLEProjectFileReader",
            return_value=reader,
        ):
            sampling_distance = batch_views._project_sampling_distance(
                "/project/workspace.simple"
            )

        self.assertEqual(sampling_distance, 1.5)
        reader.read_records.assert_called_once_with("out")

    def test_project_sampling_distance_falls_back_to_stack_records(self):
        reader = Mock()
        reader.read_records.side_effect = lambda oritype: {
            "out": [{"imgkind": "frc2D"}],
            "stk": [{"smpd": 1.3}],
            "mic": [{"smpd": 0.885}],
        }[oritype]

        with patch.object(
            batch_views,
            "SIMPLEProjectFileReader",
            return_value=reader,
        ):
            sampling_distance = batch_views._project_sampling_distance(
                "/project/workspace.simple"
            )

        self.assertEqual(sampling_distance, 1.3)
        self.assertEqual(
            [call.args[0] for call in reader.read_records.call_args_list],
            ["out", "stk"],
        )

    def test_project_picked_particle_count_sums_micrograph_records(self):
        reader = Mock()
        reader.read_records.return_value = [
            {"nptcls": 229.0},
            {"nptcls": 0.0},
            {"nptcls": 241},
        ]

        with patch.object(
            batch_views,
            "SIMPLEProjectFileReader",
            return_value=reader,
        ):
            particle_count = batch_views._project_picked_particle_count(
                "/project/workspace.simple"
            )

        self.assertEqual(particle_count, 470)
        reader.read_records.assert_called_once_with("mic")

    def test_project_picked_particle_count_rejects_incomplete_counts(self):
        reader = Mock()
        reader.read_records.return_value = [
            {"nptcls": 229.0},
            {},
        ]

        with patch.object(
            batch_views,
            "SIMPLEProjectFileReader",
            return_value=reader,
        ):
            particle_count = batch_views._project_picked_particle_count(
                "/project/workspace.simple"
            )

        self.assertIsNone(particle_count)

    def test_class_overlay_uses_pixels_until_mask_and_sampling_are_available(self):
        jobmodel = SimpleNamespace(args={"mskdiam": "190"})
        no_sampling = batch_views._class_overlay_settings(
            jobmodel,
            SimpleNamespace(sampling_distance=None),
        )
        no_mask = batch_views._class_overlay_settings(
            SimpleNamespace(args={}),
            SimpleNamespace(sampling_distance=1.3),
        )
        project_sampling = batch_views._class_overlay_settings(
            jobmodel,
            SimpleNamespace(sampling_distance=None),
            project_sampling_distance=1.3,
        )

        self.assertEqual(no_sampling["overlay_unit"], "pixels")
        self.assertIsNone(no_sampling["overlay_size"])
        self.assertEqual(no_mask["overlay_unit"], "pixels")
        self.assertIsNone(no_mask["overlay_size"])
        self.assertEqual(no_mask["sampling_distance"], 1.3)
        self.assertEqual(project_sampling["overlay_unit"], "angstroms")
        self.assertAlmostEqual(
            project_sampling["overlay_size"],
            190 / 1.3,
            places=5,
        )

    def test_batch_class_selector_is_loaded_only_by_explicit_query_key(self):
        project = SimpleNamespace(name="project", dirc="/project")
        workspace = SimpleNamespace(name="workspace", proj=project)
        jobmodel = SimpleNamespace(
            id=7,
            disp=2,
            name="2D classification",
            desc="",
            status="finished",
            cdat="created",
            args={"mskdiam": "190"},
            pckg="simple",
            prog="abinitio2D",
            dset=workspace,
        )
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = "/project/workspace.simple"
        batch_job.get_artifact_summary.return_value = {"counts": [], "images": []}
        batch_job.get_safe_job_dir.return_value = "/project/2_cluster2D"
        batch_job.get_log_tails.return_value = []
        selection = SimpleNamespace(
            classes=({"class_id": 1, "stack_index": 1},),
            stack_name="classes.mrcs",
            width=128,
            height=128,
            sampling_distance=1.3,
            initial_selected_class_ids=(1,),
            browser_data=lambda: {
                "classes": [{"class_id": 1, "stack_index": 1}],
                "initial_selected_class_ids": [1],
                "selection_storage_key": "selection-key",
            },
        )

        with (
            patch.object(batch_views, "SIMPLEProjFile") as project_reader_cls,
            patch.object(
                batch_views,
                "load_batch_class_selection",
                return_value=selection,
            ) as load_selection,
        ):
            project_reader_cls.return_value.getGlobalStats.return_value = {
                "cls2D": {"n": 1}
            }
            closed_context = batch_views._batch_detail_context(batch_job, jobmodel)
            open_context = batch_views._batch_detail_context(
                batch_job,
                jobmodel,
                class_selector_requested=True,
            )

        self.assertTrue(closed_context["class_selector_available"])
        self.assertIsNone(closed_context["batch_class_selector"])
        self.assertFalse(
            closed_context["class_selector_replaces_artifact_previews"]
        )
        load_selection.assert_called_once_with(
            "/project/workspace.simple",
            "/project",
            7,
        )
        self.assertEqual(open_context["batch_class_selector"]["class_count"], 1)
        self.assertEqual(open_context["batch_class_selector"]["width"], 128)
        self.assertEqual(open_context["batch_class_selector"]["height"], 128)
        self.assertAlmostEqual(
            open_context["batch_class_selector"]["sampling_distance"],
            1.3,
            places=6,
        )
        self.assertAlmostEqual(
            open_context["batch_class_selector"]["overlay_size"],
            190 / 1.3,
            places=5,
        )
        self.assertEqual(
            open_context["batch_class_selector"]["overlay_display_size"],
            190,
        )
        self.assertEqual(
            open_context["batch_class_selector"]["overlay_unit"],
            "angstroms",
        )
        self.assertTrue(
            open_context["class_selector_replaces_artifact_previews"]
        )
        self.assertTrue(open_context["output_dimensions_available"])

    def test_batch_class_thumbnail_renders_owned_stack_slice(self):
        jobmodel = SimpleNamespace(
            id=7,
            status="finished",
            dset=SimpleNamespace(proj=SimpleNamespace(dirc="/project")),
        )
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = "/project/workspace.simple"
        selection = SimpleNamespace(
            stack_path="/project/classes.mrcs",
            classes=({"class_id": 4, "stack_index": 2},),
        )

        with (
            patch.object(
                batch_views,
                "_get_accessible_batch_job",
                return_value=(batch_job, jobmodel),
            ),
            patch.object(
                batch_views,
                "load_batch_class_selection",
                return_value=selection,
            ),
            patch.object(
                batch_views,
                "render_mrc_particle_png",
                return_value=b"png",
            ) as render_thumbnail,
        ):
            response = batch_views.view_batch_class_thumbnail(
                self._get_request("/batchclass/7/2"),
                7,
                2,
            )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "image/png")
        render_thumbnail.assert_called_once_with(
            "/project/classes.mrcs",
            2,
            max_size=512,
        )

    def test_abinitio3d_volume_viewer_is_explicit_and_omits_server_path(self):
        jobmodel = SimpleNamespace(
            id=7,
            disp=9,
            name="Initial 3D Reconstruction",
            desc="",
            status="finished",
            cdat="created",
            args={},
            pckg="simple",
            prog="abinitio3D",
            master_stats={
            },
            dset=SimpleNamespace(
                name="workspace",
                proj=SimpleNamespace(name="project", dirc="/workspace"),
            ),
        )
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = "/workspace/9_abinitio3D/workspace.simple"
        batch_job.get_artifact_summary.return_value = {"counts": [], "images": []}
        batch_job.get_safe_job_dir.return_value = "/workspace/9_abinitio3D"
        batch_job.get_log_tails.return_value = []
        batch_job.get_volume_outputs.return_value = [{
            "path": "/workspace/9_abinitio3D/recvol_state01.mrc",
            "name": "recvol_state01.mrc",
            "state": 1,
            "population": 5542,
            "width": 256,
            "height": 256,
            "depth": 256,
            "voxel_size": (1.3, 1.3, 1.3),
            "minimum": -2.0,
            "maximum": 8.0,
        }]

        with patch.object(batch_views, "SIMPLEProjFile") as projfile:
            projfile.return_value.getGlobalStats.return_value = {}
            closed_context = batch_views._batch_detail_context(batch_job, jobmodel)
            open_context = batch_views._batch_detail_context(
                batch_job,
                jobmodel,
                volume_viewer_requested=True,
            )

        self.assertTrue(closed_context["volume_viewer_available"])
        self.assertFalse(closed_context["volume_viewer_requested"])
        self.assertEqual(closed_context["batch_volume_viewer"], [])
        self.assertTrue(open_context["volume_viewer_requested"])
        self.assertEqual(
            open_context["batch_volume_viewer"][0]["name"],
            "recvol_state01.mrc",
        )
        self.assertNotIn("path", open_context["batch_volume_viewer"][0])

    def test_batch_volume_data_streams_the_owned_mrc_file(self):
        temporary_job_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_job_dir.cleanup)
        volume_path = os.path.join(
            temporary_job_dir.name,
            "recvol_state01.mrc",
        )
        volume_bytes = b"MRC volume bytes"
        with open(volume_path, "wb") as volume_file:
            volume_file.write(volume_bytes)

        jobmodel = SimpleNamespace(
            status="finished",
            prog="abinitio3D",
        )
        batch_job = Mock()
        batch_job.get_volume_outputs.return_value = [{
            "name": "recvol_state01.mrc",
            "path": volume_path,
        }]

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

        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.streaming)
        self.assertEqual(response["Content-Type"], "application/octet-stream")
        self.assertEqual(response["Content-Length"], str(len(volume_bytes)))
        self.assertIn(
            "recvol_state01.mrc",
            response["Content-Disposition"],
        )
        self.assertEqual(b"".join(response.streaming_content), volume_bytes)
        response.close()

    def test_batch_volume_data_rejects_an_undeclared_filename(self):
        jobmodel = SimpleNamespace(
            status="finished",
            prog="abinitio3D",
        )
        batch_job = Mock()
        batch_job.get_volume_outputs.return_value = [{
            "name": "recvol_state01.mrc",
            "path": "/workspace/9_abinitio3D/recvol_state01.mrc",
        }]

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

    def test_batch_class_export_returns_project_ordered_state_infile(self):
        jobmodel = SimpleNamespace(
            id=7,
            status="finished",
            dset=SimpleNamespace(proj=SimpleNamespace(dirc="/project")),
        )
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = "/project/workspace.simple"
        selection = SimpleNamespace(classes=(
            {"class_id": 3},
            {"class_id": 2},
            {"class_id": 1},
        ))
        request = self.factory.post(
            "/batchclass/7/infile",
            {"selected_class_ids": "[3, 1]"},
        )
        request.user = _AuthUser()

        with (
            patch.object(
                batch_views,
                "_get_accessible_batch_job",
                return_value=(batch_job, jobmodel),
            ),
            patch.object(
                batch_views,
                "load_batch_class_selection",
                return_value=selection,
            ),
        ):
            response = batch_views.view_batch_class_selection_export(
                request,
                7,
            )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content, b"1\n0\n1\n")
        self.assertIn("batch_7_class_selection.txt", response["Content-Disposition"])

    def test_batch_class_selection_replaces_infile_and_opens_prefilled_builder(self):
        temporary_job_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_job_dir.cleanup)
        infile_path = os.path.join(temporary_job_dir.name, "class_selection.txt")
        with open(infile_path, "w", encoding="utf-8") as infile:
            infile.write("stale\n")

        project = SimpleNamespace(id=3, dirc="/project")
        jobmodel = SimpleNamespace(
            id=7,
            disp=4,
            name="Create 2D Class Averages",
            status="finished",
            pckg="simple",
            prog="abinitio2D",
            dset=SimpleNamespace(proj=project, proj_id=project.id),
            dset_id=9,
        )
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = (
            "/project/.workspace_9/4_abinitio2D/workspace.simple"
        )
        batch_job.get_safe_job_dir.return_value = temporary_job_dir.name
        selection = SimpleNamespace(classes=(
            {"class_id": 3},
            {"class_id": 1},
            {"class_id": 2},
        ))
        request = self.factory.post(
            "/batchclass/7/selection",
            {"selected_class_ids": "[2, 3]"},
        )
        request.user = _AuthUser()

        with (
            patch.object(
                batch_views,
                "_get_accessible_batch_job",
                return_value=(batch_job, jobmodel),
            ),
            patch.object(
                batch_views,
                "load_batch_class_selection",
                return_value=selection,
            ),
            patch.object(batch_views.messages, "add_message") as add_message,
        ):
            response = batch_views.view_batch_class_selection_run(request, 7)

        self.assertEqual(response.status_code, 302)
        self.assertEqual(urlparse(response.url).path, "/workspace")
        with open(infile_path, encoding="utf-8") as infile:
            self.assertEqual(infile.read(), "1\n0\n1\n")
        query = parse_qs(urlparse(response.url).query)
        self.assertEqual(
            query,
            {
                "selected_job_id": ["7"],
                "class_selection": ["1"],
                "selected_project_id": ["3"],
                "selected_workspace_id": ["9"],
            },
        )
        self.assertIn(
            "class selection infile saved; review and start the selection job",
            add_message.call_args.args[2],
        )

    def test_batch_class_selection_rejects_empty_selection(self):
        jobmodel = SimpleNamespace(
            id=7,
            status="finished",
            pckg="simple",
            prog="abinitio2D",
            dset=SimpleNamespace(proj=SimpleNamespace(dirc="/project")),
        )
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = "/project/workspace.simple"
        selection = SimpleNamespace(classes=({"class_id": 1},))
        request = self.factory.post(
            "/batchclass/7/selection",
            {"selected_class_ids": "[]"},
        )
        request.user = _AuthUser()

        with (
            patch.object(
                batch_views,
                "_get_accessible_batch_job",
                return_value=(batch_job, jobmodel),
            ),
            patch.object(
                batch_views,
                "load_batch_class_selection",
                return_value=selection,
            ),
            patch.object(
                batch_views,
                "_save_batch_class_selection_infile",
            ) as save_infile,
            patch.object(batch_views.messages, "add_message") as add_message,
        ):
            response = batch_views.view_batch_class_selection_run(request, 7)

        self.assertEqual(response.status_code, 302)
        self.assertIn("?class_selector=1#batch_class_selector", response.url)
        save_infile.assert_not_called()
        self.assertIn(
            "Select at least one class before running a selection job.",
            add_message.call_args.args[2],
        )

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
            pckg="unknown",
            prog="extract",
            dset=workspace,
        )
        particle_stack_page = {
            "stacks": [{"name": "particles.mrcs", "count": 80}],
            "particles": [{
                "number": 41,
                "stack_name": "particles.mrcs",
                "stack_index": 41,
                "width": 128,
                "height": 128,
            }],
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
        self.assertTrue(context["output_dimensions_available"])

    def test_reproject_batch_context_paginates_projection_stack_headers(self):
        project = SimpleNamespace(name="project")
        workspace = SimpleNamespace(name="workspace", proj=project)
        jobmodel = SimpleNamespace(
            id=7,
            disp=1,
            name="Reproject Volume",
            desc="",
            status="finished",
            cdat="created",
            args={},
            pckg="simple",
            prog="reproject",
            dset=workspace,
        )
        stack_page = {
            "stacks": [{"name": "reprojs.mrcs", "count": 60}],
            "particles": [{
                "number": 1,
                "stack_name": "reprojs.mrcs",
                "stack_index": 1,
                "width": 240,
                "height": 240,
            }],
            "total": 60,
        }
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = None
        batch_job.get_artifact_summary.return_value = {"counts": [], "images": []}
        batch_job.get_particle_stack_page.return_value = stack_page
        batch_job.get_safe_job_dir.return_value = "/project/workspace/1_reproject"
        batch_job.get_log_tails.return_value = []

        context = batch_views._batch_detail_context(batch_job, jobmodel)

        batch_job.get_particle_stack_page.assert_called_once_with(page=1, page_size=40)
        self.assertIs(context["particle_stack_page"], stack_page)
        self.assertEqual(context["artifact_counts"], [{"extension": "MRCS", "count": 1}])
        self.assertTrue(context["output_dimensions_available"])

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

        with patch.object(
            batch_views,
            "read_movie_dimensions",
            return_value=(4096, 4096),
        ) as read_dimensions:
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
        self.assertEqual(movie_page["movies"][0]["width"], 4096)
        self.assertEqual(movie_page["movies"][0]["height"], 4096)
        self.assertEqual(read_dimensions.call_count, 40)
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
            pckg="unknown",
            prog="ctf_estimate",
            master_stats={
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
        self.assertFalse(context["output_dimensions_available"])

    def test_batch_detail_uses_picker_generated_denoised_thumbnails(self):
        jobmodel = SimpleNamespace(
            id=7,
            disp=3,
            name="Pick Particles from Micrographs",
            desc="",
            status="finished",
            cdat="created",
            args={},
            pckg="unknown",
            prog="pick",
            master_stats={
            },
            dset=SimpleNamespace(
                name="workspace",
                proj=SimpleNamespace(name="project"),
            ),
        )
        batch_job = Mock()
        batch_job.get_result_project_path.return_value = "/workspace/3_pick/workspace.simple"
        batch_job.get_artifact_summary.return_value = {
            "counts": [{"extension": "JPEG", "count": 1}],
            "images": [{
                "name": "pickrefs_source.jpeg",
                "path": "/workspace/3_pick/pickrefs_source.jpeg",
                "previews": [],
            }],
        }
        batch_job.get_safe_job_dir.return_value = "/workspace/3_pick"
        batch_job.get_log_tails.return_value = []
        batch_job.get_pick_micrograph_previews.return_value = [{
            "path": "/workspace/3_pick/movie_intg_den.jpg",
            "number": 11,
            "xdim": 4096,
            "ydim": 3072,
            "width": 512,
            "height": 512,
            "boxes": [{"x": 101, "y": 202, "width": 180, "height": 180}],
        }]

        with (
            patch.object(batch_views, "SIMPLEProjFile") as projfile,
            patch.object(
                batch_views,
                "_project_picked_particle_count",
                return_value=5708,
            ) as picked_particle_count,
        ):
            project_reader = projfile.return_value
            project_reader.getGlobalStats.return_value = {"mic": {"n": 30}}

            context = batch_views._batch_detail_context(batch_job, jobmodel)

        batch_job.get_pick_micrograph_previews.assert_called_once_with(
            max_previews=20,
            max_coordinates=1500,
        )
        picked_particle_count.assert_called_once_with(
            "/workspace/3_pick/workspace.simple"
        )
        project_reader.getFieldStats.assert_not_called()
        self.assertEqual(context["pick_micrographs"], [{
            "path": "/workspace/3_pick/movie_intg_den.jpg",
            "number": 11,
            "xdim": 4096,
            "ydim": 3072,
            "width": 512,
            "height": 512,
            "boxes": [{"x": 101, "y": 202, "width": 180, "height": 180}],
        }])
        self.assertEqual(context["artifact_images"], [])
        self.assertEqual(
            context["artifact_counts"],
            [{"extension": "JPEG", "count": 1}],
        )
        self.assertEqual(context["pick_particle_count"], 5708)
        self.assertTrue(context["pick_box_overlay_available"])
        self.assertTrue(context["output_dimensions_available"])

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
            prog="cluster2D",
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
