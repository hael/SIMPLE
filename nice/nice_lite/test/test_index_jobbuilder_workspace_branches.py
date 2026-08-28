from types import SimpleNamespace
from unittest.mock import Mock
from unittest.mock import patch

from django.http import HttpResponse
from django.template.loader import render_to_string
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

    def test_stream_empty_project_selector_keeps_create_new_selectable(self):
        html = render_to_string(
            "index.html",
            {
                "current_project_id": None,
                "current_workspace_id": None,
                "projects": [],
                "workspaces": [],
                "iframeurl": None,
            },
        )

        select_index = html.index('<option disabled selected value="">select</option>')
        create_index = html.index('<option value="-1">')
        self.assertLess(select_index, create_index)

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

    def test_collect_programs_preserves_cross_field_requirements(self):
        requirements = [{
            "label": "Project location",
            "help": "Supply a project name, a project directory, or both.",
            "min_selected": 1,
            "max_selected": 2,
            "keys": ["projname", "dir"],
        }]
        batchui = {
            "new_project": {
                "program": {
                    "executable": "simple_exec",
                    "requirements": requirements,
                },
                "inputs": [{"key": "projname"}, {"key": "dir"}],
            },
        }

        _, program_inputs = job_builder_views._collect_programs(batchui, "simple_exec")

        self.assertEqual(program_inputs[0]["requirements"], requirements)

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

    def test_create_batch_allowlists_program_and_arguments_from_ui(self):
        request = self.factory.post("/createbatch", {
            "package": "single",
            "program": "pick",
            "mode": "fast",
            "unknown": "discard-me",
        })
        request.user = _AuthUser()

        workspace = Mock()
        launcher = Mock()
        launcher.get_ui.return_value = {
            "pick": {
                "program": {"executable": "single_exec"},
                "inputs": [
                    {"key": "mode", "options": ["fast", "slow"]},
                    {"key": "projfile", "required": True},
                ],
            },
        }
        batchjob = Mock()
        batchjob.new.return_value = True

        with patch.object(job_builder_views, "get_workspace_id", return_value=4), patch.object(job_builder_views, "get_project_id", return_value=3), patch.object(job_builder_views, "Workspace", return_value=workspace), patch.object(job_builder_views, "_is_workspace_accessible", return_value=True), patch.object(job_builder_views, "SIMPLEBatch", return_value=launcher), patch.object(job_builder_views, "BatchJob", return_value=batchjob), patch.object(job_builder_views.messages, "add_message"):
            response = job_builder_views.view_create_batch(request)

        self.assertEqual(response.status_code, 302)
        batchjob.new.assert_called_once_with(workspace, "single", "pick", {"mode": "fast"})

    def test_collect_batch_args_accepts_free_form_directory_with_empty_options(self):
        program_cfg = {
            "program": {"executable": "simple_exec"},
            "files": [{
                "key": "dir_movies",
                "keytype": "dir",
                "options": [],
            }],
        }

        args, error = job_builder_views._collect_batch_args(
            {"dir_movies": "/data/movies"},
            program_cfg,
        )

        self.assertIsNone(error)
        self.assertEqual(args, {"dir_movies": "/data/movies"})

    def test_create_batch_rejects_program_for_wrong_executable(self):
        request = self.factory.post("/createbatch", {"package": "single", "program": "pick"})
        request.user = _AuthUser()
        launcher = Mock()
        launcher.get_ui.return_value = {
            "pick": {"program": {"executable": "simple_exec"}},
        }

        with patch.object(job_builder_views, "get_workspace_id", return_value=4), patch.object(job_builder_views, "get_project_id", return_value=3), patch.object(job_builder_views, "Workspace"), patch.object(job_builder_views, "_is_workspace_accessible", return_value=True), patch.object(job_builder_views, "SIMPLEBatch", return_value=launcher), patch.object(job_builder_views, "BatchJob") as batchjob_cls, patch.object(job_builder_views.messages, "add_message"):
            response = job_builder_views.view_create_batch(request)

        self.assertEqual(response.status_code, 302)
        batchjob_cls.assert_not_called()

    def test_create_batch_applies_nthr_default_but_omits_other_untouched_defaults(self):
        request = self.factory.post("/createbatch", {
            "package": "simple",
            "program": "cluster2D",
        })
        request.user = _AuthUser()

        workspace = Mock()
        launcher = Mock()
        launcher.get_ui.return_value = {
            "cluster2D": {
                "program": {"executable": "simple_exec"},
                "compute": [
                    {"key": "nthr", "keytype": "int", "has_default": True, "default": 8.0},
                    {"key": "scale", "keytype": "float", "has_default": True, "default": 1.5},
                ],
            },
        }
        batchjob = Mock()
        batchjob.new.return_value = True

        with patch.object(job_builder_views, "get_workspace_id", return_value=4), patch.object(job_builder_views, "get_project_id", return_value=3), patch.object(job_builder_views, "Workspace", return_value=workspace), patch.object(job_builder_views, "_is_workspace_accessible", return_value=True), patch.object(job_builder_views, "SIMPLEBatch", return_value=launcher), patch.object(job_builder_views, "BatchJob", return_value=batchjob), patch.object(job_builder_views.messages, "add_message"):
            response = job_builder_views.view_create_batch(request)

        self.assertEqual(response.status_code, 302)
        batchjob.new.assert_called_once_with(
            workspace,
            "simple",
            "cluster2D",
            {"nthr": "8"},
        )

    def test_batch_requirements_reject_missing_and_excess_alternatives(self):
        program_cfg = {
            "program": {
                "requirements": [{
                    "label": "Input data",
                    "help": "Supply either a volume or an image stack.",
                    "min_selected": 1,
                    "max_selected": 1,
                    "keys": ["vol1", "stk"],
                }],
            },
        }

        self.assertEqual(
            job_builder_views._validate_batch_requirements(program_cfg, {}),
            "Supply either a volume or an image stack.",
        )
        self.assertIsNone(job_builder_views._validate_batch_requirements(
            program_cfg,
            {"vol1": "map.mrc"},
        ))
        self.assertEqual(
            job_builder_views._validate_batch_requirements(
                program_cfg,
                {"vol1": "map.mrc", "stk": "particles.mrcs"},
            ),
            "Supply either a volume or an image stack.",
        )

    def test_batch_requirements_include_launcher_project_file(self):
        program_cfg = {
            "program": {
                "requirements": [{
                    "label": "Input data",
                    "help": "Supply an input source.",
                    "min_selected": 1,
                    "max_selected": 1,
                    "keys": ["vol1", "projfile"],
                }],
            },
        }

        self.assertIsNone(job_builder_views._validate_batch_requirements(program_cfg, {}))

    def test_create_batch_rejects_unsatisfied_ui_requirement_before_launch(self):
        request = self.factory.post("/createbatch", {
            "package": "simple",
            "program": "new_project",
        })
        request.user = _AuthUser()
        launcher = Mock()
        launcher.get_ui.return_value = {
            "new_project": {
                "program": {
                    "executable": "simple_exec",
                    "requirements": [{
                        "label": "Project location",
                        "help": "Supply a project name, a project directory, or both.",
                        "min_selected": 1,
                        "max_selected": 2,
                        "keys": ["projname", "dir"],
                    }],
                },
                "inputs": [
                    {"key": "projname", "keytype": "str"},
                    {"key": "dir", "keytype": "dir"},
                ],
            },
        }

        with patch.object(job_builder_views, "get_workspace_id", return_value=4), patch.object(job_builder_views, "get_project_id", return_value=3), patch.object(job_builder_views, "Workspace"), patch.object(job_builder_views, "_is_workspace_accessible", return_value=True), patch.object(job_builder_views, "SIMPLEBatch", return_value=launcher), patch.object(job_builder_views, "BatchJob") as batchjob_cls, patch.object(job_builder_views.messages, "add_message"):
            response = job_builder_views.view_create_batch(request)

        self.assertEqual(response.status_code, 302)
        batchjob_cls.assert_not_called()

    def test_import_movies_does_not_materialize_unselected_directory_default(self):
        request = self.factory.post("/createbatch", {
            "package": "simple",
            "program": "import_movies",
            "filetab": "/data/movies.txt",
        })
        request.user = _AuthUser()
        workspace = Mock()
        launcher = Mock()
        launcher.get_ui.return_value = {
            "import_movies": {
                "program": {"executable": "simple_exec"},
                "files": [
                    {
                        "key": "filetab",
                        "keytype": "file",
                        "has_default": True,
                        "default": "",
                    },
                    {
                        "key": "dir_movies",
                        "keytype": "dir",
                        "has_default": True,
                        "default": "preprocess/",
                    },
                ],
            },
        }
        batchjob = Mock()
        batchjob.new.return_value = True

        with (
            patch.object(job_builder_views, "get_workspace_id", return_value=4),
            patch.object(job_builder_views, "get_project_id", return_value=3),
            patch.object(job_builder_views, "Workspace", return_value=workspace),
            patch.object(job_builder_views, "_is_workspace_accessible", return_value=True),
            patch.object(job_builder_views, "SIMPLEBatch", return_value=launcher),
            patch.object(job_builder_views, "BatchJob", return_value=batchjob),
            patch.object(job_builder_views.messages, "add_message"),
        ):
            response = job_builder_views.view_create_batch(request)

        self.assertEqual(response.status_code, 302)
        batchjob.new.assert_called_once_with(
            workspace,
            "simple",
            "import_movies",
            {"filetab": "/data/movies.txt"},
        )

    def test_import_movies_requires_exactly_one_movie_source(self):
        self.assertEqual(
            job_builder_views._validate_batch_program_args("import_movies", {}),
            "choose a movie-file list or an input movies directory",
        )
        self.assertEqual(
            job_builder_views._validate_batch_program_args(
                "import_movies",
                {"filetab": "movies.txt", "dir_movies": "/data/movies"},
            ),
            "choose either a movie-file list or an input movies directory, not both",
        )
        self.assertIsNone(job_builder_views._validate_batch_program_args(
            "import_movies",
            {"dir_movies": "/data/movies"},
        ))

    def test_create_batch_rejects_nonfinite_numeric_input(self):
        request = self.factory.post("/createbatch", {
            "package": "simple",
            "program": "cluster2D",
            "nthr": "nan",
        })
        request.user = _AuthUser()
        launcher = Mock()
        launcher.get_ui.return_value = {
            "cluster2D": {
                "program": {"executable": "simple_exec"},
                "compute": [{"key": "nthr", "keytype": "int", "required": True}],
            },
        }

        with patch.object(job_builder_views, "get_workspace_id", return_value=4), patch.object(job_builder_views, "get_project_id", return_value=3), patch.object(job_builder_views, "Workspace"), patch.object(job_builder_views, "_is_workspace_accessible", return_value=True), patch.object(job_builder_views, "SIMPLEBatch", return_value=launcher), patch.object(job_builder_views, "BatchJob") as batchjob_cls, patch.object(job_builder_views.messages, "add_message"):
            response = job_builder_views.view_create_batch(request)

        self.assertEqual(response.status_code, 302)
        batchjob_cls.assert_not_called()

    def test_create_batch_passes_resolved_snapshot_to_batch_job(self):
        request = self.factory.post("/createbatch", {
            "package": "simple",
            "program": "cluster2D",
            "batch_source": "snapshot:12:3",
        })
        request.user = _AuthUser()
        workspace = Mock()
        launcher = Mock()
        launcher.get_ui.return_value = {
            "cluster2D": {
                "program": {"executable": "simple_exec"},
                "compute": [],
            },
        }
        batchjob = Mock()
        batchjob.new.return_value = True
        source = {
            "type": "stream_snapshot",
            "stream_job_id": 12,
            "particle_set_id": 3,
            "filename": "snapshot_3.simple",
        }

        with (
            patch.object(job_builder_views, "get_workspace_id", return_value=4),
            patch.object(job_builder_views, "get_project_id", return_value=3),
            patch.object(job_builder_views, "Workspace", return_value=workspace),
            patch.object(job_builder_views, "_is_workspace_accessible", return_value=True),
            patch.object(job_builder_views, "SIMPLEBatch", return_value=launcher),
            patch.object(
                job_builder_views,
                "_resolve_batch_project_source",
                return_value=("/workspace/snapshot_3.simple", source, None),
            ),
            patch.object(job_builder_views, "BatchJob", return_value=batchjob),
            patch.object(job_builder_views.messages, "add_message"),
        ):
            response = job_builder_views.view_create_batch(request)

        self.assertEqual(response.status_code, 302)
        batchjob.new.assert_called_once_with(
            workspace,
            "simple",
            "cluster2D",
            {},
            parent_proj="/workspace/snapshot_3.simple",
            source=source,
        )


class WorkspaceAccessBranchTests(SimpleTestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def test_view_workspace_inaccessible_returns_204(self):
        request = self.factory.get("/workspace")
        request.user = _AuthUser()

        with patch.object(workspace_views, "get_workspace_id", return_value=1), patch.object(workspace_views, "get_project_id", return_value=2), patch.object(workspace_views, "Workspace", return_value=SimpleNamespace()), patch.object(workspace_views, "_is_workspace_accessible", return_value=False), patch.object(workspace_views.messages, "add_message"), patch.object(workspace_views, "print_error"):
            response = workspace_views.view_workspace(request)

        self.assertEqual(response.status_code, 204)
