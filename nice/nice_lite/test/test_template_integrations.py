from pathlib import Path

from django.test import SimpleTestCase
from django.template.loader import render_to_string


class TemplateIntegrationTests(SimpleTestCase):
    @staticmethod
    def _read_template(relative_path):
        base_dir = Path(__file__).resolve().parents[1]
        template_path = base_dir / "templates" / relative_path
        return template_path.read_text(encoding="utf-8")

    def test_shared_2d_class_view_sorts_final_stage_by_resolution_on_load(self):
        viewer = self._read_template("includes/_cls2D_viewer.html")

        self.assertIn("const sortClsByResolution = (scope) => {", viewer)
        self.assertIn(
            'templates.sort((a, b) => '
            'getNumber(a, "res") - getNumber(b, "res"));',
            viewer,
        )
        self.assertIn("sortClsByResolution(getActiveClsScope());", viewer)
        self.assertIn("const finalStage = stageGroups.at(-1);", viewer)
        self.assertIn("sortClsByResolution(finalStage || slider);", viewer)
        self.assertNotIn("stageGroups.forEach(sortClsByResolution);", viewer)

    def test_zoom_log_template_uses_dataset_not_workspace_property(self):
        content = self._read_template("nice_stream/includes/_zoom_log_errors_parts_body.html")
        self.assertNotIn("panel.workspace", content)
        self.assertIn("panel.dataset.expanded", content)

    def test_workspace_reload_guard_skips_when_job_builder_visible(self):
        content = self._read_template("workspace.html")
        self.assertIn("isJobBuilderVisible", content)
        self.assertIn("&& !isJobBuilderVisible()", content)

    def test_job_builder_button_closes_visible_builder_and_opens_hidden_builder(self):
        workspace = self._read_template("workspace.html")

        self.assertIn('onsubmit="return toggleJobBuilder(event)"', workspace)
        self.assertIn('title="toggle job builder" aria-label="toggle job builder"', workspace)
        self.assertIn("const toggleJobBuilder = (event) => {", workspace)
        self.assertIn("if (!isJobBuilderVisible()) return true;", workspace)
        self.assertIn("event.preventDefault();", workspace)
        self.assertIn("closeWorkspaceJobBuilder();", workspace)
        self.assertIn('iframe.setAttribute("src", "");', workspace)
        self.assertIn('iframe.setAttribute("hidden", "hidden");', workspace)
        self.assertIn('iframe.classList.add("hidden");', workspace)

    def test_workspace_refresh_button_posts_missing_directory_cleanup(self):
        content = self._read_template("workspace.html")

        self.assertIn('<form method="POST" action="{% url \'nice_lite:refresh_workspace_jobs\' %}">', content)
        self.assertIn('title="refresh job cards"', content)
        self.assertIn("{% csrf_token %}", content)

    def test_workspace_has_job_style_title_bar_with_parent_project_back_link(self):
        workspace = self._read_template("workspace.html")
        nav_header = self._read_template("includes/_nav_header.html")

        self.assertIn(
            "{% include 'includes/_nav_header.html' with ",
            workspace,
        )
        self.assertIn("header_type='workspace'", workspace)
        self.assertIn("back_project_id=current_project_id", workspace)
        self.assertIn("back_target='_parent'", workspace)
        self.assertIn("back_aria_label='back to project'", workspace)
        self.assertIn("workspace_title=current_workspace_name", workspace)
        self.assertNotIn('<header class="bg-streambar', workspace)
        self.assertNotIn('viewBox="0 0 12 12"', workspace)
        self.assertIn("bg-streambar border-b border-streamline px-4 h-10", nav_header)
        self.assertIn("?selected_project_id={{ back_project_id }}", nav_header)
        self.assertIn('target="{{ back_target }}"', nav_header)
        self.assertIn('{% if header_type == "workspace" %}', nav_header)
        self.assertIn("{{ workspace_title }}</span>", nav_header)

        rendered_header = render_to_string(
            "includes/_nav_header.html",
            {
                "header_type": "workspace",
                "workspace_title": "workspace 2",
                "back_url": "/",
                "back_project_id": 7,
                "back_target": "_parent",
                "back_aria_label": "back to project",
            },
        )
        self.assertIn('href="/?selected_project_id=7"', rendered_header)
        self.assertIn('target="_parent"', rendered_header)
        self.assertIn('aria-label="back to project"', rendered_header)
        self.assertIn(">workspace 2</span>", rendered_header)

    def test_browser_messages_remain_visible(self):
        content = self._read_template("messages.html")
        index = self._read_template("index.html")
        project = self._read_template("project.html")
        workspace = self._read_template("workspace.html")
        stream = self._read_template("nice_stream/streamview.html")
        batch = self._read_template("nice_batch/batchview.html")
        manual_picker = self._read_template("nice_batch/manual_picker.html")
        jobbuilder = self._read_template("jobbuilder.html")

        self.assertIn('id="message_alert"', content)
        self.assertIn("w-full flex-shrink-0", content)
        self.assertNotIn("absolute", content)
        self.assertNotIn("setTimeout", content)
        self.assertNotIn('style.display="none"', content)
        message_include = "{% include 'messages.html' %}"
        for template, footer_marker in (
            (project, "<footer class="),
            (workspace, "<!-- Workspace details footer"),
            (stream, "{% include 'includes/_job_details.html' %}"),
            (batch, "{% include 'includes/_job_details.html' %}"),
            (manual_picker, "{% include 'includes/_job_details.html' %}"),
        ):
            self.assertLess(template.index(message_include), template.index(footer_marker))

        self.assertIn("{% if not project or current_workspace_id %}", index)
        self.assertIn("data-nice-message-footer", jobbuilder)
        self.assertEqual(jobbuilder.count("data-nice-message-footer"), 3)
        self.assertIn("{% include 'messages.html' %}", jobbuilder)
        self.assertIn("function positionBuilderMessages()", jobbuilder)
        self.assertEqual(jobbuilder.count("positionBuilderMessages();"), 2)
        self.assertNotIn("MutationObserver", content)
        self.assertNotIn("MutationObserver", jobbuilder)

    def test_jobbuilder_submits_to_named_workspace_iframe(self):
        jobbuilder = self._read_template("jobbuilder.html")
        index = self._read_template("index.html")
        self.assertIn('target="workspace_iframe"', jobbuilder)
        self.assertIn('name="workspace_iframe"', index)

    def test_project_selector_restores_its_most_recent_workspace(self):
        index = self._read_template("index.html")

        self.assertIn('onchange="selectProjectWorkspace(this)"', index)
        self.assertIn('data-recent-workspace-input', index)
        self.assertIn('const workspaceStorageKey = (projectId) => `${WORKSPACE_KEY}:${projectId}`;', index)
        self.assertIn('const workspaceId = rememberedWorkspace(select.value);', index)
        self.assertIn('workspaceInput.value = workspaceId || "";', index)
        self.assertIn(
            "const wantsWorkspaceId = !params.has(PROJECT_KEY) && !params.has(WORKSPACE_KEY)",
            index,
        )
        self.assertIn(
            'sessionStorage.setItem(workspaceStorageKey(currentProjectId), currentWorkspaceId);',
            index,
        )

    def test_project_page_reuses_nav_header_with_previous_project_back_control(self):
        project = self._read_template("project.html")
        nav_header = self._read_template("includes/_nav_header.html")
        index = self._read_template("index.html")

        self.assertIn(
            "{% include 'includes/_nav_header.html' with "
            "header_type='project' project_title=project.name %}",
            project,
        )
        self.assertIn('data-project-back disabled aria-disabled="true"', nav_header)
        self.assertIn('title="no previous project"', nav_header)
        self.assertIn('{% if header_type == "project" %}', nav_header)
        self.assertIn("{{ project_title }}</span>", nav_header)
        self.assertIn('const PREVIOUS_PROJECT_KEY = "previous_project_id";', index)
        self.assertIn("sessionStorage.setItem(PREVIOUS_PROJECT_KEY, storedProjectId);", index)
        self.assertIn('document.querySelector("[data-project-back]")', index)
        self.assertIn("projectBackButton.disabled = false;", index)
        self.assertIn("previousProjectId !== currentProjectId", index)

        rendered_header = render_to_string(
            "includes/_nav_header.html",
            {
                "header_type": "project",
                "project_title": "project 2",
            },
        )
        self.assertIn("data-project-back", rendered_header)
        self.assertIn(" disabled", rendered_header)
        self.assertIn('aria-disabled="true"', rendered_header)
        self.assertIn(">project 2</span>", rendered_header)

    def test_projects_page_uses_project_details_and_workspace_cards(self):
        projects = self._read_template("project.html")

        self.assertIn("Project details", projects)
        self.assertIn("{% for workspace in workspaces %}", projects)
        self.assertIn("{{ project.name }}", projects)
        self.assertIn("selected_project_id={{ project.id }}", projects)
        self.assertIn("{% url 'nice_lite:create_workspace' %}", projects)
        self.assertIn("border-streamring/20 hover:border-streamring", projects)
        self.assertIn("w-[260px] h-[361px]", projects)
        self.assertIn("border-b border-streamdivider", projects)
        self.assertIn("border-t border-streamdivider", projects)
        self.assertIn("workspace folder", projects)
        self.assertIn("includes/_nav_header.html", projects)
        self.assertNotIn('<header class="bg-streambar', projects)
        self.assertNotIn("<!DOCTYPE html>", projects)
        self.assertNotIn("<body", projects)
        self.assertNotIn('aria-label="project and workspace navigation"', projects)
        self.assertNotIn("workspace_jobs", projects)
        self.assertNotIn("job_builder_iframe", projects)

        index = self._read_template("index.html")
        self.assertIn("{% include 'project.html' %}", index)

    def test_file_browser_openers_forward_current_input_path(self):
        jobbuilder = self._read_template("jobbuilder.html")
        newproject = self._read_template("newproject.html")
        classic_newjob = self._read_template("nice_classic/newjob.html")

        self.assertIn('params.set("selectedpath", selectedPath)', jobbuilder)
        self.assertIn('jobsIframe.setAttribute("src", browserUrl + "?" + params.toString())', jobbuilder)
        self.assertIn('onclick="openProjectDirectoryBrowser(this)"', newproject)
        self.assertIn('params.set("selectedpath", selectedPath)', newproject)
        self.assertIn('onclick="openClassicFileBrowser(this, \'dir\')"', classic_newjob)
        self.assertIn('onclick="openClassicFileBrowser(this, \'file\')"', classic_newjob)
        self.assertIn('new URLSearchParams({ selectedpath: selectedPath })', classic_newjob)

    def test_jobbuilder_converts_dropped_file_uris_to_paths(self):
        jobbuilder = self._read_template("jobbuilder.html")

        context = {
            "stream_user_inputs": [
                {"key": "movies", "keytype": "dir", "label": "movies"},
                {"key": "gain", "keytype": "file", "label": "gain"},
                {"key": "threads", "keytype": "int", "label": "threads"},
            ],
            "simple_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "simple_program_inputs": [{
                "prg": "demo",
                "disp": "Demo",
                "sections": [{
                    "name": "inputs",
                    "inputs": [
                        {"key": "directory", "keytype": "dir", "label": "directory"},
                        {"key": "document", "keytype": "file", "label": "document"},
                        {"key": "count", "keytype": "int", "label": "count"},
                    ],
                }],
            }],
            "single_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "single_program_inputs": [{
                "prg": "demo",
                "disp": "Demo",
                "sections": [{
                    "name": "inputs",
                    "inputs": [
                        {"key": "directory", "keytype": "dir", "label": "directory"},
                        {"key": "document", "keytype": "file", "label": "document"},
                        {"key": "count", "keytype": "int", "label": "count"},
                    ],
                }],
            }],
            "default_batch_project_file": "/workspace/workspace.simple",
        }
        rendered = render_to_string("jobbuilder.html", context)

        self.assertEqual(jobbuilder.count(" data-path-input"), 8)
        for input_id in (
            "field_movies",
            "field_gain",
            "simple_batch_project_file",
            "batch_demo_directory",
            "batch_demo_document",
            "single_batch_project_file",
            "single_demo_directory",
            "single_demo_document",
        ):
            self.assertRegex(rendered, rf'<input id="{input_id}"[^>]*data-path-input')
        for input_id in ("field_threads", "batch_demo_count", "single_demo_count"):
            self.assertNotRegex(rendered, rf'<input id="{input_id}"[^>]*data-path-input')
        self.assertIn('dataTransfer.getData("text/uri-list")', jobbuilder)
        self.assertIn('dataTransfer.getData("text/plain")', jobbuilder)
        self.assertIn('dataTransfer.getData("text/x-moz-url-data")', jobbuilder)
        self.assertIn('decodeURIComponent(fileUrl.pathname)', jobbuilder)
        self.assertIn('document.addEventListener("drop"', jobbuilder)
        self.assertIn('input.matches("[data-path-input]")', jobbuilder)
        self.assertIn('function applyDroppedFilePath(input, path)', jobbuilder)
        self.assertIn('window.setTimeout(function()', jobbuilder)
        self.assertIn('const insertedPath = fileUriToPath(input.value);', jobbuilder)
        self.assertIn('input.dispatchEvent(new Event("input", {bubbles: true}));', jobbuilder)
        self.assertEqual(jobbuilder.count(" data-batch-project-drop-target"), 2)
        self.assertIn('const BATCH_PROJECT_DRAG_TYPE = "application/x-nice-batch-project";', jobbuilder)
        self.assertIn('input.matches("[data-batch-project-drop-target]")', jobbuilder)
        self.assertIn('event.dataTransfer.getData(BATCH_PROJECT_DRAG_TYPE)', jobbuilder)
        self.assertIn('applyDroppedFilePath(input, batchProjectPath);', jobbuilder)

    def test_new_project_back_button_closes_form_in_parent_shell(self):
        newproject = self._read_template("newproject.html")

        self.assertIn("{% url 'nice_lite:close_new_project' %}", newproject)
        self.assertIn('target="_parent" title="back"', newproject)
        self.assertNotIn("href={% url 'nice_lite:workspace' %}", newproject)

    def test_file_browser_navigation_passes_paths_as_query_parameters(self):
        filebrowser = self._read_template("filebrowser.html")

        self.assertNotIn("{% url 'nice_lite:file_browser' type parentdir %}", filebrowser)
        self.assertNotIn("{% url 'nice_lite:file_browser' type path|add:'/'|add:dir %}", filebrowser)
        self.assertIn('name="selectedpath" value="{{ parentdir }}"', filebrowser)

    def test_batch_file_browser_restore_does_not_show_loading_overlay(self):
        filebrowser = self._read_template("filebrowser.html")

        self.assertIn("const restoreJobsIframe = (iframe) => {", filebrowser)
        self.assertIn('restoreUrl.searchParams.set("force", "1");', filebrowser)
        self.assertNotIn("parent.toggleJobsIframeLoading(true);", filebrowser)

    def test_file_browser_openers_share_last_directory(self):
        filebrowser = self._read_template("filebrowser.html")
        jobbuilder = self._read_template("jobbuilder.html")
        newproject = self._read_template("newproject.html")
        classic_newjob = self._read_template("nice_classic/newjob.html")

        self.assertIn('localStorage.setItem(FILE_BROWSER_LAST_DIRECTORY_KEY, directoryToRemember)', filebrowser)
        for opener in (jobbuilder, newproject, classic_newjob):
            self.assertIn('localStorage.getItem("niceFileBrowserLastDirectory")', opener)
            self.assertIn('params.set("remembered", "1")', opener)

    def test_batch_program_selectors_have_search_filters(self):
        jobbuilder = self._read_template("jobbuilder.html")

        self.assertIn('id="simple_program_search"', jobbuilder)
        self.assertIn('placeholder="search simple programs"', jobbuilder)
        self.assertIn('oninput="filterSimplePrograms(this.value)"', jobbuilder)
        self.assertIn('data-simple-program-search=', jobbuilder)
        self.assertIn('id="simple_program_no_results"', jobbuilder)
        self.assertIn('function filterSimplePrograms(query)', jobbuilder)

        self.assertNotIn("select a single program", jobbuilder)
        self.assertIn('id="single_program_search"', jobbuilder)
        self.assertIn('placeholder="search single programs"', jobbuilder)
        self.assertIn('oninput="filterSinglePrograms(this.value)"', jobbuilder)
        self.assertIn('data-single-program-search=', jobbuilder)
        self.assertIn('id="single_program_no_results"', jobbuilder)
        self.assertIn('function filterSinglePrograms(query)', jobbuilder)

        self.assertGreaterEqual(jobbuilder.count('absolute right-2'), 2)
        self.assertIn('program.classList.toggle("hidden", !matches)', jobbuilder)

    def test_batch_developer_mode_uses_commander_and_argument_names(self):
        jobbuilder = self._read_template("jobbuilder.html")
        context = {
            "stream_user_inputs": [{
                "key": "nthr",
                "keytype": "int",
                "label": "Stream thread count",
            }],
            "simple_programs": [{"prg": "simple_demo", "disp": "Friendly SIMPLE", "desc": ""}],
            "simple_program_inputs": [{
                "prg": "simple_demo",
                "disp": "Friendly SIMPLE",
                "sections": [{
                    "name": "compute",
                    "inputs": [{
                        "key": "nthr",
                        "keytype": "int",
                        "label": "Number of threads",
                    }],
                }],
            }],
            "single_programs": [{"prg": "single_demo", "disp": "Friendly SINGLE", "desc": ""}],
            "single_program_inputs": [{
                "prg": "single_demo",
                "disp": "Friendly SINGLE",
                "sections": [{
                    "name": "input_output",
                    "inputs": [{
                        "key": "mskdiam",
                        "keytype": "float",
                        "label": "Mask diameter",
                    }],
                }],
            }],
            "default_batch_project_file": "/workspace/workspace.simple",
        }
        rendered = render_to_string("jobbuilder.html", context)

        self.assertEqual(jobbuilder.count('id="batch_visibility_control"'), 1)
        self.assertLess(
            jobbuilder.index('id="batch_visibility_control"'),
            jobbuilder.index('id="tab_stream"'),
        )
        self.assertIn(
            'data-batch-commander-label data-friendly-name="Friendly SIMPLE" '
            'data-commander-name="simple_demo">Friendly SIMPLE</span>',
            rendered,
        )
        self.assertIn(
            'data-batch-commander-label data-friendly-name="Friendly SINGLE" '
            'data-commander-name="single_demo">Friendly SINGLE</span>',
            rendered,
        )
        self.assertEqual(
            rendered.count('<span class="block text-xs font-medium text-streamtext" data-batch-commander-label'),
            2,
        )
        self.assertIn('for="field_nthr">Stream thread count</label>', rendered)
        self.assertIn(
            'data-friendly-label="Number of threads" data-argument-key="nthr">Number of threads</span>',
            rendered,
        )
        self.assertIn(
            'data-friendly-label="Mask diameter" data-argument-key="mskdiam">Mask diameter</span>',
            rendered,
        )
        self.assertEqual(rendered.count('>input project</label>'), 2)
        self.assertEqual(rendered.count("<span data-batch-argument-label"), 2)
        self.assertIn('const showDeveloperNames = mode === "developer";', jobbuilder)
        self.assertIn('document.querySelectorAll("[data-batch-commander-label]")', jobbuilder)
        self.assertIn('? label.dataset.commanderName', jobbuilder)
        self.assertIn(': label.dataset.friendlyName;', jobbuilder)
        self.assertIn('? label.dataset.argumentKey', jobbuilder)
        self.assertIn(': label.dataset.friendlyLabel;', jobbuilder)
        self.assertIn('updateSelectedProgramTitle("simple", selectedSimpleProgramKey);', jobbuilder)
        self.assertIn('updateSelectedProgramTitle("single", selectedSingleProgramKey);', jobbuilder)

    def test_batch_tabs_restore_selected_commanders(self):
        jobbuilder = self._read_template("jobbuilder.html")

        self.assertIn("let selectedSimpleProgramKey = null;", jobbuilder)
        self.assertIn("showSimpleProgramConfig(selectedSimpleProgramKey);", jobbuilder)
        self.assertIn("selectedSimpleProgramKey = programKey;", jobbuilder)
        self.assertIn("selectedSimpleProgramKey = null;", jobbuilder)
        self.assertIn("let selectedSingleProgramKey = null;", jobbuilder)
        self.assertIn("showSingleProgramConfig(selectedSingleProgramKey);", jobbuilder)
        self.assertIn("selectedSingleProgramKey = programKey;", jobbuilder)
        self.assertIn("selectedSingleProgramKey = null;", jobbuilder)

    def test_batch_start_submits_selected_registered_program_form(self):
        jobbuilder = self._read_template("jobbuilder.html")

        self.assertEqual(jobbuilder.count("action=\"{% url 'nice_lite:create_batch' %}\""), 2)
        self.assertIn('data-batch-package="simple"', jobbuilder)
        self.assertIn('data-batch-package="single"', jobbuilder)
        self.assertIn('name="package" value="simple"', jobbuilder)
        self.assertIn('name="package" value="single"', jobbuilder)
        self.assertIn('name="program" value="{{ program.prg }}"', jobbuilder)
        self.assertIn('onclick="submitSelectedBatch(\'simple\')"', jobbuilder)
        self.assertIn('onclick="submitSelectedBatch(\'single\')"', jobbuilder)
        self.assertIn('function submitSelectedBatch(packageName)', jobbuilder)
        self.assertIn('if (!form) return;', jobbuilder)
        self.assertIn('sourceInput.value = sourceSelector.value;', jobbuilder)
        self.assertIn('projectFile.toLowerCase().endsWith(".simple")', jobbuilder)
        self.assertIn('projectFileInput.value = projectFile;', jobbuilder)
        self.assertIn('form.requestSubmit();', jobbuilder)
        self.assertIn('function validateSelectedBatchForm(packageName, programKey, form)', jobbuilder)
        self.assertIn('form.reportValidity();', jobbuilder)
        self.assertIn('complete required fields:', jobbuilder)
        self.assertIn('data-batch-requirement', jobbuilder)
        self.assertIn('requirement.dataset.requirementKeys', jobbuilder)
        self.assertIn('programKey === "import_movies"', jobbuilder)
        self.assertIn('choose a movie-file list or an input movies directory', jobbuilder)
        self.assertIn('id="simple_batch_feedback" role="alert"', jobbuilder)
        self.assertIn('id="single_batch_feedback" role="alert"', jobbuilder)
        self.assertEqual(jobbuilder.count('class="text-streamerror" aria-hidden="true"> *</span>'), 4)
        self.assertIn('const batchSubmissionPending = {simple: false, single: false};', jobbuilder)
        self.assertIn('form.addEventListener("submit", markBatchSubmissionPending);', jobbuilder)
        self.assertIn('button.textContent = "starting...";', jobbuilder)

    def test_batch_project_validation_uses_visible_selection(self):
        jobbuilder = self._read_template("jobbuilder.html")

        self.assertIn('data-requires-project="{{ program.requires_project|yesno:', jobbuilder)
        self.assertIn('form.dataset.requiresProject === "true"', jobbuilder)
        self.assertIn('if (key === "projfile") return Boolean(projectFile);', jobbuilder)

    def test_batch_submission_guard_does_not_change_stream_submission(self):
        jobbuilder = self._read_template("jobbuilder.html")

        self.assertIn("action=\"{% url 'nice_lite:create_stream' %}\"", jobbuilder)
        self.assertIn('<button type="submit"', jobbuilder)
        self.assertEqual(jobbuilder.count('data-batch-package='), 2)
        self.assertIn('document.querySelectorAll("form[data-batch-package]")', jobbuilder)

    def test_batch_project_file_selector_starts_empty_without_prefill(self):
        rendered = render_to_string("jobbuilder.html", {
            "stream_user_inputs": [],
            "simple_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "simple_program_inputs": [{"prg": "demo", "disp": "Demo", "sections": []}],
            "single_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "single_program_inputs": [{"prg": "demo", "disp": "Demo", "sections": []}],
            "default_batch_project_file": "",
        })

        self.assertIn('id="simple_batch_project_file"', rendered)
        self.assertIn('id="single_batch_project_file"', rendered)
        self.assertNotIn('id="simple_batch_source"', rendered)
        self.assertNotIn('id="single_batch_source"', rendered)
        self.assertEqual(
            rendered.count('name="batch_project_file" value=""'),
            2,
        )
        self.assertEqual(rendered.count('data-required-extension=".simple"'), 2)
        self.assertEqual(
            rendered.count("openWorkspaceBrowser('file', 'simple_batch_project_file', '.simple', 'project_file')"),
            1,
        )
        self.assertEqual(
            rendered.count("openWorkspaceBrowser('file', 'single_batch_project_file', '.simple', 'project_file')"),
            1,
        )

    def test_file_browser_preserves_and_enforces_extension_filter(self):
        filebrowser = self._read_template("filebrowser.html")

        self.assertGreaterEqual(
            filebrowser.count('name="extension" value="{{ required_extension }}"'),
            3,
        )
        self.assertIn('const requiredExtension = "{{ required_extension|escapejs }}";', filebrowser)
        self.assertIn('!selectedpath.value.toLowerCase().endsWith(requiredExtension)', filebrowser)

    def test_batch_inputs_render_current_ui_json_schema(self):
        context = {
            "stream_user_inputs": [],
            "simple_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "simple_program_inputs": [{
                "prg": "demo",
                "disp": "Demo",
                "sections": [{
                    "name": "inputs",
                    "inputs": [
                        {"key": "mode", "keytype": "multi", "label": "mode", "options": ["one", "two"], "has_default": True, "default": "two"},
                        {"key": "count", "keytype": "int", "label": "count", "has_default": True, "default": 4},
                        {"key": "offset", "keytype": "num", "label": "offset", "has_default": True, "default": -1.5},
                    ],
                }],
            }],
            "single_programs": [],
            "single_program_inputs": [],
        }

        rendered = render_to_string("jobbuilder.html", context)

        self.assertIn('<form id="batch_form_simple_demo"', rendered)
        self.assertIn('<option value="two"\n                                                                    selected', rendered)
        self.assertIn('placeholder="4"', rendered)
        self.assertIn('name="offset" type="number" step="any"', rendered)
        self.assertIn('placeholder="-1.5"', rendered)
        self.assertNotIn('option value=""', rendered)

    def test_batch_visibility_hides_sections_without_visible_fields(self):
        context = {
            "stream_user_inputs": [],
            "simple_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "simple_program_inputs": [{
                "prg": "demo",
                "disp": "Demo",
                "sections": [
                    {
                        "name": "standard",
                        "inputs": [{"key": "visible", "keytype": "str", "label": "visible"}],
                    },
                    {
                        "name": "advanced",
                        "inputs": [{
                            "key": "hidden_in_standard",
                            "keytype": "str",
                            "label": "hidden in standard",
                            "visibility": "advanced",
                        }],
                    },
                ],
            }],
            "single_programs": [],
            "single_program_inputs": [],
        }

        rendered = render_to_string("jobbuilder.html", context)
        jobbuilder = self._read_template("jobbuilder.html")

        self.assertEqual(rendered.count("data-batch-section"), 2)
        self.assertIn('data-batch-field-visibility="advanced"', rendered)
        self.assertIn('document.querySelectorAll("[data-batch-section]")', jobbuilder)
        self.assertIn('section.querySelectorAll("[data-batch-field-visibility]")', jobbuilder)
        self.assertIn('section.classList.toggle("hidden", !hasVisibleField);', jobbuilder)

    def test_batch_required_argument_labels_are_bold(self):
        context = {
            "stream_user_inputs": [],
            "simple_programs": [{"prg": "simple_demo", "disp": "Simple Demo", "desc": ""}],
            "simple_program_inputs": [{
                "prg": "simple_demo",
                "disp": "Simple Demo",
                "sections": [{
                    "name": "inputs",
                    "inputs": [
                        {"key": "required_arg", "keytype": "str", "label": "required simple", "required": True},
                        {"key": "optional_arg", "keytype": "str", "label": "optional simple", "required": False},
                    ],
                }],
            }],
            "single_programs": [{"prg": "single_demo", "disp": "Single Demo", "desc": ""}],
            "single_program_inputs": [{
                "prg": "single_demo",
                "disp": "Single Demo",
                "sections": [{
                    "name": "inputs",
                    "inputs": [
                        {"key": "required_arg", "keytype": "str", "label": "required single", "required": True},
                        {"key": "optional_arg", "keytype": "str", "label": "optional single", "required": False},
                    ],
                }],
            }],
            "batch_project_sources": [],
            "default_batch_source": "workspace",
        }

        rendered = render_to_string("jobbuilder.html", context)

        self.assertEqual(
            rendered.count('class="text-xs font-bold text-streamtext whitespace-nowrap"'),
            4,
        )
        self.assertIn(
            'for="batch_simple_demo_required_arg"><span data-batch-argument-label '
            'data-friendly-label="required simple" data-argument-key="required_arg">required simple</span>',
            rendered,
        )
        self.assertIn(
            'for="single_single_demo_required_arg"><span data-batch-argument-label '
            'data-friendly-label="required single" data-argument-key="required_arg">required single</span>',
            rendered,
        )
        self.assertIn(
            'class="text-xs font-medium text-streamtext whitespace-nowrap" '
            'for="batch_simple_demo_optional_arg"><span data-batch-argument-label '
            'data-friendly-label="optional simple" data-argument-key="optional_arg">optional simple</span>',
            rendered,
        )
        self.assertIn(
            'class="text-xs font-medium text-streamtext whitespace-nowrap" '
            'for="single_single_demo_optional_arg"><span data-batch-argument-label '
            'data-friendly-label="optional single" data-argument-key="optional_arg">optional single</span>',
            rendered,
        )

    def test_batch_input_requirements_are_rendered_for_client_validation(self):
        context = {
            "stream_user_inputs": [],
            "simple_programs": [{"prg": "new_project", "disp": "New project", "desc": ""}],
            "simple_program_inputs": [{
                "prg": "new_project",
                "disp": "New project",
                "requirements": [{
                    "label": "Project location",
                    "help": "Supply a project name, a project directory, or both.",
                    "min_selected": 1,
                    "max_selected": 2,
                    "keys": ["projname", "dir"],
                }],
                "sections": [],
            }],
            "single_programs": [],
            "single_program_inputs": [],
            "batch_project_sources": [],
            "default_batch_source": "workspace",
        }

        rendered = render_to_string("jobbuilder.html", context)

        self.assertIn("data-batch-requirement", rendered)
        self.assertIn('data-requirement-keys="projname,dir"', rendered)
        self.assertIn("Supply a project name, a project directory, or both.", rendered)

    def test_batch_panels_match_stream_spacing(self):
        jobbuilder = self._read_template("jobbuilder.html")

        heading_classes = "text-xs font-semibold text-streamtext text-left mb-6 mx-4 uppercase tracking-widest"
        footer_classes = "flex items-center justify-end mt-3 p-2 border-t border-streamline gap-2"
        button_classes = "text-xs px-6 py-1.5 rounded-lg font-medium bg-streamaccent text-white hover:bg-streamring transition-colors"
        self.assertEqual(jobbuilder.count(heading_classes), 3)
        self.assertEqual(jobbuilder.count(footer_classes), 3)
        self.assertEqual(jobbuilder.count(button_classes), 6)
        self.assertIn('id="panel_batch" class="flex-1 w-full flex flex-col hidden"', jobbuilder)
        self.assertIn('id="card_simple_batch" class="flex-1 w-full flex flex-col"', jobbuilder)
        self.assertIn('id="card_single_batch" class="flex-1 w-full flex flex-col"', jobbuilder)
        self.assertIn('id="simple_config_panel" class="hidden flex flex-1 flex-col mx-4 pb-16"', jobbuilder)
        self.assertIn('id="single_config_panel" class="hidden flex flex-1 flex-col mx-4 pb-16"', jobbuilder)
        fixed_footer_classes = f"{footer_classes} fixed bottom-0 left-0 w-full z-10 bg-streambar hidden"
        self.assertIn(f'id="simple_batch_actions" class="{fixed_footer_classes}"', jobbuilder)
        self.assertIn(f'id="single_batch_actions" class="{fixed_footer_classes}"', jobbuilder)
        self.assertIn('id="simple_batch_start" type="button"', jobbuilder)
        self.assertIn('id="single_batch_start" type="button"', jobbuilder)
        self.assertEqual(jobbuilder.count('onclick="closeJobBuilder()"'), 3)
        self.assertEqual(jobbuilder.count('if (actions) actions.classList.add("hidden");'), 2)
        self.assertEqual(jobbuilder.count('if (actions) actions.classList.remove("hidden");'), 2)

    def test_batch_cards_reuse_stream_stop_and_delete_controls(self):
        batch_card = self._read_template("nice_classic/_batch_card.html")
        batch_footer = self._read_template("includes/_job_card_footer.html")
        jobs = self._read_template("jobs_cards.html")
        jobs_table = self._read_template("jobs_table.html")

        self.assertIn("{% url 'nice_lite:stop_batch' %}", batch_card)
        self.assertIn('onclick="stopBatchJob(this)"', batch_card)
        self.assertIn('<circle cx="8" cy="8" r="6"/>', batch_card)
        self.assertIn('<rect x="5.5" y="5.5" width="5" height="5"', batch_card)
        self.assertIn("{% url 'nice_lite:delete_batch' %}", batch_card)
        self.assertIn('onclick="deleteBatchJob(this)"', batch_card)
        self.assertIn('<polyline points="2,4 14,4"/>', batch_card)
        self.assertIn('<path d="M5 4V3a1 1 0 0 1 1-1h4a1 1 0 0 1 1 1v1M6 7v5M10 7v5M3 4l1 9a1 1 0 0 0 1 1h6a1 1 0 0 0 1-1l1-9"/>', batch_card)
        self.assertIn("const stopBatchJob = (button) => {", jobs)
        self.assertIn("{% url 'nice_lite:finish_batch' %}", batch_footer)
        self.assertIn('onclick="markBatchJobFinished(this)"', batch_footer)
        self.assertIn("const markBatchJobFinished = (button) => {", jobs)
        self.assertIn("{% url 'nice_lite:clone_batch' %}", batch_card)
        self.assertIn('target="job_builder_iframe"', batch_card)
        self.assertIn('<button type="submit"', batch_card)
        self.assertNotIn("const rerunBatchJob = (button) => {", jobs)
        self.assertIn("{% url 'nice_lite:clone_batch' %}", jobs_table)
        self.assertIn('target="job_builder_iframe"', jobs_table)
        self.assertNotIn("const rerunBatchJob = (button) => {", jobs_table)
        self.assertIn("const deleteBatchJob = (button) => {", jobs)
        self.assertIn("permanently delete batch job", jobs)
        self.assertIn("This cannot be undone.", jobs)

    def test_batch_footer_marks_only_queued_or_failed_jobs_finished(self):
        job = {
            "id": 7,
            "pckg": "simple",
            "prog": "import_movies",
            "status": "queued",
        }

        queued = render_to_string(
            "includes/_job_card_footer.html",
            {"job": job, "card_type": "batch"},
        )
        job["status"] = "failed"
        failed = render_to_string(
            "includes/_job_card_footer.html",
            {"job": job, "card_type": "batch"},
        )
        self.assertIn('title="mark batch job finished"', queued)
        self.assertIn('action="/finishbatch"', queued)
        self.assertIn('name="mark_finished" value="1"', queued)
        self.assertIn('title="mark batch job finished"', failed)

        for status in ("running", "finished", "stopped"):
            with self.subTest(status=status):
                job["status"] = status
                rendered = render_to_string(
                    "includes/_job_card_footer.html",
                    {"job": job, "card_type": "batch"},
                )
                self.assertNotIn('title="mark batch job finished"', rendered)

    def test_batch_rerun_form_opens_selected_program_with_saved_values(self):
        rendered = render_to_string("jobbuilder.html", {
            "stream_user_inputs": [],
            "simple_programs": [{
                "prg": "cluster2D",
                "disp": "Create 2D Classes",
                "desc": "",
            }],
            "simple_program_inputs": [{
                "prg": "cluster2D",
                "disp": "Create 2D Classes",
                "requirements": [],
                "rerun_of": 17,
                "sections": [{
                    "name": "compute",
                    "inputs": [{
                        "key": "nthr",
                        "keytype": "int",
                        "label": "Threads",
                        "value": "8",
                    }],
                }],
            }],
            "single_programs": [],
            "single_program_inputs": [],
            "default_batch_project_file": "/workspace/input.simple",
            "batch_prefill": {
                "job_id": 17,
                "package": "simple",
                "program": "cluster2D",
            },
        })

        self.assertIn('name="rerun_of" value="17"', rendered)
        self.assertIn(
            'id="batch_cluster2D_nthr" name="nthr" type="number" min="0" step="1"',
            rendered,
        )
        self.assertIn('value="8"', rendered)
        self.assertIn('value="/workspace/input.simple"', rendered)
        self.assertIn('const batchPrefillPackage = "simple";', rendered)
        self.assertIn('const batchPrefillProgram = "cluster2D";', rendered)
        self.assertIn('selectBuilderTab("simple_batch");', rendered)

    def test_batch_card_controls_follow_job_status(self):
        job = {
            "id": 7,
            "disp": 1,
            "name": "import movies",
            "dirc": "1_import_movies",
            "args": {},
            "pckg": "simple",
            "prog": "import_movies",
            "master_stats": {},
            "status": "running",
        }

        context = {"job": job}
        running = render_to_string("nice_classic/_batch_card.html", context)
        job["status"] = "queued"
        queued = render_to_string("nice_classic/_batch_card.html", context)
        job["status"] = "finished"
        finished = render_to_string("nice_classic/_batch_card.html", context)
        self.assertIn('title="stop batch job"', running)
        self.assertNotIn('title="delete batch job"', running)
        self.assertIn("rounded-full bg-streamaction/10 text-streamaction", running)
        self.assertNotIn("w-1.5 h-1.5", running)
        self.assertLess(
            running.index('title="stop batch job"'),
            running.index('job directory'),
        )
        self.assertNotIn('title="stop batch job"', queued)
        self.assertNotIn('title="rerun batch job"', queued)
        self.assertIn('title="delete batch job"', queued)
        self.assertIn("rounded-full bg-streamline/40 text-streamlabel", queued)
        self.assertNotIn('title="stop batch job"', finished)
        self.assertIn('title="rerun batch job"', finished)
        self.assertIn('title="delete batch job"', finished)
        self.assertIn("rounded-full bg-streamring/10 text-streamaccent", finished)
        self.assertLess(
            finished.index('title="delete batch job"'),
            finished.index('job directory'),
        )

    def test_batch_card_always_links_to_detail_view(self):
        job = {
            "id": 7,
            "disp": 1,
            "name": "Import Movie Data",
            "dirc": "1_import_movies",
            "args": {},
            "pckg": "simple",
            "prog": "import_movies",
            "master_stats": {},
            "status": "finished",
        }

        rendered = render_to_string("nice_classic/_batch_card.html", {"job": job})

        self.assertIn('href="/viewbatch/7"', rendered)
        self.assertIn('target="_parent"', rendered)
        self.assertIn('aria-label="view batch job 1 Import Movie Data"', rendered)
        self.assertNotIn('<circle cx="8" cy="8" r="3"', rendered)

        job["name"] = "Create 2D Class Averages"
        job["prog"] = "abinitio2D"
        class_average_rendered = render_to_string(
            "nice_classic/_batch_card.html",
            {"job": job},
        )

        self.assertIn(
            'href="/viewbatch/7?class_selector=1"',
            class_average_rendered,
        )

        job["name"] = "Initial 3D Reconstruction"
        job["prog"] = "abinitio3D"
        volume_rendered = render_to_string(
            "nice_classic/_batch_card.html",
            {"job": job},
        )

        self.assertIn(
            'href="/viewbatch/7?volume_viewer=1#batch_volume_viewer"',
            volume_rendered,
        )

    def test_active_abinitio3d_batch_card_opens_volume_output(self):
        job = {
            "id": 7,
            "disp": 9,
            "name": "Initial 3D Reconstruction",
            "dirc": "9_abinitio3D",
            "args": {},
            "pckg": "simple",
            "prog": "abinitio3D",
            "master_stats": {},
            "status": "finished",
        }

        rendered = render_to_string(
            "nice_batch/includes/_batch_card.html",
            {"job": job},
        )

        self.assertIn('action="/viewbatch/7"', rendered)
        self.assertIn('name="volume_viewer" value="1"', rendered)

        job["prog"] = "abinitio2D"
        other_output = render_to_string(
            "nice_batch/includes/_batch_card.html",
            {"job": job},
        )
        self.assertNotIn('name="volume_viewer"', other_output)

    def test_finished_batch_card_defers_drag_state_until_builder_is_open(self):
        project_path = "/workspace/7_abinitio2D/workspace.simple"
        job = {
            "id": 7,
            "disp": 7,
            "name": "Create 2D Class Averages",
            "dirc": "7_abinitio2D",
            "args": {},
            "pckg": "simple",
            "prog": "abinitio2D",
            "master_stats": {},
            "status": "finished",
            "batch_project_drag_path": project_path,
        }

        rendered = render_to_string(
            "nice_batch/includes/_batch_card.html",
            {"job": job},
        )
        opening_tag = rendered.split(">", 1)[0]

        self.assertIn("cursor cursor-pointer", opening_tag)
        self.assertNotIn("cursor-grab", opening_tag)
        self.assertNotIn("draggable=", opening_tag)
        self.assertIn(f'data-batch-project-path="{project_path}"', rendered)
        self.assertIn('ondragstart="dragBatchProject(event)"', rendered)

        job["batch_project_drag_path"] = ""
        unavailable = render_to_string(
            "nice_batch/includes/_batch_card.html",
            {"job": job},
        )
        self.assertNotIn('data-batch-project-path=', unavailable)
        self.assertNotIn('ondragstart="dragBatchProject(event)"', unavailable)

    def test_job_card_project_drag_keeps_artifact_payload_separate(self):
        jobs = self._read_template("jobs_cards.html")

        self.assertIn('const dragBatchProject = (event) => {', jobs)
        self.assertIn('source.closest(".artifact-drag-row [draggable=\'true\']")', jobs)
        self.assertIn('event.dataTransfer.setData(BATCH_PROJECT_DRAG_TYPE, projectPath);', jobs)
        self.assertIn('event.dataTransfer.setData("application/json", payload);', jobs)
        self.assertIn('event.dataTransfer.setData("text/plain", payload);', jobs)
        self.assertIn('const isJobBuilderActive = () => {', jobs)
        self.assertIn('card.draggable = active;', jobs)
        self.assertIn('card.classList.toggle("cursor-pointer", !active);', jobs)
        self.assertIn('card.classList.toggle("cursor-grab", active);', jobs)

    def test_job_builder_artifact_icons_overlay_only_job_card_footers(self):
        stream_card = self._read_template("nice_stream/includes/_stream_card.html")
        batch_card = self._read_template("nice_batch/includes/_batch_card.html")
        footer = self._read_template("includes/_job_card_footer.html")
        artifact_row = self._read_template("includes/_artifact_drag_row.html")

        self.assertNotIn("artifact-drag-row", stream_card)
        self.assertNotIn("artifact-drag-row", batch_card)
        self.assertEqual(
            footer.count("{% include 'includes/_artifact_drag_row.html' %}"),
            2,
        )
        self.assertIn(
            'class="artifact-drag-row hidden absolute inset-0',
            artifact_row,
        )
        self.assertEqual(artifact_row.count(" title="), 6)
        self.assertEqual(artifact_row.count(" aria-label="), 6)
        self.assertNotIn("group-hover:", artifact_row)

    def test_scale_volume_and_reproject_stack_icons_carry_owned_output_paths(self):
        scale_output = render_to_string(
            "includes/_artifact_drag_row.html",
            {"job": {
                "id": 7,
                "artifact_drag_paths": {"volume3D": "/workspace/7_scale/outvol.mrc"},
            }},
        )
        reproject_output = render_to_string(
            "includes/_artifact_drag_row.html",
            {"job": {
                "id": 8,
                "artifact_drag_paths": {"particles": "/workspace/8_reproject/reprojs.mrcs"},
            }},
        )

        self.assertIn('data-artifact-path="/workspace/7_scale/outvol.mrc" ondragstart="dragArtifact(event, \'7\', \'volume3D\')"', scale_output)
        self.assertNotIn('data-artifact-path="/workspace/7_scale/outvol.mrc" ondragstart="dragArtifact(event, \'7\', \'particles\')"', scale_output)
        self.assertIn('data-artifact-path="/workspace/8_reproject/reprojs.mrcs" ondragstart="dragArtifact(event, \'8\', \'particles\')"', reproject_output)

    def test_job_builder_accepts_typed_volume_and_particle_stack_drops(self):
        jobbuilder = self._read_template("jobbuilder.html")

        self.assertEqual(jobbuilder.count('data-artifact-drop-type="volume3D"'), 2)
        self.assertEqual(jobbuilder.count('data-artifact-drop-type="particles"'), 2)
        self.assertEqual(
            jobbuilder.count('input.key == "stk" or input.key == "pickrefs"'),
            2,
        )
        self.assertIn('function getDroppedJobArtifact(dataTransfer)', jobbuilder)
        self.assertIn('dataTransfer.getData("application/json")', jobbuilder)
        self.assertIn('dataTransfer.getData("text/plain")', jobbuilder)
        self.assertIn('input.dataset.artifactDropType === artifact.dataType', jobbuilder)
        self.assertIn('applyDroppedFilePath(input, artifact.path);', jobbuilder)

    def test_batch_detail_template_has_common_result_and_log_panels(self):
        batch_view = self._read_template("nice_classic/batchview.html")

        rendered = render_to_string("nice_classic/batchview.html", {
            "jobid": 7,
            "disp": 1,
            "name": "Import Movie Data",
            "desc": "",
            "status": "finished",
            "created": "today",
            "package": "simple",
            "program": "import_movies",
            "project": "project",
            "workspace": "workspace",
            "job_dir": "/project/workspace/1_import_movies",
            "source_label": "workspace project",
            "arguments": [
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
            ],
            "submitted_argument_count": 1,
            "result_project": None,
            "project_sections": [],
            "project_summary_available": False,
            "logs": [
                {
                    "name": "stdout.log",
                    "label": "standard output",
                    "exists": True,
                    "size": 12,
                    "truncated": False,
                    "text": "output",
                },
                {
                    "name": "stderr.log",
                    "label": "standard error",
                    "exists": False,
                    "size": 0,
                    "truncated": False,
                    "text": "",
                },
                {
                    "name": "nice_status.log",
                    "label": "NICE status callbacks",
                    "exists": True,
                    "size": 8,
                    "truncated": False,
                    "text": "running",
                },
            ],
            "artifact_counts": [],
            "artifact_images": [],
            "auto_refresh": False,
        })

        self.assertIn("job overview", batch_view)
        self.assertIn("submitted arguments", batch_view)
        self.assertIn("project results", batch_view)
        self.assertIn("<span>output</span>", batch_view)
        self.assertIn("{% for log in logs %}", batch_view)
        self.assertIn("{% if auto_refresh %}", batch_view)
        self.assertIn("Import Movie Data", rendered)
        self.assertIn("Number of threads", rendered)
        self.assertNotIn(">nthr</dt>", rendered)
        self.assertIn('id="argument_view_toggle"', rendered)
        self.assertIn('type="checkbox" role="switch"', rendered)
        self.assertIn('aria-label="show all arguments"', rendered)
        self.assertIn('data-argument-view-label="all"', rendered)
        self.assertNotIn('data-argument-view-label="submitted"', rendered)
        self.assertLess(
            rendered.index('data-argument-view-label="all"'),
            rendered.index('id="argument_view_toggle"'),
        )
        self.assertNotIn('id="submitted_arguments_button"', rendered)
        self.assertNotIn('id="all_arguments_button"', rendered)
        self.assertIn('batchArgumentViewToggle.addEventListener("change"', rendered)
        self.assertIn('batchArgumentViewToggle.checked ? "all" : "submitted"', rendered)
        self.assertIn('sessionStorage.getItem(batchArgumentViewStorageKey)', rendered)
        self.assertIn('data-submitted="false"', rendered)
        self.assertIn("Scale", rendered)
        self.assertIn("Mask diameter", rendered)
        self.assertIn("not set", rendered)
        self.assertIn("standard output", rendered)
        self.assertIn("standard error", rendered)
        self.assertIn("NICE status callbacks", rendered)

    def test_active_batch_detail_has_opt_in_molstar_volume_viewer(self):
        batch_view = self._read_template("nice_batch/batchview.html")
        volume_viewer = self._read_template("includes/_volume_viewer.html")
        viewer_3d = self._read_template("includes/_3D_viewer.html")
        context = {
            "jobid": 7,
            "disp": 9,
            "name": "Initial 3D Reconstruction",
            "desc": "",
            "proj": "project",
            "dset": "workspace",
            "args": {},
            "created": "today",
            "folder": "/workspace/9_abinitio3D",
            "jobstats": {},
            "log": [],
            "error": None,
            "arguments": [],
            "submitted_argument_count": 0,
            "volume_viewer_requested": True,
            "volume_outputs": [{
                "name": "recvol_state01.mrc",
                "state": 1,
                "width": 256,
                "height": 256,
                "depth": 256,
                "voxel_size": (1.3, 1.3, 1.3),
                "minimum": -2.0,
                "maximum": 8.0,
            }],
        }

        rendered = render_to_string("nice_batch/batchview.html", context)

        self.assertIn("{% include 'includes/_volume_viewer.html'", batch_view)
        self.assertIn('data-panel-target="volume3D"', rendered)
        self.assertIn("<span>3D volume</span>", rendered)
        self.assertIn('data-panel="volume3D"', rendered)
        self.assertIn('id="batch_volume_viewer"', rendered)
        self.assertIn("{% include 'includes/_3D_viewer.html' %}", volume_viewer)
        self.assertIn("{% include 'includes/_slider.html'", viewer_3d)
        self.assertIn("data-basic-3d-viewer", viewer_3d)
        self.assertIn("data-volume-molstar", rendered)
        self.assertIn("data-volume-isovalue", rendered)
        self.assertIn("data-volume-isovalue-text", rendered)
        self.assertIn('value="/batchvolume/7/recvol_state01.mrc"', rendered)
        self.assertIn("molstar@5.11.0/build/viewer/molstar.css", rendered)
        self.assertIn("molstar@5.11.0/build/viewer/molstar.js", rendered)
        self.assertIn("nice_lite/molstar_volume_viewer.js", rendered)
        self.assertNotIn("molstar_volume_viewer.js?v=", rendered)

        context["volume_viewer_requested"] = False
        context["volume_outputs"] = []
        default_off = render_to_string("nice_batch/batchview.html", context)

        self.assertNotIn('data-panel-target="volume3D"', default_off)
        self.assertNotIn('id="batch_volume_viewer"', default_off)
        self.assertNotIn("molstar@5.11.0", default_off)
