from pathlib import Path

from django.test import SimpleTestCase
from django.template.loader import render_to_string


class TemplateIntegrationTests(SimpleTestCase):
    @staticmethod
    def _read_template(relative_path):
        base_dir = Path(__file__).resolve().parents[1]
        template_path = base_dir / "templates" / relative_path
        return template_path.read_text(encoding="utf-8")

    def test_zoom_log_template_uses_dataset_not_workspace_property(self):
        content = self._read_template("nice_stream/includes/_zoom_log_errors_parts_body.html")
        self.assertNotIn("panel.workspace", content)
        self.assertIn("panel.dataset.expanded", content)

    def test_workspace_reload_guard_skips_when_job_builder_visible(self):
        content = self._read_template("workspace.html")
        self.assertIn("isJobBuilderVisible", content)
        self.assertIn("&& !isJobBuilderVisible()", content)

    def test_browser_messages_remain_visible(self):
        content = self._read_template("messages.html")

        self.assertIn('id="message_alert"', content)
        self.assertNotIn("setTimeout", content)
        self.assertNotIn('style.display="none"', content)

    def test_jobbuilder_submits_to_named_workspace_iframe(self):
        jobbuilder = self._read_template("jobbuilder.html")
        index = self._read_template("index.html")
        self.assertIn('target="workspace_iframe"', jobbuilder)
        self.assertIn('name="workspace_iframe"', index)

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

    def test_file_browser_navigation_passes_paths_as_query_parameters(self):
        filebrowser = self._read_template("filebrowser.html")

        self.assertNotIn("{% url 'nice_lite:file_browser' type parentdir %}", filebrowser)
        self.assertNotIn("{% url 'nice_lite:file_browser' type path|add:'/'|add:dir %}", filebrowser)
        self.assertIn('name="selectedpath" value="{{ parentdir }}"', filebrowser)

    def test_batch_file_browser_restore_does_not_show_loading_overlay(self):
        filebrowser = self._read_template("filebrowser.html")

        self.assertIn("const restoreJobsIframe = (iframe) => {", filebrowser)
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

    def test_batch_submission_guard_does_not_change_stream_submission(self):
        jobbuilder = self._read_template("jobbuilder.html")

        self.assertIn("action=\"{% url 'nice_lite:create_stream' %}\"", jobbuilder)
        self.assertIn('<button type="submit"', jobbuilder)
        self.assertEqual(jobbuilder.count('data-batch-package='), 2)
        self.assertIn('document.querySelectorAll("form[data-batch-package]")', jobbuilder)

    def test_batch_project_source_defaults_to_workspace_and_lists_snapshots(self):
        rendered = render_to_string("jobbuilder.html", {
            "stream_user_inputs": [],
            "simple_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "simple_program_inputs": [{"prg": "demo", "disp": "Demo", "sections": []}],
            "single_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "single_program_inputs": [{"prg": "demo", "disp": "Demo", "sections": []}],
            "batch_snapshot_sources": [{
                "key": "snapshot:12:3",
                "label": "stream 2 - particle set 3",
            }],
        })

        self.assertEqual(rendered.count('value="workspace" selected>workspace project'), 2)
        self.assertEqual(rendered.count('value="snapshot:12:3"'), 2)
        self.assertEqual(rendered.count('name="batch_source" value="workspace"'), 2)
        self.assertIn('id="simple_batch_source"', rendered)
        self.assertIn('id="single_batch_source"', rendered)

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
            "batch_snapshot_sources": [],
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
