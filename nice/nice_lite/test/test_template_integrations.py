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

    def test_workspace_refresh_button_posts_missing_directory_cleanup(self):
        content = self._read_template("workspace.html")

        self.assertIn("workspace_job_refresh_enabled", content)
        self.assertIn('<form method="POST" action="{% url \'nice_lite:refresh_workspace_jobs\' %}">', content)
        self.assertIn('title="refresh job cards"', content)
        self.assertIn("{% csrf_token %}", content)

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
            "batch_project_sources": [{
                "key": "snapshot:12:3",
                "label": "stream 2 - particle set 3",
            }],
            "default_batch_source": "workspace",
        })

        self.assertEqual(rendered.count('value="workspace" selected>workspace project'), 2)
        self.assertEqual(rendered.count('value="snapshot:12:3"'), 2)
        self.assertEqual(rendered.count('name="batch_source" value="workspace"'), 2)
        self.assertIn('id="simple_batch_source"', rendered)
        self.assertIn('id="single_batch_source"', rendered)

    def test_batch_project_source_defaults_to_latest_completed_job(self):
        rendered = render_to_string("jobbuilder.html", {
            "stream_user_inputs": [],
            "simple_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "simple_program_inputs": [{"prg": "demo", "disp": "Demo", "sections": []}],
            "single_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "single_program_inputs": [{"prg": "demo", "disp": "Demo", "sections": []}],
            "batch_project_sources": [{
                "key": "job:8",
                "label": "job 1 - Import Movie Data",
            }],
            "default_batch_source": "job:8",
        })

        self.assertEqual(rendered.count('value="job:8" selected'), 2)
        self.assertEqual(rendered.count('name="batch_source" value="job:8"'), 2)
        self.assertEqual(rendered.count('value="workspace" >workspace project'), 2)

    def test_batch_project_file_selector_uses_inherited_simple_project(self):
        project_path = "/workspace/1_import_movies/workspace.simple"
        rendered = render_to_string("jobbuilder.html", {
            "stream_user_inputs": [],
            "simple_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "simple_program_inputs": [{"prg": "demo", "disp": "Demo", "sections": []}],
            "single_programs": [{"prg": "demo", "disp": "Demo", "desc": ""}],
            "single_program_inputs": [{"prg": "demo", "disp": "Demo", "sections": []}],
            "batch_project_file_selector_enabled": True,
            "default_batch_project_file": project_path,
        })

        self.assertIn('id="simple_batch_project_file"', rendered)
        self.assertIn('id="single_batch_project_file"', rendered)
        self.assertNotIn('id="simple_batch_source"', rendered)
        self.assertNotIn('id="single_batch_source"', rendered)
        self.assertEqual(
            rendered.count(f'name="batch_project_file" value="{project_path}"'),
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
        self.assertIn('for="batch_simple_demo_required_arg">required simple', rendered)
        self.assertIn('for="single_single_demo_required_arg">required single', rendered)
        self.assertIn(
            'class="text-xs font-medium text-streamtext whitespace-nowrap" '
            'for="batch_simple_demo_optional_arg">optional simple',
            rendered,
        )
        self.assertIn(
            'class="text-xs font-medium text-streamtext whitespace-nowrap" '
            'for="single_single_demo_optional_arg">optional single',
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
        jobs = self._read_template("jobs.html")

        self.assertIn("batch_job_controls_enabled", batch_card)
        self.assertIn("{% url 'nice_lite:stop_batch' %}", batch_card)
        self.assertIn('onclick="stopBatchJob(this)"', batch_card)
        self.assertIn('<circle cx="8" cy="8" r="6"/>', batch_card)
        self.assertIn('<rect x="5.5" y="5.5" width="5" height="5"', batch_card)
        self.assertIn("{% url 'nice_lite:delete_batch' %}", batch_card)
        self.assertIn('onclick="deleteBatchJob(this)"', batch_card)
        self.assertIn('<polyline points="2,4 14,4"/>', batch_card)
        self.assertIn('<path d="M5 4V3a1 1 0 0 1 1-1h4a1 1 0 0 1 1 1v1M6 7v5M10 7v5M3 4l1 9a1 1 0 0 0 1 1h6a1 1 0 0 0 1-1l1-9"/>', batch_card)
        self.assertIn("const stopBatchJob = (button) => {", jobs)
        self.assertIn("const deleteBatchJob = (button) => {", jobs)
        self.assertIn("permanently delete batch job", jobs)
        self.assertIn("This cannot be undone.", jobs)

    def test_batch_card_controls_follow_opt_in_and_job_status(self):
        job = {
            "id": 7,
            "disp": 1,
            "name": "import movies",
            "dirc": "1_import_movies",
            "args": {},
            "master_stats": {"package": "simple"},
            "status": "running",
        }

        disabled = render_to_string("nice_classic/_batch_card.html", {
            "job": job,
            "batch_job_controls_enabled": False,
        })
        running = render_to_string("nice_classic/_batch_card.html", {
            "job": job,
            "batch_job_controls_enabled": True,
        })
        job["status"] = "queued"
        queued = render_to_string("nice_classic/_batch_card.html", {
            "job": job,
            "batch_job_controls_enabled": True,
        })
        job["status"] = "finished"
        finished = render_to_string("nice_classic/_batch_card.html", {
            "job": job,
            "batch_job_controls_enabled": True,
        })

        self.assertNotIn('title="stop batch job"', disabled)
        self.assertNotIn('title="delete batch job"', disabled)
        self.assertIn('title="stop batch job"', running)
        self.assertNotIn('title="delete batch job"', running)
        self.assertIn("rounded-full bg-streamaction/10 text-streamaction", running)
        self.assertNotIn("w-1.5 h-1.5", running)
        self.assertLess(
            running.index('title="stop batch job"'),
            running.index('job directory'),
        )
        self.assertNotIn('title="stop batch job"', queued)
        self.assertIn('title="delete batch job"', queued)
        self.assertIn("rounded-full bg-streamline/40 text-streamlabel", queued)
        self.assertNotIn('title="stop batch job"', finished)
        self.assertIn('title="delete batch job"', finished)
        self.assertIn("rounded-full bg-streamring/10 text-streamaccent", finished)
        self.assertLess(
            finished.index('title="delete batch job"'),
            finished.index('job directory'),
        )
