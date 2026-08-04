from pathlib import Path

from django.test import SimpleTestCase


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
