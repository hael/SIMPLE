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
