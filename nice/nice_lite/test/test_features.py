from django.core.exceptions import ImproperlyConfigured
from django.test import SimpleTestCase, override_settings

from ..features import (
    batch_job_controls_enabled,
    batch_project_file_selector_enabled,
    batch_project_inheritance_enabled,
    batch_status_callbacks_enabled,
    workspace_job_refresh_enabled,
)


class FeatureFlagTests(SimpleTestCase):
    @override_settings(NICE_LITE_BATCH_JOB_CONTROLS=False)
    def test_batch_job_controls_default_off_path(self):
        self.assertFalse(batch_job_controls_enabled())

    @override_settings(NICE_LITE_BATCH_JOB_CONTROLS=True)
    def test_batch_job_controls_opt_in_path(self):
        self.assertTrue(batch_job_controls_enabled())

    @override_settings(NICE_LITE_BATCH_JOB_CONTROLS="true")
    def test_batch_job_controls_reject_non_boolean_setting(self):
        with self.assertRaises(ImproperlyConfigured):
            batch_job_controls_enabled()

    @override_settings(NICE_LITE_BATCH_PROJECT_INHERITANCE=False)
    def test_batch_project_inheritance_default_off_path(self):
        self.assertFalse(batch_project_inheritance_enabled())

    @override_settings(NICE_LITE_BATCH_PROJECT_INHERITANCE=True)
    def test_batch_project_inheritance_opt_in_path(self):
        self.assertTrue(batch_project_inheritance_enabled())

    @override_settings(NICE_LITE_BATCH_PROJECT_INHERITANCE="true")
    def test_batch_project_inheritance_reject_non_boolean_setting(self):
        with self.assertRaises(ImproperlyConfigured):
            batch_project_inheritance_enabled()

    @override_settings(NICE_LITE_BATCH_PROJECT_FILE_SELECTOR=False)
    def test_batch_project_file_selector_default_off_path(self):
        self.assertFalse(batch_project_file_selector_enabled())

    @override_settings(NICE_LITE_BATCH_PROJECT_FILE_SELECTOR=True)
    def test_batch_project_file_selector_opt_in_path(self):
        self.assertTrue(batch_project_file_selector_enabled())

    @override_settings(NICE_LITE_BATCH_PROJECT_FILE_SELECTOR="true")
    def test_batch_project_file_selector_reject_non_boolean_setting(self):
        with self.assertRaises(ImproperlyConfigured):
            batch_project_file_selector_enabled()

    @override_settings(NICE_LITE_BATCH_STATUS_CALLBACKS=False)
    def test_batch_status_callbacks_default_off_path(self):
        self.assertFalse(batch_status_callbacks_enabled())

    @override_settings(NICE_LITE_BATCH_STATUS_CALLBACKS=True)
    def test_batch_status_callbacks_opt_in_path(self):
        self.assertTrue(batch_status_callbacks_enabled())

    @override_settings(NICE_LITE_BATCH_STATUS_CALLBACKS="true")
    def test_batch_status_callbacks_reject_non_boolean_setting(self):
        with self.assertRaises(ImproperlyConfigured):
            batch_status_callbacks_enabled()

    @override_settings(NICE_LITE_WORKSPACE_JOB_REFRESH=False)
    def test_workspace_job_refresh_default_off_path(self):
        self.assertFalse(workspace_job_refresh_enabled())

    @override_settings(NICE_LITE_WORKSPACE_JOB_REFRESH=True)
    def test_workspace_job_refresh_opt_in_path(self):
        self.assertTrue(workspace_job_refresh_enabled())

    @override_settings(NICE_LITE_WORKSPACE_JOB_REFRESH="true")
    def test_workspace_job_refresh_rejects_non_boolean_setting(self):
        with self.assertRaises(ImproperlyConfigured):
            workspace_job_refresh_enabled()
