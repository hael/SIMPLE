from django.core.exceptions import ImproperlyConfigured
from django.test import SimpleTestCase, override_settings

from ..features import batch_job_controls_enabled, workspace_job_refresh_enabled


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
