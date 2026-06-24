import json
import os
import tempfile
from unittest.mock import Mock
from unittest.mock import patch

from django.http import HttpResponse
from django.test import RequestFactory
from django.test import SimpleTestCase

from .. import api


class _AuthUser:
    is_authenticated = True
    username = "worker"


class ApiHelperTests(SimpleTestCase):
    def test_get_valid_job_id_rejects_bool(self):
        self.assertIsNone(api._get_valid_job_id({"jobid": True}))

    def test_get_valid_job_id_accepts_positive_int(self):
        self.assertEqual(api._get_valid_job_id({"jobid": 5}), 5)

    def test_has_numeric_version_accepts_int_or_float(self):
        self.assertTrue(api._has_numeric_version({"version": 1}))
        self.assertTrue(api._has_numeric_version({"version": 1.5}))

    def test_has_numeric_version_rejects_bool_and_string(self):
        self.assertFalse(api._has_numeric_version({"version": True}))
        self.assertFalse(api._has_numeric_version({"version": "1"}))

    @patch.dict(os.environ, {"NICE_LITE_WORKER_TOKEN": "secret"}, clear=False)
    def test_worker_auth_accepts_header_token(self):
        req = RequestFactory().post("/api", HTTP_X_WORKER_TOKEN="secret")
        self.assertTrue(api._is_worker_authorized(req))

    @patch.dict(os.environ, {"NICE_LITE_WORKER_TOKEN": "secret"}, clear=False)
    def test_worker_auth_accepts_bearer_token(self):
        req = RequestFactory().post("/api", HTTP_AUTHORIZATION="Bearer secret")
        self.assertTrue(api._is_worker_authorized(req))

    @patch.dict(os.environ, {"NICE_LITE_WORKER_TOKEN": "secret"}, clear=False)
    def test_worker_auth_rejects_missing_token(self):
        req = RequestFactory().post("/api")
        self.assertFalse(api._is_worker_authorized(req))

    @patch.dict(os.environ, {}, clear=False)
    def test_worker_auth_is_opt_in_when_env_missing(self):
        original = os.environ.pop("NICE_LITE_WORKER_TOKEN", None)
        try:
            req = RequestFactory().post("/api")
            self.assertTrue(api._is_worker_authorized(req))
        finally:
            if original is not None:
                os.environ["NICE_LITE_WORKER_TOKEN"] = original

    def test_image_content_type_mapping(self):
        self.assertEqual(api._image_content_type_from_path("a.jpg"), "image/jpeg")
        self.assertEqual(api._image_content_type_from_path("a.png"), "image/png")
        self.assertEqual(api._image_content_type_from_path("a.gif"), "image/gif")
        self.assertEqual(api._image_content_type_from_path("a.webp"), "image/webp")


class ApiEndpointTests(SimpleTestCase):
    def setUp(self):
        self.factory = RequestFactory()

    def _post_json(self, path, payload, **headers):
        return self.factory.post(path, data=json.dumps(payload), content_type="application/json", **headers)

    @patch.dict(os.environ, {"NICE_LITE_WORKER_TOKEN": "secret"}, clear=False)
    def test_index_rejects_unauthorized_worker(self):
        response = api.index(self._post_json("/api", {"version": 1, "jobid": 1, "stream_heartbeat": {}}))
        self.assertEqual(response.status_code, 403)

    def test_index_rejects_non_dict_payload(self):
        response = api.index(self._post_json("/api", [1, 2, 3]))
        self.assertEqual(response.status_code, 400)

    def test_index_rejects_bool_jobid(self):
        payload = {"version": 1, "jobid": True, "stream_heartbeat": {}}
        response = api.index(self._post_json("/api", payload))
        self.assertEqual(response.status_code, 400)

    def test_index_returns_404_for_unknown_job(self):
        payload = {"version": 1, "jobid": 77, "stream_heartbeat": {}}
        with patch.object(api, "StreamJob") as mock_stream_job:
            mock_stream_job.return_value.get_jobmodel.return_value = None
            response = api.index(self._post_json("/api", payload))
        self.assertEqual(response.status_code, 404)

    def test_index_returns_master_update_on_valid_stream_heartbeat(self):
        payload = {"version": 1, "jobid": 9, "stream_heartbeat": {"master": {}}}
        stream_job = Mock()
        stream_job.get_jobmodel.return_value = object()
        stream_job.update_stats.return_value = True
        stream_job.get_master_update.return_value = {"ok": True}

        with patch.object(api, "StreamJob", return_value=stream_job):
            response = api.index(self._post_json("/api", payload))

        self.assertEqual(response.status_code, 200)
        self.assertJSONEqual(response.content.decode("utf-8"), {"ok": True})

    def test_index_returns_400_for_unknown_heartbeat_type(self):
        payload = {"version": 1, "jobid": 9}
        stream_job = Mock()
        stream_job.get_jobmodel.return_value = object()
        with patch.object(api, "StreamJob", return_value=stream_job):
            response = api.index(self._post_json("/api", payload))
        self.assertEqual(response.status_code, 400)

    def test_image_returns_404_when_path_invalid(self):
        request = self.factory.get("/image")
        request.user = _AuthUser()
        with patch.object(api, "_resolve_safe_image_path", return_value=None):
            response = api.image(request, "x.jpg")
        self.assertEqual(response.status_code, 404)

    def test_image_uses_extension_content_type(self):
        request = self.factory.get("/image")
        request.user = _AuthUser()

        with tempfile.NamedTemporaryFile(suffix=".png") as image_file:
            image_file.write(b"PNGDATA")
            image_file.flush()

            with patch.object(api, "_resolve_safe_image_path", return_value=image_file.name):
                response = api.image(request, "x.png")

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "image/png")
