import json
from types import SimpleNamespace
from unittest.mock import patch

from django.test import SimpleTestCase

from ..data_structures import simple as simple_module
from ..data_structures.simple import SIMPLEProjFile


class SIMPLEProjFileTests(SimpleTestCase):
    def test_project_json_keeps_first_duplicate_character_value(self):
        output = json.dumps({"data": []})
        output = output.replace(
            "[]",
            '[{"movie":"/data/movie.mrcs","movie":0.0}]',
        )
        completed = SimpleNamespace(stdout=output)

        with patch.object(simple_module.subprocess, "run", return_value=completed):
            result = SIMPLEProjFile("/project/workspace.simple")._run(["simple_exec"])

        self.assertEqual(result["data"][0]["movie"], "/data/movie.mrcs")
