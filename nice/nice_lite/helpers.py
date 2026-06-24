"""Shared helper utilities for NICE Lite views and API handlers.

This module centralizes small, reusable helpers for:
- relaxed request value parsing and validation
- filesystem setup helpers used by project/workspace operations
- cookie and heartbeat convenience operations
"""

# global imports
import os

# django imports
from django.http import HttpResponse


# ------------------------------------------------------------------
# HTTP Helpers
# ------------------------------------------------------------------

class HttpResponseNoContent(HttpResponse):
    """HTTP 204 response that intentionally ignores content operations."""

    def __init__(self, content="", mimetype=None, status=None, content_type=None):
        del content, mimetype, status, content_type
        super().__init__(status=204)

        if "content-type" in self.headers:
            del self.headers["content-type"]

    def _set_content(self, value):
        """Ignore body assignment for explicit no-content responses."""
        del value

    def _get_content(self):
        """Always expose empty bytes for no-content responses."""
        return b""


# ------------------------------------------------------------------
# Value Presence Validators
# ------------------------------------------------------------------

def int_present(json, key, silent=False):
    """Return True when key exists and its value is an integer."""
    if key not in json:
        if not silent:
            print_error(key + " missing from json")
        return False
    if not isinstance(json[key], int):
        if not silent:
            print_error(key + " is not an integer")
        return False
    return True


def float_present(json, key, silent=False):
    """Return True when key exists and its value is a float."""
    if key not in json:
        if not silent:
            print_error(key + " missing from json")
        return False
    if not isinstance(json[key], float):
        if not silent:
            print_error(key + " is not a float")
        return False
    return True


def dict_present(json, key, silent=False):
    """Return True when key exists and its value is a dict."""
    if key not in json:
        if not silent:
            print_error(key + " missing from json")
        return False
    if not isinstance(json[key], dict):
        if not silent:
            print_error(key + " is not a dict")
        return False
    return True


def string_present(json, key, silent=False):
    """Return True when key exists and its value is a string."""
    if key not in json:
        if not silent:
            print_error(key + " missing from json")
        return False
    if not isinstance(json[key], str):
        if not silent:
            print_error(key + " is not a string")
        return False
    return True


# ------------------------------------------------------------------
# Type Parsing Helpers
# ------------------------------------------------------------------

def get_float(json, key, silent=False):
    """Return float value for key, parsing numeric strings when possible."""
    if key not in json:
        if not silent:
            print_error(key + " missing from json")
        return None

    if isinstance(json[key], float):
        return json[key]
    if isinstance(json[key], str):
        try:
            return float(json[key])
        except ValueError:
            if not silent:
                print_error(key + " is not an float string")
    return None


def get_integer(json, key, silent=False):
    """Return int value for key, parsing integer strings when possible."""
    if key not in json:
        if not silent:
            print_error(key + " missing from json")
        return None

    if isinstance(json[key], int):
        return json[key]
    if isinstance(json[key], str):
        try:
            return int(json[key])
        except ValueError:
            if not silent:
                print_error(key + " is not an integer string")
    return None


def get_string(json, key, silent=False):
    """Return string value for key when present and correctly typed."""
    if key not in json:
        if not silent:
            print_error(key + " missing from json")
        return None
    if isinstance(json[key], str):
        return json[key]
    return None


def print_error(error):
    """Emit a standardized helper error message to stdout."""
    print("ERROR: " + error)


# ------------------------------------------------------------------
# Filesystem Helpers
# ------------------------------------------------------------------

def directory_exists(dirc):
    """Return True when directory path exists and is a directory."""
    if not os.path.exists(dirc):
        return False
    if not os.path.isdir(dirc):
        return False
    return True


def ensure_directory(dirc):
    """Create directory path recursively when missing."""
    try:
        os.makedirs(dirc, exist_ok=True)
    except OSError:
        print_error("Directory " + dirc + " can not be created")
        return False
    return True


def create_symlink(src, dest):
    """Create symlink from src to dest, returning False on handled errors."""
    try:
        os.symlink(src, dest)
    except FileExistsError:
        print_error("Symlink already exists.")
        return False
    except PermissionError:
        print_error("Permission denied: You might need admin rights.")
        return False
    except OSError as error:
        print_error("OS error occurred: " + str(error))
        return False
    return True


# ------------------------------------------------------------------
# Request Selection Helpers
# ------------------------------------------------------------------

def get_project_id(request):
    """Resolve selected project id by GET, then POST, then cookie fallback."""
    selected_project_id_get = get_integer(request.GET, "selected_project_id", silent=True)
    selected_project_id_post = get_integer(request.POST, "selected_project_id", silent=True)
    selected_project_id_cookie = get_integer(request.COOKIES, "selected_project_id", silent=True)

    if selected_project_id_get is not None:
        return selected_project_id_get
    if selected_project_id_post is not None:
        return selected_project_id_post
    if selected_project_id_cookie is not None:
        return selected_project_id_cookie
    return None


def get_workspace_id(request):
    """Resolve selected workspace id by GET, then POST, then cookie fallback."""
    selected_workspace_id_get = get_integer(request.GET, "selected_workspace_id", silent=True)
    selected_workspace_id_post = get_integer(request.POST, "selected_workspace_id", silent=True)
    selected_workspace_id_cookie = get_integer(request.COOKIES, "selected_workspace_id", silent=True)

    if selected_workspace_id_get is not None:
        return selected_workspace_id_get
    if selected_workspace_id_post is not None:
        return selected_workspace_id_post
    if selected_workspace_id_cookie is not None:
        return selected_workspace_id_cookie
    return None


def get_job_id(request):
    """Resolve selected job id by GET, then POST."""
    selected_job_id_get = get_integer(request.GET, "selected_job_id", silent=True)
    selected_job_id_post = get_integer(request.POST, "selected_job_id", silent=True)

    if selected_job_id_get is not None:
        return selected_job_id_get
    if selected_job_id_post is not None:
        return selected_job_id_post
    return None


# ------------------------------------------------------------------
# Cookie and Heartbeat Helpers
# ------------------------------------------------------------------

def clear_checksum_cookies(request, response):
    """Remove checksum cookies used by iframe payload cache checks."""
    for cookie in request.COOKIES:
        if "checksum" in cookie:
            response.delete_cookie(key=cookie)


def analyse_heartbeat(heartbeat):
    """Return (status, timestamp) from heartbeat payload with defaults."""
    status = "unknown"
    timestamp = None
    if "status" in heartbeat:
        status = heartbeat["status"]
    if "timestamp" in heartbeat:
        timestamp = heartbeat["timestamp"]
    return status, timestamp
