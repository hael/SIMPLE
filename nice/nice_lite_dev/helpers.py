# global imports
import os
from django.http import HttpResponse

class HttpResponseNoContent(HttpResponse):
    """Special HTTP response with no content, just headers.

    The content operations are ignored.
    """

    def __init__(self, content="", mimetype=None, status=None, content_type=None):
        super().__init__(status=204)

        if "content-type" in self.headers:
            del self.headers["content-type"]

    def _set_content(self, value):
        pass

    def _get_content(self, value):
        pass

def int_present( json, key, silent=False ):
  if key not in json:
    if not silent:
      print_error(key + " missing from json")
    return False
  if not isinstance(json[key], int):
    if not silent:
      print_error(key + " is not an integer")
    return False
  return True

def float_present( json, key, silent=False ):
  if key not in json:
    if not silent:
      print_error(key + " missing from json")
    return False
  if not isinstance(json[key], float):
    if not silent:
      print_error(key + " is not a float")
    return False
  return True

def dict_present( json, key, silent=False ):
  if key not in json:
    if not silent:
      print_error(key + " missing from json")
    return False
  if not isinstance(json[key], dict):
    if not silent:
      print_error(key + " is not a dict")
    return False
  return True

def string_present( json, key, silent=False ):
  if key not in json:
    if not silent:
      print_error(key + " missing from json")
    return False
  if not isinstance(json[key], str):
    if not silent:
      print_error(key + " is not a string")
    return False
  return True

def get_float( json, key, silent=False ):
  if key not in json:
    if not silent:
      print_error(key + " missing from json")
    return None
  if isinstance(json[key], float):
    return json[key]
  elif isinstance(json[key], str):
    try:
      num = float(json[key])
      return num
    except ValueError:
      if not silent:
        print_error(key + " is not an float string")
  return None

def get_integer( json, key, silent=False ):
  if key not in json:
    if not silent:
      print_error(key + " missing from json")
    return None
  if isinstance(json[key], int):
    return json[key]
  elif isinstance(json[key], str):
    try:
      num = int(json[key])
      return num
    except ValueError:
      if not silent:
        print_error(key + " is not an integer string")
  return None

def get_string( json, key, silent=False ):
  if key not in json:
    if not silent:
      print_error(key + " missing from json")
    return None
  if isinstance(json[key], str):
    return json[key]
  return None

def print_error( error ):
  print("ERROR: " + error)

def directory_exists( dirc ):
  if not os.path.exists(dirc):
    return False
  if not os.path.isdir(dirc):
    return False
  return True

def ensure_directory( dirc ):
  try:
    os.makedirs(dirc, exist_ok=True)
  except OSError as error:
    print_error("Directory " + dirc + " can not be created")
    return False
  return True

def create_symlink( src, dest ):
  try:
      os.symlink(src, dest)
  except FileExistsError:
      print_error("Symlink already exists.")
      return False
  except PermissionError:
      print_error("Permission denied: You might need admin rights.")
      return False
  except OSError as e:
      print_error("OS error occurred:", e)
      return False
  return True

def get_project_id( request ):
  selected_project_id_get    = get_integer(request.GET,     "selected_project_id", silent=True)
  selected_project_id_post   = get_integer(request.POST,    "selected_project_id", silent=True)
  selected_project_id_cookie = get_integer(request.COOKIES, "selected_project_id", silent=True)
  if selected_project_id_get is not None:
    return selected_project_id_get
  elif selected_project_id_post is not None:
    return selected_project_id_post
  elif selected_project_id_cookie is not None:
    return selected_project_id_cookie
  else:
    return None

def get_dataset_id( request ):
  selected_dataset_id_get    = get_integer(request.GET,     "selected_dataset_id", silent=True)
  selected_dataset_id_post   = get_integer(request.POST,    "selected_dataset_id", silent=True)
  selected_dataset_id_cookie = get_integer(request.COOKIES, "selected_dataset_id", silent=True)
  if selected_dataset_id_get is not None:
    return selected_dataset_id_get
  elif selected_dataset_id_post is not None:
    return selected_dataset_id_post  
  elif selected_dataset_id_cookie is not None:
    return selected_dataset_id_cookie
  else:
    return None

def get_job_id( request ):
  selected_job_id_get  = get_integer(request.GET,  "selected_job_id", silent=True)
  selected_job_id_post = get_integer(request.POST, "selected_job_id", silent=True)
  if selected_job_id_get is not None:
    return selected_job_id_get
  elif selected_job_id_post is not None:
    return selected_job_id_post  
  else:
    return None

def clear_checksum_cookies( request, response ):
  for cookie in request.COOKIES:
    if "checksum" in cookie:
      response.delete_cookie(key=cookie)

def analyse_heartbeat( heartbeat ):
  status    = 'unknown'
  timestamp = None
  if "status" in heartbeat:
    status = heartbeat["status"]
  if "timestamp" in heartbeat:
    timestamp = heartbeat["timestamp"]
  return status, timestamp
