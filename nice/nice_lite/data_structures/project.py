# global imports
import os

from django.utils import timezone

# local imports
from ..helpers import directory_exists, ensure_directory
from ..models import ProjectModel, WorkspaceModel


class Project:
    """Top-level project directory grouping one or more workspaces."""

    def __init__(self, id=None, project_id=None, request=None):
        self.id = 0
        self.projectmodel = None
        self.absdir = None
        self.trashdir = None
        self.name = ""
        self.desc = ""
        self.dirc = ""
        self.date = None
        self.workspaces_list = []

        if project_id is not None:
            self.id = project_id
            self.load()
        elif id is not None:
            self.id = id
            self.load()
        elif request is not None:
            self.setIDFromRequest(request)
            if self.id > 0:
                self.load()

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def setIDFromRequest(self, request):
        """Populate self.id from request data/cookies when available."""
        project_id_str = None
        if "selected_project_id" in request.POST:
            project_id_str = request.POST["selected_project_id"]
        elif "selected_project_id" in request.GET:
            project_id_str = request.GET["selected_project_id"]
        else:
            project_id_str = request.COOKIES.get("selected_project_id", "none")

        if project_id_str is not None and str(project_id_str).isnumeric():
            self.id = int(project_id_str)

    def load(self):
        """Populate fields from database state, reset if record is missing."""
        self.projectmodel = ProjectModel.objects.filter(id=self.id).first()
        if self.projectmodel is None:
            self.id = 0
            self.absdir = None
            self.trashdir = None
            self.name = ""
            self.desc = ""
            self.dirc = ""
            self.date = None
            self.workspaces_list = []
            return

        self.absdir = self.projectmodel.dirc
        self.trashdir = os.path.join(self.projectmodel.dirc, "TRASH")
        self.name = self.projectmodel.name
        self.desc = self.projectmodel.desc
        self.dirc = self.projectmodel.dirc
        self.date = self.projectmodel.date
        self.workspaces_list = WorkspaceModel.objects.filter(proj=self.projectmodel)

    # ------------------------------------------------------------------
    # Mutators
    # ------------------------------------------------------------------

    def new(self, name, dirc=None):
        """
        Create a new project directory and DB record.

        Supports both signatures for compatibility:
        - new(name, dirc)
        - new(request) where POST has new_project_name/new_project_dirc
        """
        if dirc is None and hasattr(name, "POST"):
            request = name
            name = request.POST.get("new_project_name")
            dirc = request.POST.get("new_project_dirc")

        if name is None or dirc is None:
            return False
        if not directory_exists(dirc):
            return False

        new_project_path = os.path.join(dirc, str(name).replace(" ", "_"))
        if not ensure_directory(new_project_path):
            return False

        projectmodel = ProjectModel(name=name, dirc=new_project_path, date=timezone.now())
        projectmodel.save()
        self.id = projectmodel.id
        if self.id == 0:
            return False

        self.load()
        return True

    def ensureTrashfolder(self):
        """Ensure the project's TRASH folder exists on disk."""
        if self.absdir is None:
            return False
        if not os.path.isdir(self.absdir):
            return False
        if os.path.isdir(self.trashdir):
            return True
        try:
            os.makedirs(self.trashdir, exist_ok=True)
        except OSError:
            return False
        return True

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def get_projectmodel(self):
        """Return the loaded ProjectModel instance (or None)."""
        return self.projectmodel

    def get_absdir(self):
        """Return absolute filesystem path for this project."""
        return self.absdir

    def get_trashdir(self):
        """Return absolute filesystem path for this project's TRASH dir."""
        return self.trashdir
