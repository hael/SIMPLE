# global imports
import os
from django.utils import timezone

# local imports
from ..helpers import *
from ..models  import ProjectModel


class Project:
    """
    Represents a top-level project directory that groups one or more datasets.

    A project maps to a single directory on disk whose path is stored in the DB.
    A TRASH/ subdirectory is used for soft-deletes of child datasets and workspaces.

    Lifecycle:
      - Create : Project().new(name, dirc)
      - Load   : Project(id)
      - Delete : called implicitly by Dataset.delete() / Workspace.delete()
                 when the project becomes empty
    """

    def __init__(self, id=None):
        # unit tests:
        self.id           = 0
        self.projectmodel = None
        self.absdir       = None  # absolute path to the project directory
        self.trashdir     = None  # absolute path to the TRASH subdirectory
        if id is not None:
            self.id = id
            self.load()

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def load(self):
        """Populate fields from the database. Resets to empty state if not found."""
        self.projectmodel = ProjectModel.objects.filter(id=self.id).first()
        if self.projectmodel is None:
            self.id       = 0
            self.absdir   = None
            self.trashdir = None
        else:
            self.absdir   = self.projectmodel.dirc
            self.trashdir = os.path.join(self.projectmodel.dirc, "TRASH")

    # ------------------------------------------------------------------
    # Mutators
    # ------------------------------------------------------------------

    def new(self, name, dirc):
        """
        Create a new project directory inside dirc and save the DB record.
        The directory name is derived from name with spaces replaced by underscores.
        Returns True on success, False on any failure.

        unit tests:
        """
        if not directory_exists(dirc):
            return False
        new_project_path = os.path.join(dirc, name.replace(" ", "_"))
        if not ensure_directory(new_project_path):
            return False
        projectmodel = ProjectModel(name=name, dirc=new_project_path, date=timezone.now())
        projectmodel.save()
        self.id = projectmodel.id
        if self.id == 0:
            return False
        self.load()
        return True

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    # UNUSED
    # def get_id(self):
    #     return self.id

    def get_projectmodel(self):
        return self.projectmodel

    def get_absdir(self):
        return self.absdir

    def get_trashdir(self):
        return self.trashdir
