"""Workspace data structure: filesystem mapping and lifecycle helpers."""

# global imports
import os
import re

from django.utils import timezone

# local imports
from ..helpers import *
from ..models import WorkspaceModel


class Workspace:
    """
    Represents a cryo-EM workspace within a project.

    A workspace maps to:
      - a hidden directory  (.workspace_<id>/)  inside the project directory
      - a human-readable symlink (ds_<name>/) pointing to that directory
      - a TRASH/ subdirectory inside the workspace directory used for soft-deletes

    Lifecycle:
      - Create : Workspace().new(project, user)
      - Load   : Workspace(id)
      - Delete : workspace.delete()     — moves files to project TRASH, removes DB record
      - Rename : workspace.rename(name) — renames the symlink and updates the DB
    """

    def __init__(self, id=None):
        # unit tests: test_workspace_init, test_workspace_init_by_id
        self.id             = 0
        self.workspacemodel = None
        self.absdir         = None  # absolute path to the hidden workspace directory
        self.linkpath       = None  # absolute path to the human-readable symlink
        self.trashdir       = None  # absolute path to the TRASH subdirectory
        if id is not None:
            self.id = id
            self.load()

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def load(self):
        """Populate fields from the database. Resets to empty state if not found."""
        self.workspacemodel = WorkspaceModel.objects.filter(id=self.id).first()
        if self.workspacemodel is None:
            self.id       = 0
            self.absdir   = None
            self.linkpath = None
            self.trashdir = None
        else:
            self.absdir   = os.path.join(self.workspacemodel.proj.dirc, self.workspacemodel.dirc)
            self.linkpath = os.path.join(self.workspacemodel.proj.dirc, self.workspacemodel.link)
            self.trashdir = os.path.join(self.absdir, "TRASH")

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def in_project(self, projectid):
        """Return True if this workspace belongs to the given project id."""
        return self.workspacemodel is not None and self.workspacemodel.proj.id == projectid

    def get_id(self):
        """Return the workspace id."""
        return self.id

    def get_workspacemodel(self):
        if self.workspacemodel is None:
            print_error("workspacemodel is none")
            return None
        return self.workspacemodel

    def get_absdir(self):
        """Return absolute path to the hidden workspace directory."""
        return self.absdir

    def get_linkpath(self):
        """Return absolute path to the human-readable workspace symlink."""
        return self.linkpath

    def get_trashdir(self):
        """Return absolute path to the workspace TRASH directory."""
        return self.trashdir

    # ------------------------------------------------------------------
    # Mutators
    # ------------------------------------------------------------------

    def new(self, project, user):
        """
        Create a new workspace inside the given project.

        Creates the hidden directory and symlink on disk, then writes the DB record.
        The DB record is rolled back if filesystem operations fail.
        Returns True on success, False on any failure.

        unit tests: test_workspace_new
        """
        projectmodel = project.get_projectmodel()
        projectdir   = project.get_absdir()
        if projectmodel is None:
            print_error("projectmodel is none")
            return False
        if not directory_exists(projectdir):
            print_error("project directory " + projectdir + " doesn't exist")
            return False

        # derive a display index and default name from existing workspace count
        display_id = WorkspaceModel.objects.filter(proj=projectmodel).count() + 1
        new_workspace_name = "new workspace " + str(display_id)

        # save the record early to obtain the auto-assigned id, which is used for the directory name
        workspacemodel = WorkspaceModel(
            proj=projectmodel, name=new_workspace_name, disp=display_id,
            cdat=timezone.now(), mdat=timezone.now(), user=user
        )
        workspacemodel.save()
        self.id = workspacemodel.id
        if self.id == 0:
            return False

        new_workspace_dirc = ".workspace_" + str(self.id)
        new_workspace_link = "ds_" + new_workspace_name.replace(" ", "_")
        new_workspace_path = os.path.join(projectdir, new_workspace_dirc)

        if not ensure_directory(new_workspace_path):
            workspacemodel.delete()
            return False
        if not create_symlink(new_workspace_path, os.path.join(projectdir, new_workspace_link)):
            workspacemodel.delete()
            return False

        workspacemodel.dirc = new_workspace_dirc
        workspacemodel.link = new_workspace_link
        workspacemodel.save()
        self.load()
        return True

    def delete(self):
        """
        Soft-delete this workspace by moving its directory and symlink into the project TRASH.

        Also removes the DB record. If the project has no remaining workspaces
        after deletion, the project itself is also deleted.
        Returns True on success, False on any failure.

        unit tests: test_workspace_delete
        """
        if self.workspacemodel is None:
            print_error("workspacemodel is none")
            return False

        projecttrash = os.path.join(self.workspacemodel.proj.dirc, "TRASH")
        if not ensure_directory(projecttrash):
            print_error("project trash folder doesn't exist")
            return False

        trash_path = os.path.join(projecttrash, self.workspacemodel.dirc)
        trash_link = os.path.join(projecttrash, self.workspacemodel.link)

        # verify source filesystem state
        if not os.path.exists(self.absdir):
            print_error("workspace directory doesn't exist")
            return False
        if not os.path.isdir(self.absdir):
            print_error("workspace directory isn't a directory")
            return False
        if not os.path.exists(self.linkpath):
            print_error("workspace link doesn't exist")
            return False
        if not os.path.islink(self.linkpath):
            print_error("workspace link isn't a link")
            return False

        # verify trash destinations are free
        if os.path.exists(trash_path):
            print_error("workspace trash directory already exists")
            return False
        if os.path.exists(trash_link):
            print_error("trash link already exists")
            return False

        try:
            os.rename(self.absdir, trash_path)
        except OSError:
            print_error("directory cannot be renamed: " + self.absdir)
            return False

        try:
            os.remove(self.linkpath)
            os.symlink(trash_path, trash_link)
        except OSError:
            print_error("link cannot be updated: " + self.linkpath)
            return False

        # remove DB record; delete the parent project if it is now empty
        projectmodel = self.workspacemodel.proj
        self.workspacemodel.delete()
        remaining_workspaces = WorkspaceModel.objects.filter(proj=projectmodel).count()
        if remaining_workspaces == 0:
            projectmodel.delete()

        return True

    def rename(self, newworkspacename):
        """
        Rename this workspace: renames the symlink on disk and updates the DB record.
        The directory itself is not moved — only the symlink is renamed.
        Returns True on success, False on any failure.

        unit tests: test_workspace_rename
        """
        if self.workspacemodel is None:
            print_error("workspacemodel is none")
            return False
        if newworkspacename is None:
            print_error("new workspace name is none")
            return False

        # build a filesystem-safe link name: spaces become _, non-alphanumeric stripped
        new_link_name    = "ds_" + re.sub(r'\W+', '', newworkspacename.replace(" ", "_"))
        new_workspace_link = os.path.join(self.workspacemodel.proj.dirc, new_link_name)

        try:
            os.rename(self.linkpath, new_workspace_link)
        except IsADirectoryError:
            print_error("Source is a file but destination is a directory.")
            return False
        except NotADirectoryError:
            print_error("Source is a directory but destination is a file.")
            return False
        except PermissionError:
            print_error("Operation not permitted")
            return False
        except OSError as error:
            print_error(error)
            return False

        # update in-memory state and DB only after the filesystem rename succeeds
        self.linkpath = new_workspace_link
        self.workspacemodel.name = newworkspacename
        self.workspacemodel.link = new_link_name
        self.workspacemodel.save()
        return True

    def updateDescription(self, newdescription):
        """
        Update the human-readable description of this workspace in the DB.
        Returns True on success, False on any failure.

        unit tests: test_workspace_update_description
        """
        if self.workspacemodel is None:
            print_error("workspacemodel is none")
            return False
        if newdescription is None:
            print_error("new description is none")
            return False
        self.workspacemodel.desc = newdescription
        self.workspacemodel.save()
        return True
