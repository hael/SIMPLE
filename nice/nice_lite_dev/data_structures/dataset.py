# global imports
import os
import re
from django.utils import timezone

# local imports
from ..helpers import *
from ..models  import DatasetModel, WorkspaceModel


class Dataset:
    """
    Represents a cryo-EM dataset within a project.

    A dataset maps to:
      - a hidden directory  (.dataset_<id>/)  inside the project directory
      - a human-readable symlink (ds_<name>/) pointing to that directory
      - a TRASH/ subdirectory inside the dataset directory used for soft-deletes

    Lifecycle:
      - Create : Dataset().new(project, user)
      - Load   : Dataset(id)
      - Delete : dataset.delete()     — moves files to project TRASH, removes DB record
      - Rename : dataset.rename(name) — renames the symlink and updates the DB
    """

    def __init__(self, id=None):
        # unit tests: test_dataset_init, test_dataset_init_by_id
        self.id           = 0
        self.datasetmodel = None
        self.absdir       = None  # absolute path to the hidden dataset directory
        self.linkpath     = None  # absolute path to the human-readable symlink
        self.trashdir     = None  # absolute path to the TRASH subdirectory
        if id is not None:
            self.id = id
            self.load()

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def load(self):
        """Populate fields from the database. Resets to empty state if not found."""
        self.datasetmodel = DatasetModel.objects.filter(id=self.id).first()
        if self.datasetmodel is None:
            self.id       = 0
            self.absdir   = None
            self.linkpath = None
            self.trashdir = None
        else:
            self.absdir   = os.path.join(self.datasetmodel.proj.dirc, self.datasetmodel.dirc)
            self.linkpath = os.path.join(self.datasetmodel.proj.dirc, self.datasetmodel.link)
            self.trashdir = os.path.join(self.absdir, "TRASH")

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def in_project(self, projectid):
        """Return True if this dataset belongs to the given project id."""
        return self.datasetmodel is not None and self.datasetmodel.proj.id == projectid

    def get_id(self):
        return self.id

    def get_datasetmodel(self):
        return self.datasetmodel

    def get_absdir(self):
        return self.absdir

    def get_linkpath(self):
        return self.linkpath

    def get_trashdir(self):
        return self.trashdir

    # ------------------------------------------------------------------
    # Mutators
    # ------------------------------------------------------------------

    def new(self, project, user):
        """
        Create a new dataset inside the given project.

        Creates the hidden directory and symlink on disk, then writes the DB record.
        The DB record is rolled back if filesystem operations fail.
        Returns True on success, False on any failure.

        unit tests: test_dataset_new
        """
        projectmodel = project.get_projectmodel()
        projectdir   = project.get_absdir()
        if projectmodel is None:
            print_error("projectmodel is none")
            return False
        if not directory_exists(projectdir):
            print_error("project directory " + projectdir + " doesn't exist")
            return False

        # derive a display index and default name from existing dataset count
        display_id       = DatasetModel.objects.filter(proj=projectmodel).count() + 1
        new_dataset_name = "new dataset " + str(display_id)

        # save the record early to obtain the auto-assigned id, which is used for the directory name
        datasetmodel = DatasetModel(
            proj=projectmodel, name=new_dataset_name, disp=display_id,
            cdat=timezone.now(), mdat=timezone.now(), user=user
        )
        datasetmodel.save()
        self.id = datasetmodel.id
        if self.id == 0:
            return False

        new_dataset_dirc = ".dataset_" + str(self.id)
        new_dataset_link = "ds_" + new_dataset_name.replace(" ", "_")
        new_dataset_path = os.path.join(projectdir, new_dataset_dirc)

        if not ensure_directory(new_dataset_path):
            datasetmodel.delete()
            return False
        if not create_symlink(new_dataset_path, os.path.join(projectdir, new_dataset_link)):
            datasetmodel.delete()
            return False

        datasetmodel.dirc = new_dataset_dirc
        datasetmodel.link = new_dataset_link
        datasetmodel.save()
        self.load()
        return True

    def delete(self):
        """
        Soft-delete this dataset by moving its directory and symlink into the project TRASH.

        Also removes the DB record. If the project has no remaining datasets or
        workspaces after deletion, the project itself is also deleted.
        Returns True on success, False on any failure.

        unit tests: test_dataset_delete
        """
        if self.datasetmodel is None:
            print_error("datasetmodel is none")
            return False

        projecttrash = os.path.join(self.datasetmodel.proj.dirc, "TRASH")
        if not ensure_directory(projecttrash):
            print_error("project trash folder doesn't exist")
            return False

        trash_path = os.path.join(projecttrash, self.datasetmodel.dirc)
        trash_link = os.path.join(projecttrash, self.datasetmodel.link)

        # verify source filesystem state
        if not os.path.exists(self.absdir):
            print_error("dataset directory doesn't exist")
            return False
        if not os.path.isdir(self.absdir):
            print_error("dataset directory isn't a directory")
            return False
        if not os.path.exists(self.linkpath):
            print_error("dataset link doesn't exist")
            return False
        if not os.path.islink(self.linkpath):
            print_error("dataset link isn't a link")
            return False

        # verify trash destinations are free
        if os.path.exists(trash_path):
            print_error("dataset trash directory already exists")
            return False
        if os.path.exists(trash_link):
            print_error("trash link already exists")
            return False

        try:
            os.rename(self.absdir, trash_path)
        except OSError:
            print_error("Directory '%s' cannot be renamed")
            return False

        try:
            os.remove(self.linkpath)
            os.symlink(trash_path, trash_link)
        except OSError:
            print_error("Link '%s' cannot be renamed")
            return False

        # remove DB record; delete the parent project if it is now empty
        projectmodel = self.datasetmodel.proj
        self.datasetmodel.delete()
        remaining_workspaces = WorkspaceModel.objects.filter(proj=projectmodel).count()
        remaining_datasets   = DatasetModel.objects.filter(proj=projectmodel).count()
        if remaining_workspaces + remaining_datasets == 0:
            projectmodel.delete()

        return True

    def rename(self, newdatasetname):
        """
        Rename this dataset: renames the symlink on disk and updates the DB record.
        The directory itself is not moved — only the symlink is renamed.
        Returns True on success, False on any failure.

        unit tests: test_dataset_rename
        """
        if self.datasetmodel is None:
            print_error("datasetmodel is none")
            return False
        if newdatasetname is None:
            print_error("new dataset name is none")
            return False

        # build a filesystem-safe link name: spaces become _, non-alphanumeric stripped
        new_link_name    = "ds_" + re.sub(r'\W+', '', newdatasetname.replace(" ", "_"))
        new_dataset_link = os.path.join(self.datasetmodel.proj.dirc, new_link_name)

        try:
            os.rename(self.linkpath, new_dataset_link)
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
        self.linkpath          = new_dataset_link
        self.datasetmodel.name = newdatasetname
        self.datasetmodel.link = new_link_name
        self.datasetmodel.save()
        return True

    def updateDescription(self, newdescription):
        """
        Update the human-readable description of this dataset in the DB.
        Returns True on success, False on any failure.

        unit tests: test_dataset_update_description
        """
        if self.datasetmodel is None:
            print_error("datasetmodel is none")
            return False
        if newdescription is None:
            print_error("new description is none")
            return False
        self.datasetmodel.desc = newdescription
        self.datasetmodel.save()
        return True
