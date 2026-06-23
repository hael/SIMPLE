"""Shared base data structure for NICE stream/classic jobs."""

import os

from ..helpers import directory_exists, ensure_directory, print_error
from ..models import JobModel


class Job:
    """
    Shared base class for NICE job data structures.

    Subclasses are expected to implement `load()` and populate `self.jobmodel`
    and `self.absdir`.
    """

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def __init__(self, id=None):
        """Initialize this job. If id is provided, hydrate from the database."""
        self.id       = 0
        self.jobmodel = None
        self.absdir   = None
        self.desc     = ""
        if id is not None:
            self.id = id
            self.load()

    def __del__(self):
        """Release object-held references when this instance is garbage-collected."""
        self.jobmodel = None
        self.absdir   = None

    def load(self):
        """Populate fields from the database. Resets to empty state if not found."""
        self.jobmodel = JobModel.objects.filter(id=self.id).first()
        if self.jobmodel is None:
            self.id     = 0
            self.absdir = None
        else:
            self.absdir = os.path.join(
                self.jobmodel.dset.proj.dirc,
                self.jobmodel.dset.dirc,
                self.jobmodel.dirc
            )

    def delete(self):
        """
        Soft-delete this job by moving its directory into the workspace TRASH.
        The DB record is removed only after the filesystem move succeeds.
        Returns True on success, False on any failure.
        """
        if self.jobmodel is None:
            print_error("jobmodel is None")
            return False

        # derive trashdir from the already-loaded related model to avoid an extra DB query
        trashdir   = os.path.join(self.jobmodel.dset.proj.dirc, self.jobmodel.dset.dirc, "TRASH")
        trash_path = os.path.join(trashdir, self.jobmodel.dirc)

        if not ensure_directory(trashdir):
            print_error("failed to create trash directory")
            return False
        if directory_exists(trash_path):
            print_error("trash directory already exists")
            return False
        if not directory_exists(self.absdir):
            print_error("job directory doesn't exist")
            return False
        try:
            os.rename(self.absdir, trash_path)
        except OSError:
            print_error("failed to rename job directory: " + self.absdir)
            return False

        self.jobmodel.delete()
        return True

    # ------------------------------------------------------------------
    # Getters
    # ------------------------------------------------------------------

    def get_status(self):
        """Return the current status string for this job."""
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return None
        return self.jobmodel.status

    def get_absdir(self):
        """Return the absolute directory path for this job."""
        if self.absdir is None:
            print_error("absdir is none")
            return None
        return self.absdir

    def get_jobmodel(self):
        """Return the JobModel instance for this job."""
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return None
        return self.jobmodel

    def get_master_update(self):
        """Return the master_update dict, used by the running job to poll for GUI commands."""
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return None
        return self.jobmodel.master_update

    # ------------------------------------------------------------------
    # Setters
    # ------------------------------------------------------------------

    def set_description(self, description):
        """Update the human-readable description of this job in the DB."""
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        if description is None:
            print_error("description is none")
            return False
        self.desc          = description
        self.jobmodel.desc = self.desc
        self.jobmodel.save()
        return True

    def set_id_from_request(self, request, key="jobid"):
        test_id_str = None
        if key in request.POST:
            test_id_str = request.POST[key]
        elif key in request.GET:
            test_id_str = request.GET[key]
        if test_id_str is not None and test_id_str.isnumeric():
            self.id = int(test_id_str)