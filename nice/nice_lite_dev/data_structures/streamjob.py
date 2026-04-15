# global imports
import os
import copy
import time
from time import gmtime, strftime
from django.utils import timezone

# local imports
from ..helpers import *
from ..models  import JobModel
from .simple   import SIMPLEStream


class StreamJob:
    """
    Represents a single SIMPLE stream processing job within a dataset.

    A stream job orchestrates a multi-stage cryo-EM pipeline (preprocessing,
    optics assignment, picking, 2D classification, etc.). Each stage runs as a
    separate subprocess managed by SIMPLEStream and communicates its status back
    via heartbeat updates written to the DB.

    The job directory sits inside the dataset directory and is named
    <jcnt>_simple_stream/, where jcnt is a per-dataset counter.

    Lifecycle:
      - Create  : StreamJob().new(dataset, args)
      - Load    : StreamJob(id)
      - Delete  : job.delete()          — moves directory to dataset TRASH
      - Control : terminate_master(), terminate_process(), restart_process()
      - Updates : update_stats()        — called by the running job via the API
                  update()              — called by the GUI to adjust thresholds
    """

    def __init__(self, id=None):
        self.id       = 0
        self.jobmodel = None
        self.absdir   = None  # absolute path to the job directory
        if id is not None:
            self.id = id
            self.load()

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def get_absdir(self):
        return self.absdir

    def get_jobmodel(self):
        return self.jobmodel

    def get_status(self):
        """
        Return the master status string for this job.
        If the job has been running or unknown for more than 60 seconds without
        a heartbeat, it is automatically marked as failed.
        """
        if self.jobmodel is None:
            return "unknown"
        if self.jobmodel.master_status in ("running", "unknown"):
            if self.jobmodel.master_heartbeat < int(time.time()) - 60:
                print_error("heartbeat from job " + str(self.jobmodel.id) + " not heard for 60 seconds. Failing job")
                self.jobmodel.master_status = "failed"
                self.jobmodel.save()
        return self.jobmodel.master_status

    def get_master_update(self):
        """Return the master_update dict, used by the running job to poll for GUI commands."""
        if self.jobmodel is None:
            return None
        return self.jobmodel.master_update

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def new(self, dataset, args):
        """
        Create a new stream job inside the given dataset.

        Validates the UI JSON from args, creates the job directory, initialises
        all stage statuses to 'queued', saves the DB record, then launches the
        SIMPLEStream subprocess.
        Returns True on success, False on any failure.
        """
        if dataset is None:
            print_error("dataset is None")
            return False

        datasetmodel = dataset.get_datasetmodel()
        disp         = datasetmodel.jcnt + 1
        datasetmodel.jcnt = disp
        dirc   = str(disp) + "_simple_stream"
        absdir = os.path.join(dataset.get_absdir(), dirc)

        # validate the UI JSON before touching the filesystem or DB
        simplestream = SIMPLEStream(copy.deepcopy(args))
        if not simplestream.loadUIJSON():
            return False

        if not ensure_directory(absdir):
            return False

        # initialise all stage statuses and heartbeats
        jobmodel = JobModel(dset=datasetmodel, disp=disp, dirc=dirc, cdat=timezone.now(), args=args)
        jobmodel.master_status               = "queued"
        jobmodel.preprocessing_status        = "queued"
        jobmodel.optics_assignment_status    = "queued"
        jobmodel.initial_picking_status      = "queued"
        jobmodel.generate_pickrefs_status    = "queued"
        jobmodel.reference_picking_status    = "queued"
        jobmodel.particle_sieving_status     = "queued"
        jobmodel.classification_2D_status    = "queued"
        jobmodel.master_heartbeat            = 0
        jobmodel.preprocessing_heartbeat     = 0
        jobmodel.optics_assignment_heartbeat = 0
        jobmodel.initial_picking_heartbeat   = 0
        jobmodel.generate_pickrefs_heartbeat = 0
        jobmodel.reference_picking_heartbeat = 0
        jobmodel.particle_sieving_heartbeat  = 0
        jobmodel.classification_2D_heartbeat = 0
        jobmodel.save()
        datasetmodel.save()

        self.id     = jobmodel.id
        self.absdir = absdir

        if not simplestream.start(absdir, self.id):
            return False
        return True

    def delete(self):
        """
        Soft-delete this job by moving its directory into the dataset TRASH.
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
    # GUI parameter updates
    # ------------------------------------------------------------------

    def update_description(self, description):
        """Update the human-readable description of this job in the DB."""
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        if description is None:
            print_error("description is none")
            return False
        self.jobmodel.desc = description
        self.jobmodel.save()
        return True

    def update(self, ctfres=None, astigmatism=None, icescore=None, increase_nmics=None):
        """
        Push GUI-driven parameter changes into the master_update dict so the
        running job picks them up on its next poll.

        ctfres         — CTF resolution threshold (Angstroms)
        astigmatism    — astigmatism threshold
        icescore       — ice fraction threshold
        increase_nmics — if set, increments the micrograph count for pickrefs
        """
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        master_update           = self.jobmodel.master_update
        preprocessing_stats     = self.jobmodel.preprocessing_stats
        generate_pickrefs_stats = self.jobmodel.generate_pickrefs_stats
        if ctfres is not None:
            master_update["ctfresthreshold"]      = ctfres
            preprocessing_stats["cutoff_ctf_res"] = ctfres
        if astigmatism is not None:
            master_update["astigthreshold"]           = astigmatism
            preprocessing_stats["cutoff_astigmatism"] = astigmatism
        if icescore is not None:
            master_update["icefracthreshold"]       = icescore
            preprocessing_stats["cutoff_ice_score"] = icescore
        if increase_nmics is not None:
            master_update["increase_nmics"] = master_update.get("increase_nmics", 0) + 1
            generate_pickrefs_stats["user_input"] = False
            generate_pickrefs_stats["stage"]      = "using more particles"
        self.jobmodel.master_update           = master_update
        self.jobmodel.preprocessing_stats     = preprocessing_stats
        self.jobmodel.generate_pickrefs_stats = generate_pickrefs_stats
        self.jobmodel.save()

    def select_pickrefs(self, final_selection):
        """Write the user's pickref selection into master_update for the job to consume."""
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        master_update = self.jobmodel.master_update
        master_update["pickrefs_selection"] = final_selection
        self.jobmodel.master_update = master_update
        stats = self.jobmodel.generate_pickrefs_stats
        stats["user_input"] = False
        stats["stage"]      = "saving selection"
        self.jobmodel.generate_pickrefs_stats = stats
        self.jobmodel.save()
        return True

    # UNUSED - called on Job, not StreamJob
    # def select_moldiam_initial_pick(self, diameter):
    # def update_moldiam_refine_initial_pick(self, diameter):
    # def select_refs_generate_pickrefs(self, final_selection, final_selection_source):
    # def snapshot_classification_2D(self, snapshot_selection, snapshot_iteration):
    # def selection_classification_2D(self, final_deselection, project, dataset, final_selection_ptcls):

    # ------------------------------------------------------------------
    # Stage heartbeat / stats ingestion (called by the running job via API)
    # ------------------------------------------------------------------

    def update_stats(self, stats_json):
        """
        Ingest a stats payload from the running SIMPLEStream process.

        The payload may contain any combination of:
          - "stream_heartbeat" : dict of per-stage heartbeat entries
          - "preprocessing"    : preprocessing stats dict
          - "optics_assignment": optics assignment stats dict
          - "initial_picking"  : initial picking stats dict
          - "opening2D"        : pickrefs generation stats dict
          - "reference_picking": reference picking stats dict

        When a stage transitions to "running", any pending restart command for
        that stage is cleared from master_update. When opening2D is not running,
        any pending nmics/pickref-selection commands are also cleared.

        Returns True on success, False if the job is not loaded.
        """
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False

        updated = False

        if "stream_heartbeat" in stats_json:
            updated       = True
            heartbeat     = stats_json["stream_heartbeat"]
            master_update = self.jobmodel.master_update
            self.jobmodel.master_heartbeat = int(time.time())
            self.jobmodel.master_stats     = heartbeat
            if "master" in heartbeat:
                status, _ = analyse_heartbeat(heartbeat["master"])
                self.jobmodel.status        = status
                self.jobmodel.master_status = status
            if "preprocessing" in heartbeat:
                status, _ = analyse_heartbeat(heartbeat["preprocessing"])
                self.jobmodel.preprocessing_status = status
                if status == "running":
                    master_update.pop("restart_preprocess", None)
            if "assign_optics" in heartbeat:
                status, _ = analyse_heartbeat(heartbeat["assign_optics"])
                self.jobmodel.optics_assignment_status = status
                if status == "running":
                    master_update.pop("restart_optics_assignment", None)
            if "initial_picking" in heartbeat:
                status, _ = analyse_heartbeat(heartbeat["initial_picking"])
                self.jobmodel.initial_picking_status = status
                if status == "running":
                    master_update.pop("restart_opening2D", None)
            if "opening2D" in heartbeat:
                status, _ = analyse_heartbeat(heartbeat["opening2D"])
                self.jobmodel.generate_pickrefs_status = status
                if status == "running":
                    master_update.pop("restart_opening2D", None)
                else:
                    # clear pending user inputs once the stage is no longer running
                    master_update.pop("increase_nmics",     None)
                    master_update.pop("pickrefs_selection", None)
            if "reference_picking" in heartbeat:
                status, _ = analyse_heartbeat(heartbeat["reference_picking"])
                self.jobmodel.reference_picking_status = status
                if status == "running":
                    master_update.pop("restart_reference_picking", None)
            if "particle_sieving" in heartbeat:
                status, _ = analyse_heartbeat(heartbeat["particle_sieving"])
                self.jobmodel.particle_sieving_status = status
                if status == "running":
                    master_update.pop("restart_particle_sieving", None)
            if "pool2D" in heartbeat:
                status, _ = analyse_heartbeat(heartbeat["pool2D"])
                self.jobmodel.classification_2D_status = status
                if status == "running":
                    master_update.pop("restart_pool2D", None)       
            self.jobmodel.master_update = master_update

        if "preprocessing" in stats_json:
            updated = True
            self.jobmodel.preprocessing_stats = stats_json["preprocessing"]
        if "optics_assignment" in stats_json:
            updated = True
            self.jobmodel.optics_assignment_stats = stats_json["optics_assignment"]
        if "initial_picking" in stats_json:
            updated = True
            self.jobmodel.initial_picking_stats = stats_json["initial_picking"]
        if "opening2D" in stats_json:
            updated = True
            self.jobmodel.generate_pickrefs_stats = stats_json["opening2D"]
        if "reference_picking" in stats_json:
            updated = True
            self.jobmodel.reference_picking_stats = stats_json["reference_picking"]
        if "particle_sieving" in stats_json:
            updated = True
            self.jobmodel.particle_sieving_stats = stats_json["particle_sieving"]
            if "initial_ref_selection" in stats_json["particle_sieving"] and "ref_selection" not in self.jobmodel.master_update:
                self.jobmodel.master_update["ref_selection"] = stats_json["particle_sieving"]["initial_ref_selection"]
        if "pool2D" in stats_json:
            updated = True
            if "snapshot" in stats_json["pool2D"]:
                snapshot    = stats_json["pool2D"]["snapshot"]
                snapshot_id = snapshot["id"]
                particle_set = next((x for x in self.jobmodel.particle_sets_stats["particle_sets"] if x["id"] == snapshot_id), None)
                if particle_set is not None and "time" not in particle_set:
                    particle_set["nptcls"]    = snapshot["snapshot_nptcls"]
                    particle_set["time"]      = snapshot["snapshot_time"]
                    particle_set["filename"]  = snapshot["snapshot_filename"]
            pool2D_stats = {k: v for k, v in stats_json["pool2D"].items() if k != "snapshot"}
            self.jobmodel.classification_2D_stats = pool2D_stats
        if updated:
            self.jobmodel.save()
        return True

    # ------------------------------------------------------------------
    # Process control
    # ------------------------------------------------------------------

    def terminate_master(self):
        """Signal the master process to terminate cleanly."""
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        if self.jobmodel.master_status != "running":
            print_error("master is not running")
            return False
        self.jobmodel.master_status = "terminating"
        master_update = self.jobmodel.master_update
        master_update["terminate"] = True
        self.jobmodel.master_update = master_update
        self.jobmodel.save()
        return True

    def terminate_process(self, term_preprocess, term_optics_assignment, term_generate_pickrefs):
        """
        Signal one or more sub-processes to terminate.
        Flags are boolean; only flagged processes are affected.
        """
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        master_update = self.jobmodel.master_update
        if term_preprocess:
            master_update["terminate_preprocess"] = True
        if term_optics_assignment:
            master_update["terminate_optics_assignment"] = True
        if term_generate_pickrefs:
            master_update["terminate_opening2D"] = True
        self.jobmodel.master_update = master_update
        self.jobmodel.save()
        return True

    def restart_process(self, restart_preprocess, restart_optics_assignment, restart_generate_pickrefs):
        """
        Signal one or more terminated sub-processes to restart.
        Clears the corresponding terminate flag before setting the restart flag.
        Flags are boolean; only flagged processes are affected.
        """
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        master_update = self.jobmodel.master_update
        if restart_preprocess:
            master_update.pop("terminate_preprocess", None)
            master_update["restart_preprocess"] = True
        if restart_optics_assignment:
            master_update.pop("terminate_optics_assignment", None)
            master_update["restart_optics_assignment"] = True
        if restart_generate_pickrefs:
            master_update.pop("terminate_opening2D", None)
            master_update["restart_opening2D"] = True
        self.jobmodel.master_update = master_update
        self.jobmodel.save()
        return True
    
    def select_sieve_particles(self, accepted_cls2D):
        """Store the user's 2D class selection for the particle-sieving stage.
        accepted_cls2D: list of class indices accepted by the user.
        Writes ref_selection into master_update so the stream picks it up.
        """
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        master_update = self.jobmodel.master_update
        master_update["ref_selection"] = accepted_cls2D
        self.jobmodel.master_update = master_update
        self.jobmodel.save()
        return True
    
    def update_mskdiam(self, mskdiam):
        """Store the mask diameter for 2D classification.
        Writes mskdiam2D into master_update so the stream picks it up.
        """
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        master_update = self.jobmodel.master_update
        master_update["mskdiam2D"] = mskdiam
        self.jobmodel.master_update = master_update
        self.jobmodel.save()
        return True
    
    def snapshot_classification_2D(self, snapshot_selection, snapshot_iteration):
        """Record a 2D-classification snapshot particle set and queue it for the stream.

        Appends a new snapshot entry to particle_sets_stats and writes the
        snapshot2D key into master_update so the stream picks it up on the
        next cycle.

        snapshot_selection: list of class indices included in the snapshot.
        snapshot_iteration: the 2D classification iteration to snapshot.
        """
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        master_update       = self.jobmodel.master_update
        particle_sets_stats = self.jobmodel.particle_sets_stats
        if not "particle_sets" in particle_sets_stats:
            particle_sets_stats["particle_sets"] = []
        setid    = len(particle_sets_stats["particle_sets"]) + 1
        projfile = "snapshot_" + str(setid) + ".simple"
        newset = {
            "id"       : setid,
            "name"     : "particle set " + str(setid),
            "type"     : "snapshot",
            "filename" : projfile
        }
        particle_sets_stats["particle_sets"].insert(0, newset)
        master_update["snapshot2D"] = {
            "id"        : setid,
            "iteration" : snapshot_iteration,
            "selection" : snapshot_selection,
            "filename"  : projfile
        }
        self.jobmodel.master_update       = master_update
        self.jobmodel.particle_sets_stats = particle_sets_stats
        self.jobmodel.save()
        return True
    
    def selection_classification_2D(self, final_deselection, final_selection_ptcls):
        """Record the final 2D-classification particle selection and persist it to disk.

        Appends a new final-selection entry to particle_sets_stats and writes
        the deselected particle indices to a per-set text file under self.absdir.

        final_deselection:     list of particle indices rejected by the user.
        final_selection_ptcls: count of particles retained in the final selection.
        """
        if self.jobmodel is None:
            print_error("jobmodel is none")
            return False
        particle_sets_stats = self.jobmodel.particle_sets_stats
        if not "particle_sets" in particle_sets_stats:
            particle_sets_stats["particle_sets"] = []
        setid    = len(particle_sets_stats["particle_sets"]) + 1
        deselfile = "particle_set_" + str(setid) + "_deselected.txt"
        newset = {
            "id"       : setid,
            "name"     : "particle set " + str(setid),
            "type"     : "final",
            "filename" : deselfile,
            "nptcls"   : final_selection_ptcls,
            "ctime"    : strftime("%Y/%m/%d %H:%M", gmtime())
        }
        particle_sets_stats["particle_sets"].insert(0, newset)
        with open(os.path.join(self.absdir, deselfile), "w") as f:
            for deselected in final_deselection:
                f.write(str(deselected) + '\n')
        self.jobmodel.particle_sets_stats = particle_sets_stats 
        self.jobmodel.save()
        return True

