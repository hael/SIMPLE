# global imports
import os
import copy
import datetime
from time import gmtime, strftime
from django.utils import timezone

# local imports
from ..models import JobModel, DatasetModel
from .simple  import SIMPLEStream

class Job:

    id       = 0
    disp     = 0
    name     = ""
    desc     = ""
    dirc     = ""
    cdat     = ""
    args     = {}
    dset     = 0
    status   = ""
    preprocessing_status     = "" 
    preprocessing_stats      = {}
    optics_assignment_status = "" 
    optics_assignment_stats  = {}
    initial_picking_status   = "" 
    initial_picking_stats    = {}
    generate_pickrefs_status = "" 
    generate_pickrefs_stats  = {}
    reference_picking_status = "" 
    reference_picking_stats  = {}
    particle_sieving_status  = "" 
    particle_sieving_stats   = {}
    classification_2D_status = "" 
    classification_2D_stats  = {}
    particle_sets_stats      = {}

    def __init__(self, id=None, request=None):
        if id is not None:
            self.id = id
            self.load()
        elif request is not None:
            self.setIDFromRequest(request)
            self.load()

    def setIDFromRequest(self, request):
        if "jobid" in request.POST:
            test_id_str = request.POST["jobid"]
        elif "jobid" in request.GET:
            test_id_str = request.GET["jobid"]
        if test_id_str.isnumeric():
            self.id = int(test_id_str)
    
    def getAbsDir(self):
        return os.path.join(self.dset.proj.dirc, self.dset.dirc, self.dirc)

    def new(self, request, project, dataset):
        self.args = {} # ensure empty
        for key, value in request.POST.items():
            if "csrfmiddlewaretoken" not in key and value != "":
                self.args[key] = value
        datasetmodel = DatasetModel.objects.filter(id=dataset.id).first()
        #jobmodels = JobModel.objects.filter(dset=datasetmodel)
      #  self.disp = jobmodels.count() + 1
        # handle no downsampling
        if "smpd" in self.args and "smpd_downscale" in self.args:
            if float(self.args["smpd"]) > float(self.args["smpd_downscale"]):
                self.args["smpd_downscale"] = self.args["smpd"]
        self.disp = datasetmodel.jcnt + 1
        datasetmodel.jcnt = self.disp
        datasetmodel.save()
        jobmodel = JobModel(dset=datasetmodel, disp=self.disp, cdat=timezone.now(), args=self.args)
        jobmodel.save()
        self.id = jobmodel.id
        self.dirc = str(self.disp) + "_simple_stream"
        if not self.createDir(os.path.join(project.dirc, dataset.dirc)):
            return False
        jobmodel.dirc = self.dirc
        jobmodel.save()
        simplestream = SIMPLEStream()
        if not simplestream.start(copy.deepcopy(self.args), os.path.join(project.dirc, dataset.dirc, self.dirc), self.id):
            return False
        self.preprocessing_status     = "running"
        self.optics_assignment_status = "running"
        self.initial_picking_status   = "running"
        self.generate_pickrefs_status = "running"
        self.reference_picking_status = "running"
        self.particle_sieving_status  = "running"
        self.classification_2D_status = "running"
        jobmodel.preprocessing_status     = "running"
        jobmodel.optics_assignment_status = "running"
        if simplestream.skip_refgen:
            jobmodel.initial_picking_status   = "skipped"
            jobmodel.generate_pickrefs_status = "skipped"
        else:
            jobmodel.initial_picking_status   = "running"
            jobmodel.generate_pickrefs_status = "running"
        jobmodel.reference_picking_status = "running"
        jobmodel.particle_sieving_status  = "running"
        jobmodel.classification_2D_status = "running"
        jobmodel.preprocessing_heartbeat     = timezone.now()
        jobmodel.optics_assignment_heartbeat = timezone.now()
        jobmodel.initial_picking_heartbeat   = timezone.now()
        jobmodel.generate_pickrefs_heartbeat = timezone.now()
        jobmodel.reference_picking_heartbeat = timezone.now()
        jobmodel.particle_sieving_heartbeat  = timezone.now()
        jobmodel.classification_2D_heartbeat = timezone.now()
        jobmodel.save()
        return True
    
    def load(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            self.name = jobmodel.name
            self.desc = jobmodel.desc
            self.dirc = jobmodel.dirc
            self.cdat = jobmodel.cdat
            self.args = jobmodel.args
            self.dset = jobmodel.dset
            self.disp = jobmodel.disp
            if jobmodel.preprocessing_status != "finished" and jobmodel.preprocessing_status != "failed":
                if jobmodel.preprocessing_heartbeat + datetime.timedelta(minutes=1) < timezone.now():
                    jobmodel.preprocessing_status = "failed"
                    jobmodel.save()
            self.preprocessing_status     = jobmodel.preprocessing_status
            self.preprocessing_stats      = jobmodel.preprocessing_stats
            if jobmodel.optics_assignment_status != "finished" and jobmodel.optics_assignment_status != "failed":
                if jobmodel.optics_assignment_heartbeat + datetime.timedelta(minutes=1) < timezone.now():
                    jobmodel.optics_assignment_status = "failed"
                    jobmodel.save()
            self.optics_assignment_status = jobmodel.optics_assignment_status 
            self.optics_assignment_stats  = jobmodel.optics_assignment_stats
            if jobmodel.initial_picking_status != "finished" and jobmodel.initial_picking_status != "failed" and jobmodel.initial_picking_status != "skipped":
                if jobmodel.initial_picking_heartbeat + datetime.timedelta(minutes=1) < timezone.now():
                    jobmodel.initial_picking_status = "failed"
                    jobmodel.save()
            self.initial_picking_status   = jobmodel.initial_picking_status
            self.initial_picking_stats    = jobmodel.initial_picking_stats
            if jobmodel.generate_pickrefs_status != "finished" and jobmodel.generate_pickrefs_status != "failed" and jobmodel.initial_picking_status != "skipped":
                if jobmodel.generate_pickrefs_heartbeat + datetime.timedelta(minutes=1) < timezone.now():
                    jobmodel.generate_pickrefs_status = "failed"
                    jobmodel.save()
            self.generate_pickrefs_status = jobmodel.generate_pickrefs_status
            self.generate_pickrefs_stats  = jobmodel.generate_pickrefs_stats
            if jobmodel.reference_picking_status != "finished" and jobmodel.reference_picking_status != "failed":
                if jobmodel.reference_picking_heartbeat + datetime.timedelta(minutes=1) < timezone.now():
                    jobmodel.reference_picking_status = "failed"
                    jobmodel.save()
            self.reference_picking_status = jobmodel.reference_picking_status
            self.reference_picking_stats  = jobmodel.reference_picking_stats
            if jobmodel.particle_sieving_status != "finished" and jobmodel.particle_sieving_status != "failed":
                if jobmodel.particle_sieving_heartbeat + datetime.timedelta(minutes=1) < timezone.now():
                    jobmodel.particle_sieving_status = "failed"
                    jobmodel.save()
            self.particle_sieving_status  = jobmodel.particle_sieving_status
            self.particle_sieving_stats   = jobmodel.particle_sieving_stats
            if jobmodel.classification_2D_status != "finished" and jobmodel.classification_2D_status != "failed":
                if jobmodel.classification_2D_heartbeat + datetime.timedelta(minutes=1) < timezone.now():
                    jobmodel.classification_2D_status = "failed"
                    jobmodel.save()
            self.classification_2D_status = jobmodel.classification_2D_status
            self.classification_2D_stats  = jobmodel.classification_2D_stats
            self.particle_sets_stats      = jobmodel.particle_sets_stats
            # global status
            running  = False
            finished = False
            failed   = False
            if jobmodel.preprocessing_status == "finished":
                finished = True
            elif jobmodel.preprocessing_status == "failed":
                failed = True
            else:
                running = True
            if jobmodel.optics_assignment_status == "finished":
                finished = True
            elif jobmodel.optics_assignment_status == "failed":
                failed = True
            else:
                running = True   
            if jobmodel.initial_picking_status == "finished":
                finished = True
            elif jobmodel.initial_picking_status == "skipped":
                finished = True
            elif jobmodel.initial_picking_status == "failed":
                failed = True
            else:
                running = True
            if jobmodel.generate_pickrefs_status == "finished":
                finished = True
            elif jobmodel.generate_pickrefs_status == "skipped":
                finished = True  
            elif jobmodel.generate_pickrefs_status == "failed":
                failed = True
            else:
                running = True
            if jobmodel.reference_picking_status == "finished":
                finished = True
            elif jobmodel.reference_picking_status == "failed":
                failed = True
            else:
                running = True
            if jobmodel.particle_sieving_status == "finished":
                finished = True
            elif jobmodel.particle_sieving_status == "failed":
                failed = True
            else:
                running = True
            if jobmodel.classification_2D_status == "finished":
                finished = True
            elif jobmodel.classification_2D_status == "failed":
                failed = True
            else:
                running = True
            if not failed and running:
                self.status = "running"
            elif finished and not running and not failed:
                self.status = "finished"
            elif not running:
                self.status = "failed"
            else:
                self.status = "alarm" 
            if jobmodel.status != self.status:
                jobmodel.status = self.status
                jobmodel.save()

    def createDir(self, parent_dir):
        if self.dirc == "":
            return False
        
        if not os.path.exists(parent_dir):
            return False
        
        if not os.path.isdir(parent_dir):
            return False
    
        new_dir_path = os.path.join(parent_dir, self.dirc)

        try:
            os.makedirs(new_dir_path, exist_ok=False)
        except OSError as error:
            print("Directory '%s' can not be created", new_dir_path)
            return False
        return True
    
    def delete(self, project, dataset):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            jobmodel.delete()
        if dataset.ensureTrashfolder(project):
            job_path   = os.path.join(project.dirc, dataset.dirc, self.dirc)
            trash_path = os.path.join(dataset.trashfolder, self.dirc)
            if not os.path.exists(job_path):
                return 
            if not os.path.isdir(job_path):
                return 
            if os.path.exists(trash_path):
                return 
            if os.path.isdir(trash_path):
                return
            try:
                os.rename(job_path, trash_path)
            except OSError as error:
                print("Directory '%s' can not be renamed")
                return
        return
    
    def update_description(self, description):
        self.desc = description
        jobmodel = JobModel.objects.filter(id=self.id).first()
        jobmodel.desc = self.desc
        jobmodel.save()

    def update_ctfres(self, ctfres):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            updated = jobmodel.preprocessing_update
            updated["ctfresthreshold"] = float(ctfres)
            jobmodel.preprocessing_update = updated
            updated = jobmodel.initial_picking_update
            updated["ctfresthreshold"] = float(ctfres)
            jobmodel.initial_picking_update = updated
            updated = jobmodel.reference_picking_update
            updated["ctfresthreshold"] = float(ctfres)
            jobmodel.reference_picking_update = updated
            self.preprocessing_stats["cutoff_ctf_res"] = float(ctfres)
            jobmodel.preprocessing_stats = self.preprocessing_stats
            jobmodel.save()

    def update_astigmatism(self, astigmatism):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            updated = jobmodel.preprocessing_update
            updated["astigthreshold"] = float(astigmatism)
            jobmodel.preprocessing_update = updated
            updated = jobmodel.initial_picking_update
            updated["astigthreshold"] = float(ctfres)
            jobmodel.initial_picking_update = updated
            updated = jobmodel.reference_picking_update
            updated["astigthreshold"] = float(ctfres)
            jobmodel.reference_picking_update = updated
            self.preprocessing_stats["cutoff_astigmatism"] = float(astigmatism)
            jobmodel.preprocessing_stats = self.preprocessing_stats
            jobmodel.save()
    
    def update_icescore(self, icescore):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            updated = jobmodel.preprocessing_update
            updated["icefracthreshold"] = float(icescore)
            jobmodel.preprocessing_update = updated
            updated = jobmodel.initial_picking_update
            updated["icefracthreshold"] = float(ctfres)
            jobmodel.initial_picking_update = updated
            updated = jobmodel.reference_picking_update
            updated["icefracthreshold"] = float(ctfres)
            jobmodel.reference_picking_update = updated
            self.preprocessing_stats["cutoff_ice_score"] = float(icescore)
            jobmodel.preprocessing_stats = self.preprocessing_stats
            jobmodel.save()

    def update_mskdiam(self, mskdiam):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            updated = jobmodel.classification_2D_update
            updated["mskdiam"] = float(mskdiam)
            jobmodel.classification_2D_update = updated
            for cls2D in self.classification_2D_stats["latest_cls2D"]:
                cls2D["mskdiam"] = mskdiam
            jobmodel.classification_2D_stats["stage"] = "updating mask diameter"
            jobmodel.classification_2D_stats = self.classification_2D_stats
            jobmodel.save() 
        
    def updateStats(self, stats_json):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        response = {}
        if "preprocessing" in stats_json:
            jobmodel.preprocessing_stats     = stats_json["preprocessing"]
            jobmodel.preprocessing_heartbeat = timezone.now()
            if "terminate" in stats_json["preprocessing"]:
                jobmodel.preprocessing_status = "finished"
            else:
                jobmodel.preprocessing_status = "running"
            response = jobmodel.preprocessing_update
            jobmodel.preprocessing_update = {}
        elif "preprocessing_heartbeat" in stats_json:
            jobmodel.preprocessing_heartbeat = timezone.now()
            jobmodel.preprocessing_status = "running"
            response = jobmodel.preprocessing_update
            jobmodel.preprocessing_update = {}
        if "optics_assignment" in stats_json:
            jobmodel.optics_assignment_stats     = stats_json["optics_assignment"]
            jobmodel.optics_assignment_heartbeat = timezone.now()
            if "terminate" in stats_json["optics_assignment"]:
                jobmodel.optics_assignment_status = "finished"
            else:
                jobmodel.optics_assignment_status = "running"
            response = jobmodel.optics_assignment_update
            jobmodel.optics_assignment_update = {}
        elif "optics_assignment_heartbeat" in stats_json:
            jobmodel.optics_assignment_heartbeat = timezone.now()
            jobmodel.optics_assignment_status = "running"
            response = jobmodel.optics_assignment_update
            jobmodel.optics_assignment_update = {}
        if "initial_picking" in stats_json:
            jobmodel.initial_picking_stats     = stats_json["initial_picking"]
            jobmodel.initial_picking_heartbeat = timezone.now()
            if "terminate" in stats_json["initial_picking"]:
                jobmodel.initial_picking_status = "finished"
            else:
                jobmodel.initial_picking_status = "running"
            response = jobmodel.initial_picking_update
            jobmodel.initial_picking_update = {}
        elif "initial_picking_heartbeat" in stats_json:
            jobmodel.initial_picking_heartbeat = timezone.now()
            jobmodel.initial_picking_status = "running"
            response = jobmodel.initial_picking_update
            jobmodel.initial_picking_update = {}
        if "generate_picking_refs" in stats_json:
            jobmodel.generate_pickrefs_stats     = stats_json["generate_picking_refs"]
            jobmodel.generate_pickrefs_heartbeat = timezone.now()
            if "terminate" in stats_json["generate_picking_refs"]:
                jobmodel.generate_pickrefs_status = "finished"
            else:
                jobmodel.generate_pickrefs_status = "running"
            response = jobmodel.generate_pickrefs_update
            jobmodel.generate_pickrefs_update = {}
        elif "generate_picking_refs_heartbeat" in stats_json:
            jobmodel.generate_pickrefs_heartbeat = timezone.now()
            jobmodel.generate_pickrefs_status = "running"
            response = jobmodel.generate_pickrefs_update
            jobmodel.generate_pickrefs_update = {}
        if "pick_extract" in stats_json:
            jobmodel.reference_picking_stats     = stats_json["pick_extract"]
            jobmodel.reference_picking_heartbeat = timezone.now()
            if "terminate" in stats_json["pick_extract"]:
                jobmodel.reference_picking_status = "finished"
            else:
                jobmodel.reference_picking_status = "running"
            response = jobmodel.reference_picking_update
            jobmodel.reference_picking_update = {}
        elif "pick_extract_heartbeat" in stats_json:
            jobmodel.reference_picking_heartbeat = timezone.now()
            jobmodel.reference_picking_status = "running"
            response = jobmodel.reference_picking_update
            jobmodel.reference_picking_update = {}
        if "sieve_cavgs" in stats_json:
            jobmodel.particle_sieving_stats     = stats_json["sieve_cavgs"]
            jobmodel.particle_sieving_heartbeat = timezone.now()
            if "terminate" in stats_json["sieve_cavgs"]:
                jobmodel.particle_sieving_status = "finished"
            else:
                jobmodel.particle_sieving_status = "running"
            response = jobmodel.particle_sieving_update
            jobmodel.particle_sieving_update = {}
        elif "sieve_cavgs_heartbeat" in stats_json:
            jobmodel.particle_sieving_heartbeat = timezone.now()
            jobmodel.particle_sieving_status = "running"
            response = jobmodel.particle_sieving_update
            jobmodel.particle_sieving_update = {}
        if "classification_2D" in stats_json:
            jobmodel.classification_2D_stats     = stats_json["classification_2D"]
            jobmodel.classification_2D_heartbeat = timezone.now()
            if "terminate" in stats_json["classification_2D"]:
                jobmodel.classification_2D_status = "finished"
            else:
                jobmodel.classification_2D_status = "running"
            response = jobmodel.classification_2D_update
            jobmodel.classification_2D_update = {}
            # snapshot stats
            if "snapshot_filename" in stats_json["classification_2D"] and "snapshot_nptcls" in stats_json["classification_2D"] and "snapshot_time" in stats_json["classification_2D"]:
                for particle_set in self.particle_sets_stats["particle_sets"]:
                    if "filename" in particle_set and particle_set["filename"] == stats_json["classification_2D"]["snapshot_filename"]:
                        particle_set["nptcls"] = stats_json["classification_2D"]["snapshot_nptcls"]
                        particle_set["ctime"]  = stats_json["classification_2D"]["snapshot_time"]
                        break
                jobmodel.particle_sets_stats  = self.particle_sets_stats
        elif "classification_2D_heartbeat" in stats_json:
            jobmodel.classification_2D_heartbeat = timezone.now()
            jobmodel.classification_2D_status = "running"
            response = jobmodel.classification_2D_update
            jobmodel.classification_2D_update = {}
        jobmodel.save()
        return response

    def select_moldiam_initial_pick(self, diameter):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            updated = jobmodel.initial_picking_update
            updated["moldiam"] = int(diameter)
            updated["interactive"] = "no"
            jobmodel.initial_picking_update = updated
            self.initial_picking_stats["user_input"] = False
            self.initial_picking_stats["stage"]      = "updating picking parameters"
            self.initial_picking_stats["latest_picked_micrographs"] = []
            jobmodel.initial_picking_stats = self.initial_picking_stats
            jobmodel.save()

    def update_moldiam_refine_initial_pick(self, diameter):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            updated = jobmodel.initial_picking_update
            updated["moldiam_refine"] = int(diameter)
            jobmodel.initial_picking_update = updated
            self.initial_picking_stats["user_input"]                = False
            self.initial_picking_stats["stage"]                     = "updating picking parameters"
            self.initial_picking_stats["latest_picked_micrographs"] = []
            jobmodel.initial_picking_stats = self.initial_picking_stats
            jobmodel.save()

    def select_refs_generate_pickrefs(self, final_selection, final_selection_source):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            updated = jobmodel.generate_pickrefs_update
            updated["final_selection"]         = final_selection
            updated["final_selection_source"]  = final_selection_source
            updated["terminate"]               = True
            jobmodel.generate_pickrefs_update  = updated
            self.generate_pickrefs_stats["user_input"] = False
            self.generate_pickrefs_stats["stage"]      = "saving selection"
            jobmodel.generate_pickrefs_stats = self.generate_pickrefs_stats
            jobmodel.save()
            self.terminate_initial_pick()

    def regenerate_pickrefs(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            updated = jobmodel.generate_pickrefs_update
            updated["increase_nmics"]          = True
            updated["terminate"]               = False
            jobmodel.generate_pickrefs_update  = updated
            self.generate_pickrefs_stats["user_input"] = False
            self.generate_pickrefs_stats["stage"]      = "using more particles"
            jobmodel.generate_pickrefs_stats = self.generate_pickrefs_stats
            jobmodel.save()

    def snapshot_classification_2D(self, snapshot_selection, snapshot_iteration):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            if not "particle_sets" in self.particle_sets_stats:
                self.particle_sets_stats["particle_sets"] = []
            setid = len(self.particle_sets_stats["particle_sets"]) + 1
            updated = jobmodel.classification_2D_update
            updated["snapshot_iteration"]  = snapshot_iteration
            updated["snapshot_selection"]  = snapshot_selection
            updated["snapshot_filename"]   = "snapshot_" + str(setid) + ".simple"
            jobmodel.classification_2D_update  = updated
            self.classification_2D_stats["stage"]  = "saving selection"
            jobmodel.classification_2D_stats = self.classification_2D_stats
            newset = {
                "id"       : setid,
                "name"     : "particle set " + str(setid),
                "type"     : "snapshot",
                "filename" : updated["snapshot_filename"]  
            }
            self.particle_sets_stats["particle_sets"].insert(0, newset)
            jobmodel.particle_sets_stats  = self.particle_sets_stats
            jobmodel.save()

    def selection_classification_2D(self, final_deselection, project, dataset, final_selection_ptcls):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            if not "particle_sets" in self.particle_sets_stats:
                self.particle_sets_stats["particle_sets"] = []
            setid = len(self.particle_sets_stats["particle_sets"]) + 1
            deselfile = "particle_set_" + str(setid) + "_deselected.txt"
            with open(os.path.join(project.dirc, dataset.dirc, self.dirc, deselfile), "w") as f:
                for deselected in final_deselection:
                    f.write(str(deselected) + '\n')
            newset = {
                "id"       : setid,
                "name"     : "particle set " + str(setid),
                "type"     : "final",
                "filename" : deselfile,
                "nptcls"   : final_selection_ptcls, 
                "ctime"    : strftime("%Y/%m/%d %H:%M", gmtime())
            }
            self.particle_sets_stats["particle_sets"].insert(0, newset)
            jobmodel.particle_sets_stats = self.particle_sets_stats
            jobmodel.save()

    def select_sieve_particles(self, accepted_cls2D, rejected_cls2D):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            updated = jobmodel.particle_sieving_update
            updated["accepted_cls2D"] = accepted_cls2D
            updated["rejected_cls2D"] = rejected_cls2D
            jobmodel.particle_sieving_update  = updated
            self.particle_sieving_stats["user_input"] = False
            self.particle_sieving_stats["stage"]      = "saving selection"
            jobmodel.particle_sieving_stats = self.particle_sieving_stats
            jobmodel.save()

    def terminate(self):
        self.terminate_preprocess()
        self.terminate_optics()
        self.terminate_initial_pick()
        self.terminate_generate_pickrefs()
        self.terminate_reference_picking()
        self.terminate_sieve_particles()
        self.terminate_classification_2D()

    def terminate_preprocess(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None and jobmodel.preprocessing_status == "running":
            jobmodel.preprocessing_status     = "terminating"
            jobmodel.preprocessing_update     = {"terminate":True}
            self.preprocessing_stats["stage"] = "terminating"
            jobmodel.preprocessing_stats = self.preprocessing_stats
            jobmodel.save()

    def terminate_optics(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None and jobmodel.optics_assignment_status == "running":
            jobmodel.optics_assignment_status     = "terminating"
            jobmodel.optics_assignment_update     = {"terminate":True}
            self.optics_assignment_stats["stage"] = "terminating"
            jobmodel.optics_assignment_stats = self.optics_assignment_stats
            jobmodel.save()

    def terminate_initial_pick(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None and jobmodel.initial_picking_status == "running":
            jobmodel.initial_picking_status     = "terminating"
            jobmodel.initial_picking_update     = {"terminate":True}
            self.initial_picking_stats["stage"] = "terminating"
            jobmodel.initial_picking_stats = self.initial_picking_stats
            jobmodel.generate_pickrefs_status     = "terminating"
            jobmodel.generate_pickrefs_status     = {"terminate":True}
            self.generate_pickrefs_stats["stage"] = "terminating"
            jobmodel.generate_pickrefs_stats = self.generate_pickrefs_stats
            jobmodel.save()
    
    def terminate_generate_pickrefs(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None and jobmodel.generate_pickrefs_status == "running":
            jobmodel.generate_pickrefs_status     = "terminating"
            jobmodel.generate_pickrefs_update     = {"terminate":True}
            self.generate_pickrefs_stats["stage"] = "terminating"
            jobmodel.generate_pickrefs_stats = self.generate_pickrefs_stats
            jobmodel.initial_picking_status     = "terminating"
            jobmodel.initial_picking_update     = {"terminate":True}
            self.initial_picking_stats["stage"] = "terminating"
            jobmodel.initial_picking_stats = self.initial_picking_stats
            jobmodel.save()      
    
    def terminate_reference_picking(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None and jobmodel.reference_picking_status == "running":
            jobmodel.reference_picking_status     = "terminating"
            jobmodel.reference_picking_update     = {"terminate":True}
            self.reference_picking_stats["stage"] = "terminating"
            jobmodel.reference_picking_stats = self.reference_picking_stats
            jobmodel.save()      
    
    def terminate_sieve_particles(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None and jobmodel.particle_sieving_status == "running":
            jobmodel.particle_sieving_status     = "terminating"
            jobmodel.particle_sieving_update     = {"terminate":True}
            self.particle_sieving_stats["stage"] = "terminating"
            jobmodel.particle_sieving_stats = self.particle_sieving_stats
            jobmodel.save()     

    def terminate_classification_2D(self):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None and jobmodel.classification_2D_status == "running":
            jobmodel.classification_2D_status     = "terminating"
            jobmodel.classification_2D_update     = {"terminate":True}
            self.classification_2D_stats["stage"] = "terminating"
            jobmodel.classification_2D_stats = self.classification_2D_stats
            jobmodel.save()  

    def restart_optics(self, project, dataset):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            jobmodel.optics_assignment_status    = "restarting"
            jobmodel.optics_assignment_heartbeat = timezone.now()
            jobmodel.optics_assignment_stats     = {}
            jobmodel.optics_assignment_update    = {}
            simplestream = SIMPLEStream()
            if not simplestream.restart(os.path.join(project.dirc, dataset.dirc, self.dirc),  "optics_assignment"):
                return
            print("SAVE")
            jobmodel.save()    
    
    def restart_preprocess(self, project, dataset):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            jobmodel.preprocessing_status    = "restarting"
            jobmodel.preprocessing_heartbeat = timezone.now()
            jobmodel.preprocessing_stats     = {}
            jobmodel.preprocessing_update    = {}
            simplestream = SIMPLEStream()
            if not simplestream.restart(os.path.join(project.dirc, dataset.dirc, self.dirc),  "preprocessing"):
                return
            jobmodel.save()    
            
    def restart_initial_pick(self, project, dataset):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            jobmodel.initial_picking_status    = "restarting"
            jobmodel.initial_picking_heartbeat = timezone.now()
            jobmodel.initial_picking_stats     = {}
            jobmodel.initial_picking_update    = {}
            jobmodel.generate_pickrefs_status    = "restarting"
            jobmodel.generate_pickrefs_heartbeat = timezone.now()
            jobmodel.generate_pickrefs_stats     = {}
            jobmodel.generate_pickrefs_update    = {}
            simplestream = SIMPLEStream()
            if not simplestream.restart(os.path.join(project.dirc, dataset.dirc, self.dirc),  "opening_2D"):
                return
            jobmodel.save()
    
    def restart_generate_pickrefs(self, project, dataset):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            jobmodel.generate_pickrefs_status    = "restarting"
            jobmodel.generate_pickrefs_heartbeat = timezone.now()
            jobmodel.generate_pickrefs_stats     = {}
            jobmodel.generate_pickrefs_update    = {}
            jobmodel.initial_picking_status    = "restarting"
            jobmodel.initial_picking_heartbeat = timezone.now()
            jobmodel.initial_picking_stats     = {}
            jobmodel.initial_picking_update    = {}
            simplestream = SIMPLEStream()
            if not simplestream.restart(os.path.join(project.dirc, dataset.dirc, self.dirc),  "opening_2D"):
                return
            jobmodel.save()    
            
    def restart_reference_picking(self, project, dataset):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            jobmodel.reference_picking_status    = "restarting"
            jobmodel.reference_picking_heartbeat = timezone.now()
            jobmodel.reference_picking_stats     = {}
            jobmodel.reference_picking_update    = {}
            simplestream = SIMPLEStream()
            if not simplestream.restart(os.path.join(project.dirc, dataset.dirc, self.dirc),  "reference_based_picking"):
                return
            jobmodel.save()
    
    def restart_sieve_particles(self, project, dataset):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            jobmodel.particle_sieving_status    = "restarting"
            jobmodel.particle_sieving_heartbeat = timezone.now()
            jobmodel.particle_sieving_stats     = {}
            jobmodel.particle_sieving_update    = {}
            simplestream = SIMPLEStream()
            if not simplestream.restart(os.path.join(project.dirc, dataset.dirc, self.dirc),  "particle_sieving"):
                return
            jobmodel.save()    
            
    def restart_classification_2D(self, project, dataset):
        jobmodel = JobModel.objects.filter(id=self.id).first()
        if jobmodel is not None:
            jobmodel.classification_2D_status    = "restarting"
            jobmodel.classification_2D_heartbeat = timezone.now()
            jobmodel.classification_2D_stats     = {}
            jobmodel.classification_2D_update    = {}
            simplestream = SIMPLEStream()
            if not simplestream.restart(os.path.join(project.dirc, dataset.dirc, self.dirc),  "classification_2D"):
                return
            jobmodel.save()