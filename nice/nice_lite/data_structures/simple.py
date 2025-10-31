# global imports
import os
import stat
import json
import copy
import shutil
import subprocess
from ..models import DispatchModel, JobModel

class SIMPLEStream:

    ui_cmd            = ["simple_private_exec", "prg=print_ui_stream"]
    ui                = {}
    base_dir          = ""
    args              = {}
    jobid             = 0
    executable        = "simple_stream"
    tplt_simple_motif = "XXXSIMPLEXXX"
    tplt_nthr_motif   = "XXXNCPUXXX"
    skip_refgen       = False

    def __init__(self):
        self.loadUIJSON()

    def loadUIJSON(self):
        """populates self.ui with simple_stream UI in JSON format by running 
        command defined in ui_cmd
        Args:
            none
        Returns:
            none
        """
        try:
            ui_str = subprocess.run(
                self.ui_cmd,
                capture_output = True,
                check          = True,
                text           = True
            )
            ui_json = json.loads(ui_str.stdout)
        except subprocess.CalledProcessError as cpe:
            print(cpe.stderr, end="")
            ui_json = {}
        self.ui = ui_json
    
    def start(self, args, base_dir, jobid):
        self.base_dir    = base_dir
        self.args        = args
        self.jobid       = jobid
        if self.base_dir == "":               return False
        if not os.path.exists(self.base_dir): return False
        if not os.path.isdir(self.base_dir):  return False
        if "processes" not in self.ui:        return False
        if "pickrefs" in self.args:
            self.skip_refgen = True
        for process in self.ui["processes"]:
            if not self.dispatch(process): return False
        return True
    
    def restart(self, base_dir, processname):
        self.base_dir = base_dir
        dispatchmodel = DispatchModel.objects.filter(active=True).last()
        if dispatchmodel is None:             return False
        if self.base_dir == "":               return False
        if not os.path.exists(self.base_dir): return False
        if not os.path.isdir(self.base_dir):  return False
        for process in self.ui["processes"]:
            if process["name"] == processname:
                dispatch_script_path = os.path.join(self.base_dir, process["name"] + ".script")
                if shutil.which(dispatchmodel.scmd) is None: return False
                submit_cmd =[dispatchmodel.scmd, dispatch_script_path, '&']
                try:
                    if dispatchmodel.scmd == 'bsub':
                        subprocess.Popen(['bsub', '-L', '/bin/sh'],
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE,
                            cwd=self.base_dir,
                            stdin=open(dispatch_script_path, 'r')
                        )
                    else:
                        subprocess.Popen(submit_cmd,
                            cwd=self.base_dir,
                            start_new_session=True
                        )
                except subprocess.CalledProcessError as cpe:
                    print(cpe.stderr, end="")
                    return False
        return True
    
    def dispatch(self, process):
        dispatchmodel = DispatchModel.objects.filter(active=True).last()
        if dispatchmodel is None:                            return False
        if self.tplt_simple_motif not in dispatchmodel.tplt: return False
        if "name" not in process:                            return False
        if "prg"  not in process:                            return False
        if self.skip_refgen and (process["name"] == "opening_2D"): 
            jobmodel = JobModel.objects.filter(id=self.jobid).first()
            jobmodel.initial_picking_stats    = {"state":"skipped"}
            jobmodel.generate_pickrefs_stats  = {"state":"skipped"}
            jobmodel.initial_picking_status   = "skipped"
            jobmodel.generate_pickrefs_status = "skipped"
            jobmodel.save()
            return True
        # set path to pickrefs for reference based picking
        if not self.skip_refgen and process["name"] == "reference_based_picking":
            self.args["pickrefs"] = "../opening_2D/selected_references.mrcs"
        if "nthr_master" in process:
            nthr_master = process["nthr_master"]
        else:
            nthr_master = 1
        command_string = self.executable
        command_string += " prg=" + process["prg"]
        if "user_inputs" in process:
            for user_input in process["user_inputs"]:
                if user_input in self.args:
                    command_string += " " + user_input + "=" + self.args[user_input]
        if "static_inputs" in process:
            for static_input in process["static_inputs"]:
                command_string += " " + static_input
        # append http communication arguments
        command_string += " niceprocid=" + str(self.jobid)  + " niceserver=" + dispatchmodel.url + "/api"
        # add output redirection to command_string
        command_string += " >> " + process["name"] + ".log 2>> " + process["name"] + ".error\n" 
        # replace motifs in dispatch script template
        dispatch_script = dispatchmodel.tplt
        dispatch_script = dispatch_script.replace(self.tplt_nthr_motif,   str(nthr_master))
        dispatch_script = dispatch_script.replace(self.tplt_simple_motif, command_string)
        # replace dos line breaks
        dispatch_script = dispatch_script.replace("\r\n", "\n")
        # write dispatch script
        dispatch_script_path = os.path.join(self.base_dir, process["name"] + ".script")
        try:
            with open(dispatch_script_path, "a") as scriptfile:
                scriptfile.write(dispatch_script)
        except IOError as error:
            print("File '%s' can not be written")
            return False
        os.chmod(dispatch_script_path, stat.S_IRWXU)
        if shutil.which(dispatchmodel.scmd) is None: return False
        submit_cmd =[dispatchmodel.scmd, dispatch_script_path, '&']
        try:
            if dispatchmodel.scmd == 'bsub':
                subprocess.Popen(['bsub', '-L', '/bin/sh'],
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE,
                    cwd=self.base_dir,
                    stdin=open(dispatch_script_path, 'r')
                )
            else:
                subprocess.Popen(submit_cmd,
                    cwd=self.base_dir,
                    start_new_session=True
                )
        except subprocess.CalledProcessError as cpe:
            print(cpe.stderr, end="")
            return False
        return True

class SIMPLE:
    
    ui_cmd      = ["simple_private_exec", "prg=print_ui_json"]
    ui          = {}
    base_dir    = ""
    args        = {}
    jobtype     = ""
    executable  = ""
    tplt_simple_motif = "XXXSIMPLEXXX"
    tplt_nthr_motif   = "XXXNCPUXXX"

    def __init__(self, pckg=None):
        if pckg is not None:
            if pckg == "simple":
                self.executable = "simple_exec"
            elif pckg == "single":
                self.executable = "single_exec"
        self.loadUIJSON()

    def loadUIJSON(self):
        """populates self.ui with simple UI in JSON format by running 
        command defined in ui_cmd
        Args:
            none
        Returns:
            none
        """
        try:
            ui_str = subprocess.run(
                self.ui_cmd,
                capture_output = True,
                check          = True,
                text           = True
            )
            ui_json = json.loads(ui_str.stdout)
        except subprocess.CalledProcessError as cpe:
            print(cpe.stderr, end="")
            ui_json = {}
        self.ui = ui_json
    
    def start(self, args, base_dir, parent_dir, jobtype, jobid, parent_proj=None):
        self.base_dir    = base_dir
        self.args        = args
        self.jobtype     = jobtype
        self.jobid       = jobid
        if parent_proj == None:
            self.parent_proj = os.path.join(parent_dir, "workspace.simple")
        else:
            self.parent_proj = parent_proj
        if self.base_dir == "":                  return False
        if not os.path.exists(self.base_dir):    return False
        if not os.path.isdir(self.base_dir):     return False
        if not os.path.exists(parent_dir):       return False
        if not os.path.isdir(parent_dir):        return False
        if not os.path.exists(self.parent_proj): return False
        return self.dispatch()

    def dispatch(self):
        dispatchmodel = DispatchModel.objects.filter(active=True).last()
        if dispatchmodel is None:                            return False
        if self.tplt_simple_motif not in dispatchmodel.tplt: return False
        nthr_master = 1
        # copy project file to job dir
        command_string = "cp -v " + self.parent_proj + " workspace.simple\n"
        # update project file
        command_string += "simple_exec prg=update_project projfile=workspace.simple\n"
        command_string += self.executable
        command_string += " prg=" + self.jobtype
        for key in self.args.keys():
            command_string += " " + key + "=" + str(self.args[key])
        # append projfile and mkdir arguments
        command_string += " mkdir=no projfile=workspace.simple"
        # append http communication arguments
        command_string += " niceprocid=" + str(self.jobid)  + " niceserver=" + dispatchmodel.url + "/api_classic"
        # add output redirection to command_string
        command_string += " >> stdout.log 2>> stderr.log\n" 
        # replace motifs in dispatch script template
        dispatch_script = dispatchmodel.tplt
        dispatch_script = dispatch_script.replace(self.tplt_nthr_motif,   str(nthr_master))
        dispatch_script = dispatch_script.replace(self.tplt_simple_motif, command_string)
        # replace dos line breaks
        dispatch_script = dispatch_script.replace("\r\n", "\n")
        # write dispatch script
        dispatch_script_path = os.path.join(self.base_dir, "job.script")
        try:
            with open(dispatch_script_path, "a") as scriptfile:
                scriptfile.write(dispatch_script)
        except IOError as error:
            print("File '%s' can not be written")
            return False
        os.chmod(dispatch_script_path, stat.S_IRWXU)
        if shutil.which(dispatchmodel.scmd) is None: return False
        submit_cmd =[dispatchmodel.scmd, dispatch_script_path, '&']
        try:
            if dispatchmodel.scmd == 'bsub':
                subprocess.Popen(['bsub', '-L', '/bin/sh'],
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE,
                    cwd=self.base_dir,
                    stdin=open(dispatch_script_path, 'r')
                )
            else:
                subprocess.Popen(submit_cmd,
                    cwd=self.base_dir,
                    start_new_session=True
                )
        except subprocess.CalledProcessError as cpe:
            print(cpe.stderr, end="")
            return False
        return True
    
class SIMPLEProject:
    
    cmd = ["simple_exec", "prg=new_project", "projname=workspace"]
    dir = None

    def __init__(self, dir):
        self.dir = dir
        
    def create(self):
        """
        Args:
            none
        Returns:
            none
        """
        if self.dir is None:
            return
        cmd = copy.deepcopy(self.cmd)
        cmd.append("dir=" + self.dir)
        try:
            ui_str = subprocess.run(
                cmd
            )

        except subprocess.CalledProcessError as cpe:
            print(cpe.stderr, end="")
        
class SIMPLEProjFile:
    
    cmd_header = ["simple_exec", "prg=print_project_info",  "json=yes"]
    cmd_field  = ["simple_exec", "prg=print_project_field", "json=yes"]
    projfile   = None
    ui         = {}

    def __init__(self, projfile):
        self.projfile = projfile
        
    def getGlobalStats(self):
        """
        Args:
            none
        Returns:
            none
        """
        if self.projfile is None:
            return {}
        if not os.path.exists(self.projfile): return False
        if not os.path.isfile(self.projfile) and not os.path.islink(self.projfile): return False
        cmd = copy.deepcopy(self.cmd_header)
        cmd.append("projfile=" + self.projfile)
        print(cmd)
        try:
            ui_str = subprocess.run(
                cmd, 
                capture_output = True,
                check          = True,
                text           = True
            )
            ui_json = json.loads(ui_str.stdout)
        except subprocess.CalledProcessError as cpe:
            print(cpe.stderr, end="")
            ui_json = {}
        self.ui = ui_json
        return self.ui

    def getFieldStats(self, oritype, fromp=None, top=None, sortkey=None, sortasc=True, hist=False, boxes=False, plotkey=None):
        """
        Args:
            none
        Returns:
            none
        """
        if oritype is None:
            return {}
        if self.projfile is None:
            return {}
        if not os.path.exists(self.projfile): return False
        if not os.path.isfile(self.projfile) and not os.path.islink(self.projfile): return False
        cmd = copy.deepcopy(self.cmd_field)
        cmd.append("oritype=" + oritype)
        cmd.append("projfile=" + self.projfile)
        if fromp is not None and top is not None:
            cmd.append("fromp=" + str(fromp))
            cmd.append("top="   + str(top))
        if sortkey is not None:   
            cmd.append("sort=" + sortkey)
        if sortasc is False:   
            cmd.append("sort_asc=no")
        if hist:   
            cmd.append("hist=yes")  
        if boxes:   
            cmd.append("boxes=yes")
        if plotkey is not None and sortkey is not None:
            cmd.append("plot_key=" + plotkey)
        print(cmd)
        try:
            ui_str = subprocess.run(
                cmd, 
                capture_output = True,
                check          = True,
                text           = True
            )
            ui_json = json.loads(ui_str.stdout)
        except subprocess.CalledProcessError as cpe:
            print(cpe.stderr, end="")
            ui_json = {}
        self.ui = ui_json
        return self.ui