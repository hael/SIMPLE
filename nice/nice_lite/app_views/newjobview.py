from django.shortcuts import render

from ..data_structures.project import Project
from ..data_structures.simple  import SIMPLE

class NewJobView:

    template = "nice_classic/newjob.html"

    def __init__(self, request, parentid, package, jobtype, args=None):
        self.request  = request
        self.parentid = parentid
        self.jobtype  = jobtype
        self.package  = package
        self.args     = args

    def render(self):
        simple   = SIMPLE()
        sections = []
        if simple.ui is not None:
          if self.jobtype in simple.ui:
              for section in simple.ui[self.jobtype]:
                  if section != "program" and len(simple.ui[self.jobtype][section]) > 0:
                      if self.args is not None:
                        for input in simple.ui[self.jobtype][section]:
                            if input["key"] in self.args:
                               input["value"] = self.args[input["key"]]
                      sections.append({
                          "name"   : section,
                          "inputs" : simple.ui[self.jobtype][section]
                      })
        context = {
            "parentid" : self.parentid,
            "jobtype"  : self.jobtype,
            "jobdesc"  : self.jobtype.replace("_", " "),
            "package"  : self.package,
            "sections" : sections
        }           
        response = render(self.request, self.template, context)
        for cookie in  self.request.COOKIES:
            if "checksum" in cookie:
                response.delete_cookie(key=cookie)
        return response
  
class NewJobTypeView:

    template = "nice_classic/newjobtype.html"

    def __init__(self, request, parentid):
        self.request  = request
        self.parentid = parentid
        self.simple   = SIMPLE()
        
    def render(self):
        simple_programs = []
        if self.simple.ui is not None:
          for prg in self.simple.ui.keys():
              if "program" in self.simple.ui[prg] and "executable" in self.simple.ui[prg]["program"]:
                if self.simple.ui[prg]["program"]["executable"] == "simple_exec" or self.simple.ui[prg]["program"]["executable"] == "all":
                    simple_programs.append({
                        "prg"  : prg,
                        "disp" : prg.replace("_", " "),
                        "desc" : self.simple.ui[prg]["program"]["descr_short"]
                    })
        single_programs = []
        if self.simple.ui is not None:
          for prg in self.simple.ui.keys():
              if "program" in self.simple.ui[prg] and "executable" in self.simple.ui[prg]["program"]:
                if self.simple.ui[prg]["program"]["executable"] == "single_exec" or self.simple.ui[prg]["program"]["executable"] == "all":
                    single_programs.append({
                        "prg"  : prg,
                        "disp" : prg.replace("_", " "),
                        "desc" : self.simple.ui[prg]["program"]["descr_short"]
                    })           
        context = {
            "simple_programs" : simple_programs,
            "single_programs" : single_programs,
            "parentid"        : self.parentid,   
        }
        response = render(self.request, self.template, context)
        for cookie in  self.request.COOKIES:
            if "checksum" in cookie:
                response.delete_cookie(key=cookie)
        return response