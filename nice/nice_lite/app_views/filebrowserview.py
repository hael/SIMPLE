import os

from django.shortcuts import render

class FileBrowserView:

    template = "filebrowser.html"
    known_file_extensions=[ ".spi", ".mrc", ".mrcs", ".dm3", ".bin", ".gain", ".eer", ".tiff", ".tif", ".jpeg", 
                            ".jpg", ".txt", ".simple", ".dat", ".img", ".map", ".head", ".ctf", ".raw", ".sbin",
                            ".dbin", ".asc", ".box", ".dat", ".pdb", ".star", ".hdf", ".pdf", ".ps"
                          ]

    def __init__(self, request, type, path):
        self.request     = request
        self.path        = path
        self.type        = type
        self.path_isdir  = None
        self.parentdir   = None
        self.error       = False
        self.errortext   = ""
        self.files       = []
        self.dirs        = []
        
    def testPath(self):
        self.path = self.path.replace('//', '/') # deal with multiple / in path
        if self.path[0] != '/':
            self.path = '/' + self.path # fix missing leading / when using proxy
        try:
            if not os.path.exists(self.path):
                self.error     = True
                self.errortext = "path does not exist"
                return
            if os.path.isdir(self.path):
                self.path_isdir = True
                self.parentdir  = os.path.dirname(self.path)
            else:
                self.path_isdir = False
        except OSError:
            self.error     = True
            self.errortext = "error"

    def listPath(self):
        try:
            contents = os.listdir(self.path)
            for entry in contents:
                # ignore hidden files/folders
                if entry[0] == '.':
                    continue
                ext = os.path.splitext(entry)[1].lower()
                if ext in self.known_file_extensions:
                    self.files.append(entry)
                elif os.path.isdir(os.path.join(self.path, entry)):
                    self.dirs.append(entry)
                else:
                    self.files.append(entry)
        except OSError:
            self.error     = True
            self.errortext = "permission denied"
        
    def render(self):
        self.testPath()
        if self.path_isdir:
            self.listPath()
        self.files.sort()
        self.dirs.sort()
        context = {
            "type"      : self.type,
            "path"      : self.path,
            "parentdir" : self.parentdir,
            "error"     : self.error,
            "errortext" : self.errortext,
            "files"     : self.files,
            "dirs"      : self.dirs
        }
        response = render(self.request, self.template, context)
        return response