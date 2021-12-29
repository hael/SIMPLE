const fs = require(global.base + '/node_modules/fs-extra');
const pug = require(global.base + '/node_modules/pug')
const simpleexec = require('./simpleExec')
const path = require('path')
//const spawn = require(global.base + '/node_modules/child-process-promise').spawn
const { spawn } = require('child_process')

class Widget{

  constructor(){
	  this.filetabgeneratorwidget = pug.compileFile(global.simplepath + '/gui_data/client/src/filetabgeneratorwidget.pug')
	  this.guinierwidget = pug.compileFile(global.simplepath + '/gui_data/client/src/guinierplotwidget.pug')
  }

   getFileTabGeneratorWidget(arg){
    console.log(arg)
    var filter = arg['filter'] ? arg['filter'] : ""
    var folder = arg['folder'] ? arg['folder'] : "/tmp"
    var projectFolder = arg['projectfolder'] ? arg['projectfolder'] : "/tmp" 
    var view = {
      files : [],
      filter : filter,
      folder : folder,
      projectfolder : projectFolder
    }
    var files
    return fs.readdir(view['folder'])
    .then(contents => {
      files = contents
      var promises = []
      for (var object of contents){
        promises.push(fs.stat(view['folder'] + "/" + object))
      }
      return Promise.all(promises)
    })
    .then(statarray => {
      for(var i in statarray){
        if(statarray[i].isFile()){
          if(view['filter'] != "" && files[i].includes(view['filter'])){
            view['files'].push(files[i])
          }else if(view['filter'] == ""){
            view['files'].push(files[i])
          }
        }
      }
      view['files'].sort()
      return({html : this.filetabgeneratorwidget(view)})
    })
  }
  
  saveFileTabGeneratorWidget(arg){
    var selection = ""
    for (var file of arg['files']){
      selection += arg['folder'] + "/" + file + '\n'
    }

    return fs.writeFile(arg['filename'], selection)
    .then(() => {
      return({status : "success"})
    })
  }
  
  getGuinierWidget(arg){
	return simpleexec.getProjectVals(arg['file'], 'out', 'vol,smpd,imgkind')
	.then(result => {
		var volumes = []
		for(var line of result.split("\n")){
			var elements = line.split((/[ ]+/))
			if(elements.length == 7 && elements[3] != "0.0000" && elements[5] == "vol"){
				volumes.push({vol:elements[3], name : path.basename(elements[3]), smpd:elements[4]})
			}
		}
		return({html : this.guinierwidget({volumes:volumes})})
	})
  }
  
  calculateGuinier(arg){
    var command = "simple_exec"
    var commandargs = []
    var ptcls = []
    var stks = []
    commandargs.push("prg=volops")
    commandargs.push("guinier=yes")
    commandargs.push("smpd=" + arg['smpd'])
    commandargs.push("vol1=" + arg['volume'])
    commandargs.push("lp=" + arg['lp'])
    commandargs.push("hp=" + arg['hp'])
    commandargs.push("nthr=1")

    return new Promise((resolve, reject) => {
	var res = ''
	var execprocess = spawn(command, commandargs, {cwd: path.dirname(arg['projfile'])})
	
	execprocess.stdout.on('data', function(_data) {
        	try {
                        var data=new Buffer(_data,'utf-8');
                        res+=data.toString();
                } catch(error) {} // ignore
        })
        execprocess.on('close', function(_) {
               return resolve(res);
        });
        execprocess.on('error', function(error) {
               return reject(error);
        });
    
    })	
    .then((result) => {
      var lines = result.split("\n")
      var bfac
      var plot = []
      for(var line of lines){
        if(line.includes("RECIPROCAL SQUARE RES")){
          var elements = line.split((/[ ]+/))
          plot.push({x : elements[4], y : elements[7]})
        }else if(line.includes("B-FACTOR DETERMINED TO")){
          var elements = line.split((/[ ]+/))
          bfac = elements[4]
        }
      }
      return({bfac:bfac, plot:plot})
    })
  }
  
}

module.exports = Widget
