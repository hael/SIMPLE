const pug = require(global.base + '/node_modules/pug')
const path = require('path');
const sharp = require(global.base + '/node_modules/sharp');
const fs = require(global.base + '/node_modules/fs-extra');
const {getHeader, toPixels} = require(global.base + '/node_modules/mrchandler')
//const spawn = require(global.base + '/node_modules/child-process-promise').spawn
const { spawn } = require('child_process')

const pack = require(global.simplepath + '/gui_data/server/src/DensityServer/pack/main.js')
const sqlite = require('./sqlite')
const Task = require('./task')
const simpleexec = require('./simpleExec')
const filesystem = require('./fileSystem')
const Relion = require('./relion')

class View{

  constructor(){
	this.task = new Task()
	this.relion = new Relion()
	this.viewlog = pug.compileFile(global.simplepath + '/gui_data/client/src/viewlog.pug')
    this.simpleview = pug.compileFile(global.simplepath + '/gui_data/client/src/simpleview.pug')
    this.simpleviewraw = pug.compileFile(global.simplepath + '/gui_data/client/src/simpleviewraw.pug')
    this.simpleviewcluster2D = pug.compileFile(global.simplepath + '/gui_data/client/src/simpleviewcluster2D.pug')
    this.simpleviewparticles = pug.compileFile(global.simplepath + '/gui_data/client/src/simpleviewparticles.pug')
    this.simpleviewini3d = pug.compileFile(global.simplepath + '/gui_data/client/src/simpleviewini3d.pug')
    this.simpleviewrefine3d = pug.compileFile(global.simplepath + '/gui_data/client/src/simpleviewrefine3d.pug')
    this.simpleviewpostprocess = pug.compileFile(global.simplepath + '/gui_data/client/src/simpleviewpostprocess.pug')
    this.simpleviewmanualpicking = pug.compileFile(global.simplepath + '/gui_data/client/src/simpleviewmanualpicking.pug')
    this.view2d = pug.compileFile(global.simplepath + '/gui_data/client/src/view2d.pug')
    setInterval(this.mdbCleaner, 3600000)
  }
  
  mdbCleaner(){
	console.log('Cleaning MDB files')
	var filestmp = []
	var folder = '/tmp'
    return fs.readdir(folder)
    .then(contents => {
      filestmp = contents
      var promises = []
      for (var object of contents){
        promises.push(fs.stat(folder + "/" + object))
      }
      //return Promise.all(promises)
      return Promise.all(promises.map(p => p.catch(e => e)))
    })
    .then(statarray => {
	  var now = new Date().getTime();
	  var promises = []
      for(var i in statarray){
		  if (statarray[i] instanceof Error == false){
			if(filestmp[i].includes('.mdb')){
				var now = new Date().getTime();
				var endTime = new Date(statarray[i].ctime).getTime() + 3600000;
				if(now > endTime){
					console.log('Removing ' + folder + '/' + filestmp[i])
					promises.push(fs.remove(folder + "/" + filestmp[i]))
				}
			}
		}
      }
      return Promise.all(promises)
    })
  }

  getMRCJPEG(arg){
    var json = toPixels(arg['stackfile'], arg['frame'])
    if(arg['pick']){
	return sharp(json['pixbuf'], { raw : {
	    width : json['nx'],
            height : json['ny'],
            channels : 1
        }})
        .resize(Number(arg['width']))
        .blur(3)
        .normalize()
        .jpeg()
        .toBuffer()
        .then(function (image) {
            return({image : image})
        })

    }else{
        return sharp(json['pixbuf'], { raw : {
            width : json['nx'],
            height : json['ny'],
            channels : 1
        }})
        .resize(Number(arg['width']))
        .jpeg()
        .toBuffer()
        .then(function (image) {
            return({image : image})
        })
    }
  }
  
  getJPEG(arg){
	if(arg['spriteframe'] == "0"){
		return sharp(arg['stackfile'])
		.extract({left:0, top:0, width:512, height:512})
		.resize(Number(arg['width']))
		.jpeg()
		.toBuffer()
		.then(image => {
			return({image : image})
		})
	}else if(arg['spriteframe'] == "1"){
		return sharp(arg['stackfile'])
		.extract({left:511, top:0, width:512, height:512})
		.resize(Number(arg['width']))
		.jpeg()
		.toBuffer()
		.then(image => {
			return({image : image})
		})
			
	}else if(arg['stackfile'] != undefined && arg['boxfile'] != undefined && arg['intg'] != undefined){
			var scale
			var padding
			var svg = "<svg width='512' height='512'>"
			var image
			return this.readMRCHeader(arg['intg'])
			.then(json => {
			  scale = 512 / json['header']['nx']
			  padding = Math.floor((json['header']['nx'] - json['header']['ny']) * scale / 2)
			  return fs.readFile(arg['boxfile'], {encoding : 'utf8'})
			})
			.then(boxes => {
			  var lines = boxes.split("\n")
			  for(var line of lines){
				var elements = line.split((/[ , \t]+/))
				if(elements.length > 2){
				// svg += "<rect style='fill:none;stroke:green;stroke-width:2' x='" + Number(elements[1]) * scale + "' y='" + ((Number(elements[2]) * scale) + padding) + "' height='" + Number(elements[3]) * scale + "' width='" + Number(elements[3]) * scale + "'/>"
				  svg += "<circle style='fill:none;stroke:green;stroke-width:2' cx='" + (Number(elements[1]) + Number(elements[3])/2 ) * scale +  "' cy='" + (((Number(elements[2]) + Number(elements[3])/2 )* scale) + padding) + "' r='2' />"
				}
			  }
			  svg += '</svg>'
			  return sharp(arg['stackfile'])
			})
			.then(sharpimage => {
				image = sharpimage
			  return image.metadata()
			})
			.then(metadata => {
				if(metadata.width > 1000){
					return image.extract({left:511, top:0, width:512, height:512})
					//	.overlayWith(Buffer.from(svg), {})
						.composite([{input:Buffer.from(svg)}])
						.sharpen()
						.toBuffer()
				}else{
					//	return image.overlayWith(Buffer.from(svg), {})
					return image.composite([{input:Buffer.from(svg)}])
					.sharpen()
					.toBuffer()
				}
			})
			.then(image => {
				return sharp(image)
				.resize(Number(arg['width']))
				.toBuffer()
			})	
			.then(image => {	
				return({image : image})
			})
	}else if(arg['trackfile'] != undefined){

		var trackfilebase = arg['trackfile'].replace("_intg.mrc", "_shifts")
		
		if (fs.existsSync(trackfilebase + ".jpg")){
			return sharp(trackfilebase + ".jpg")
			.jpeg()
			.toBuffer()
			.then(image => {
			  return({image : image})
			})
		} else if (! fs.existsSync(trackfilebase + ".jpg") && fs.existsSync(trackfilebase + ".eps")){
			var commandargs = ["-q", "-sDEVICE=jpeg", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-dDEVICEWIDTHPOINTS=760", "-dDEVICEHEIGHTPOINTS=760", "-sOutputFile=" + trackfilebase + ".jpg", trackfilebase + ".eps"]
			return new Promise((resolve, reject) => {

				var execprocess = spawn("gs", commandargs)
				execprocess.on('close', function(_) {
                   			return resolve();
                		});
                		execprocess.on('error', function(error) {
                    			return reject(error);
                		});
			})
			.then(() => {
				return sharp(trackfilebase + ".jpg")
				.jpeg()
				.toBuffer()
				.then(image => {
					return({image : image})
				})
			})
		} else if (! fs.existsSync(trackfilebase + ".jpg") && fs.existsSync(trackfilebase + ".pdf")){
			var commandargs = ["-q", "-sDEVICE=jpeg", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-dDEVICEWIDTHPOINTS=760", "-dDEVICEHEIGHTPOINTS=760", "-sOutputFile=" + trackfilebase + ".jpg", trackfilebase + ".pdf"]
			return new Promise((resolve, reject) => {

                                var execprocess = spawn("gs", commandargs)
                                execprocess.on('close', function(_) {
                                        return resolve();
                                });
                                execprocess.on('error', function(error) {
                                        return reject(error);
                                });
                        })
			.then(() => {
				return sharp(trackfilebase + ".jpg")
				.jpeg()
				.toBuffer()
				.then(image => {
					return({image : image})
				})
			})
		}
		
		
	}else{
		return sharp(arg['stackfile'])
		.resize(Number(arg['width']))
		.jpeg()
		.toBuffer()
		.then(image => {
		  return({image : image})
		})
	}
  }
  
  readMRCHeader(file) {
    return new Promise((resolve, reject) => {
      resolve({header : getHeader(file), file : file});
    })
  }

  getViewSimple(arg){
	return simpleexec.getProjectInfo(arg['projfile'])
	.then(projinfo => {
		var micrographs = (projinfo.includes('mic')) ? true : false
		var particles = (projinfo.includes('stk')) ? true : false
		var cls2d = (projinfo.includes('cls2D')) ? true : false
	//	var ini3d = (projinfo.stdout.includes('cls3D')) ? true : false
	//	var refine3d = (projinfo.stdout.includes('cls3D')) ? true : false
		var ini3d = false
		var refine3d = false
		var postprocess = false
		var manualpick = (arg['projfile'].includes('manualpick')) ? true : false
		var autopick = (arg['projfile'].includes('_pick')) ? true : false
		
		if(arg['projfile'].includes("_initial_3Dmodel")){
			ini3d = true
			cls2d = false
		}
		
		if(arg['projfile'].includes("_refine3D")){
			refine3d = true
			ini3d = false
			cls2d = false
		}
		
		if(arg['projfile'].includes("_postprocess")){
			refine3d = false
			ini3d = false
			cls2d = false
			postprocess = true
		}
		
		if(arg['projfile'].includes("export_relion")){
			return({projfile:arg['projfile'], jobs:this.relion.analyse(path.dirname(arg['projfile']))})
		}else{
			return ({html : this.simpleview({projfile : arg['projfile'], micrographs:micrographs, particles:particles, cls2d:cls2d, ini3d:ini3d, refine3d:refine3d, postprocess:postprocess, manualpick:manualpick, autopick:autopick})})
		}
  
    })
  }
  
  getViewSimpleRaw(arg){
	 return simpleexec.getProjectInfo(arg['projfile'])
	 .then(projinfo => {
		//rawtext =  rawtext.normalize().replace(/[^A-Za-z 0-9 \.,\?""!@#\$%\^&\*\(\)-_=\+;:<>\/\\\|\}\{\[\]`~]*/g, '').replace(/[[].[m]/gm, '')
		var rawtext =  projinfo.normalize().replace(/[^A-Za-z 0-9 \.,\?""!@#\$%\^&\*\(\)-_=\+;:<>\/\\\|\}\{\[\]`~]*/g, '').replace(/[[].[m]/gm, '\n')
		 return ({html : this.simpleviewraw({projinfo : rawtext})})
	 })
  }
  
  getViewLog(arg){
	 var loginfo
	 var subinfo
	 return fs.readFile(arg['logfile'], 'utf8')
	 .then(loginfo => {
		 var subinfo = ''
		 if(fs.pathExistsSync(arg['logfile'].replace('simple.log', 'SIMPLE_SUBPROC_OUTPUT'))){
			 subinfo = fs.readFileSync(arg['logfile'].replace('simple.log', 'SIMPLE_SUBPROC_OUTPUT'), 'utf8')
		 }
		 return ({html : this.viewlog({loginfo:loginfo, subinfo:subinfo, logfile:arg['logfile']})})
	 })
  }
  
  getViewSimpleManualPick(arg){
	var jobdir = path.dirname(arg['projfile'])
	return simpleexec.getProjectField(arg['projfile'], "mic")
	 .then(projmicinfo => {
		 return simpleexec.valsToArray(projmicinfo)
	 })
	 .then(projmicinfovals => {
		 var micrographs = []
		 var boxes
		 for(var mic in projmicinfovals){
			 if(projmicinfovals[mic]['intg'] && projmicinfovals[mic]['state'] > 0){
				 micrographs.push({path:this.replaceRelativePath(jobdir, projmicinfovals[mic]['intg']), name:path.basename(projmicinfovals[mic]['intg'], '.mrc'), xdim:projmicinfovals[mic]['xdim'], ydim:projmicinfovals[mic]['ydim']})
			 }
		 }
		 if(fs.existsSync(jobdir + '/boxes.json')){
			 boxes = fs.readJsonSync(jobdir + '/boxes.json')
		 }
		return ({html:this.simpleviewmanualpicking({micrographs:micrographs}), boxes:boxes})
	 })
  }
  
  getViewSimplePick(arg){
	var jobdir = path.dirname(arg['projfile'])
	return simpleexec.getProjectField(arg['projfile'], "mic")
	 .then(projmicinfo => {
		 return simpleexec.valsToArray(projmicinfo)
	 })
	 .then(projmicinfovals => {
		 var json = {boxes:{}}
		 var micrographs = []
		 for(var mic in projmicinfovals){
			 if(projmicinfovals[mic]['intg'] && projmicinfovals[mic]['boxfile']){
				 var coordinates = []
				 micrographs.push({path:this.replaceRelativePath(jobdir, projmicinfovals[mic]['intg']), name:path.basename(projmicinfovals[mic]['intg'], '.mrc'), xdim:projmicinfovals[mic]['xdim'], ydim:projmicinfovals[mic]['ydim']})
				 var boxes = fs.readFileSync(this.replaceRelativePath(jobdir, projmicinfovals[mic]['boxfile']), 'utf8')
				 var boxlines = boxes.split('\n')
				 for(var line of boxlines){
					 var elements = line.replace(/\s\s+/g, ' ').split(' ')
					 if(elements.length > 1){
						 coordinates.push([Number(elements[1]) + (Number(elements[3]) / 2), Number(elements[2]) + (Number(elements[3]) / 2)])
						 json.boxsize = Number(elements[3])
						 json.ptclradius =  Math.floor(Number(elements[3]) / 4)
					 }
				}
				json['boxes'][path.basename(projmicinfovals[mic]['intg'], '.mrc')] = coordinates
			 }
		 }

		return ({html:this.simpleviewmanualpicking({micrographs:micrographs}), boxes:json, })
	 })
  }
  
  saveManualPickBoxes(arg){
	var jobdir = path.dirname(arg['projfile'])
	var promises = []
	var micrographs = Object.keys(arg['boxes'])
	for(var boxfile of micrographs){
		var contents = ''
		for(var coordinates of arg['boxes'][boxfile]){
			contents += '\t'
			contents += coordinates[0] - arg['boxsize']/2
			contents += '\t'
			contents += coordinates[1] - arg['boxsize']/2
			contents += '\t' + arg['boxsize'] + '\t' + arg['boxsize'] + '\t-1\n'
		}
		promises.push(fs.writeFile(jobdir + '/' + boxfile + '.box', contents, 'utf8'))
	}
	return Promise.all(promises)
	.then(() => {
		return fs.writeJson(jobdir + '/boxes.json', arg)
	})
	.then(() => {
		return fs.appendFile(jobdir + '/simple.log', "Wrote coordinates to box files\n" )
	})
	.then(() => {
		return ({})
	})
  }
  
  getViewSimpleMicrographs(arg){
	var jobdir = path.dirname(arg['projfile'])
	return simpleexec.getProjectField(arg['projfile'], "mic")
	 .then(projmicinfo => {
		 return simpleexec.valsToArray(projmicinfo)
	 })
	 .then(projmicinfovals => {
		 for(var mic in projmicinfovals){
			 if(projmicinfovals[mic]['thumb']){
				 projmicinfovals[mic]['thumb'] = this.replaceRelativePath(jobdir, projmicinfovals[mic]['thumb'])
			}
			 if(projmicinfovals[mic]['intg']){
				 projmicinfovals[mic]['intg'] = this.replaceRelativePath(jobdir, projmicinfovals[mic]['intg'])
			 }
			 if(projmicinfovals[mic]['ctfjpg']){
				 projmicinfovals[mic]['ctfjpg'] = this.replaceRelativePath(jobdir, projmicinfovals[mic]['ctfjpg'])
			 }
			 if(projmicinfovals[mic]['boxfile']){
				 projmicinfovals[mic]['boxfile'] = this.replaceRelativePath(jobdir, projmicinfovals[mic]['boxfile'])
			 }
			 projmicinfovals[mic]['number'] = Number(mic) + 1
		 }
		return ({stats:projmicinfovals})
	 })
  }
  
  getViewSimpleCluster2D(arg){
	var jobdir = path.dirname(arg['projfile'])
	var iterations
	var wfilt = false
	return filesystem.getFolderContents(jobdir)
	.then(foldercontents => {
		var promises = []
		for(var file of foldercontents['files']){
			if(file.includes("cavgs_iter") && !file.includes("odd") && !file.includes("even") && !file.includes("ranked") && !file.includes("wfilt")){
				promises.push(this.readMRCHeader(jobdir + "/" + file))
			}
			if(file.includes("wfilt")){
				wfilt = true
			}
		}
		return Promise.all(promises)
	})
	.then(headers => {
		iterations = headers
		for (var iteration of iterations){
			iteration['name'] = path.basename(iteration['file'], ".mrc").replace("cavgs_iter", "")
		}
		return simpleexec.getProjectField(arg['projfile'], "out")
	})
	 .then(projoutinfo => {
		 return simpleexec.valsToArray(projoutinfo)
	 })
	 .then(projoutinfovals => {
		var jobfolder = path.basename(jobdir)
		var complete = false
		var selection = false
		var finalfilename
		
		if(jobfolder.includes('_selection')){
			for(var element of projoutinfovals){
				if(element['imgkind'] == 'cavg'){
					finalfilename = this.replaceRelativePath(jobdir, element['stk'])
					selection = true
				}
			}
		}else{
			for(var element of projoutinfovals){
				if(element['imgkind'] == 'cavg' && element['stk'].includes(jobfolder)){
					complete = true
					finalfilename = this.replaceRelativePath(jobdir, element['stk'])
					break
				}
			}
		}
		if(complete){
			var finalfilearray
			return this.readMRCHeader(finalfilename)
			.then(header => {
				finalfilearray = header
				return simpleexec.getProjectField(arg['projfile'], "cls2D")
			})
			.then(projcls2dinfo => {
				return simpleexec.valsToArray(projcls2dinfo)
			})
			.then(projcls2dinfovals => {
				finalfilearray['stats'] = projcls2dinfovals
				finalfilearray['name'] = "Final"
				iterations.push(finalfilearray)
				return ({html : this.simpleviewcluster2D({iterations:iterations, wfilt:wfilt})})
			})
		}else if(selection){
			var finalfilearray
			return this.readMRCHeader(finalfilename)
			.then(header => {
				finalfilearray = header
				return simpleexec.getProjectField(arg['projfile'], "cls2D")
			})
			.then(projcls2dinfo => {
				return simpleexec.valsToArray(projcls2dinfo)
			})
			.then(projcls2dinfovals => {
				finalfilearray['stats'] = projcls2dinfovals
				finalfilearray['name'] = "Selection"
				iterations.push(finalfilearray)
				return ({html : this.simpleviewcluster2D({iterations:iterations, wfilt:wfilt})})
			})
		}else{
			return ({html : this.simpleviewcluster2D({iterations:iterations, wfilt:wfilt})})
		}
	 })
  }
  
  getViewSimpleParticles(arg){
	var jobdir = path.dirname(arg['projfile'])
	var stks = []
	
	return simpleexec.getProjectField(arg['projfile'], "stk")
	 .then(projstkinfo => {
		 return simpleexec.valsToArray(projstkinfo)
	 })
	 .then(projstkinfovals => {
		for(var element of projstkinfovals){
			var stkarrayelement = {}
			stkarrayelement['name'] = path.basename(element['stk'], ".mrc")
			stkarrayelement['header'] = {nx:element['box'], ny:element['box'], nz:element['nptcls']}
			stkarrayelement['file'] = this.replaceRelativePath(jobdir, element['stk'])
			stkarrayelement['fromp'] = element['fromp']
			stks.push(stkarrayelement)
		}
		return simpleexec.getProjectField(arg['projfile'], "ptcl2D")
	 })
	 .then(projptclinfo => {
		 return simpleexec.valsToArray(projptclinfo)
	 })
	 .then(projptclinfovals => {
		var stats = []
		for(var element of projptclinfovals){
			for(var key of Object.keys(element)){
				stats.push(key)
			}
		}
		stats = stats.filter(function(elem, pos) {
			return stats.indexOf(elem) == pos;
		})
		return ({html : this.simpleviewparticles({stats:stats}), ptcls:projptclinfovals, stks:stks})
	 })
  }
  
  readResFile(file){
	  var returnobject = {}
	  returnobject['fsc'] = []
	  return fs.readFile(path.dirname(file) + "/" + path.basename(file).replace(".mrc","").replace("recvol", "resolution").replace("_pproc", "").toUpperCase(), "UTF8")
	  .then(resfile => {
		for(var line of resfile.split("\n")){
			var splitline = line.replace(/\s\s+/g, ' ').split(" ")
			if(splitline[4] == "CORRELATION:"){
				returnobject['fsc'].push({resolution : splitline[2], correlation : splitline[5]})
			}else if(splitline[3] == "FSC=0.500"){
				returnobject['fsc500'] = splitline[6]
			}else if(splitline[3] == "FSC=0.143"){
				returnobject['fsc143'] = splitline[6]
			}
		}
		return ({fsc:returnobject, file:file, name:path.basename(file).replace(".mrc", "").replace("recvol_state01_iter", "")})
	  })
	  .catch(err =>{
		  return ({file:file, name:path.basename(file).replace(".mrc", "").replace("recvol_state01_", "")})
	  })
  }
  
  getViewSimpleIni3D(arg){
	var jobdir = path.dirname(arg['projfile'])
	var iterations = []
	return filesystem.getFolderContents(jobdir)
	.then(foldercontents => {
		var promises = []
		for(var file of foldercontents['files']){
			if(file.includes("_iter") && !file.includes("odd") && !file.includes("centered") && !file.includes("even") && !file.includes("odd")){
				promises.push(this.readResFile(jobdir + "/" + file))
			}
		}
		return Promise.all(promises)
	})
	.then(stats => {
		iterations = stats
		return simpleexec.getProjectField(arg['projfile'], "out")
	})
	 .then(projoutinfo => {
		 return simpleexec.valsToArray(projoutinfo)
	 })
	 .then(projoutinfovals => {
		var jobfolder = path.basename(jobdir)
		var complete = false
		var finalfilename
		for(var element of projoutinfovals){
			if(element['imgkind'] == 'vol_cavg' && element['vol'].includes(jobfolder)){
				complete = true
				finalfilename = this.replaceRelativePath(jobdir, element['vol'])
				break
			}
		}
		if(complete){
			var lastelement = iterations.pop()
			var recfinal = Object.assign({}, lastelement)
			var recfinalpproc = Object.assign({}, lastelement)
			iterations.push(lastelement)
			recfinal['file'] = jobdir + "/rec_final.mrc"
			recfinal['name'] = "Final"
			recfinal['projections'] = jobdir + "/cavgs_reprojs.mrc"
			recfinalpproc['file'] = jobdir + "/rec_final_pproc.mrc"
			recfinalpproc['name'] = "Final pproc"
			recfinalpproc['projections'] = jobdir + "/cavgs_reprojs.mrc"
			return simpleexec.getProjectField(arg['projfile'], "cls3D")
			.then(projcls3dinfo => {
				return simpleexec.valsToArray(projcls3dinfo)
			})
			.then(projcls3dinfovals => {
				recfinal['stats'] = projcls3dinfovals
				recfinalpproc['stats'] = projcls3dinfovals
				return this.readMRCHeader(jobdir + "/cavgs_reprojs.mrc")
			})
			.then(header => {
				recfinal['header'] = header
				recfinalpproc['header'] = header
				iterations.push(recfinal)
				iterations.push(recfinalpproc)
				return({html : this.simpleviewini3d({iterations:iterations})})
			})
		}else{
			return ({html : this.simpleviewini3d({iterations:iterations})})
		}
	 })
  }
  
  getViewSimpleRefine3D(arg){
	var jobdir = path.dirname(arg['projfile'])
	return filesystem.getFolderContents(jobdir)
	.then(foldercontents => {
		var promises = []
		for(var file of foldercontents['files']){
			if(file.includes("_iter") && !file.includes("odd") && !file.includes("centered") && !file.includes("even") && !file.includes("odd")){
				promises.push(this.readResFile(jobdir + "/" + file))
			}
		}
		return Promise.all(promises)
	})
	.then(stats => {
		return ({html : this.simpleviewrefine3d({iterations:stats})})
	})
  }
  
  getViewSimplePostprocess(arg){
	var jobdir = path.dirname(arg['projfile'])
	return filesystem.getFolderContents(jobdir)
	.then(foldercontents => {
		var iterations = []
		for(var file of foldercontents['files']){
			if(file.includes("pproc.mrc")){
				iterations.push({file:jobdir + '/' + file, name:'pproc'})
			}else if(file.includes("pproc_mirr.mrc")){
				iterations.push({file:jobdir + '/' + file, name:'pproc_mirr'})
			}else if(file.includes("automask.mrc")){
				iterations.push({file:jobdir + '/' + file, name:'automask'})
			}
		}
		return ({html : this.simpleviewpostprocess({iterations:iterations})})
	})
  }
  
  saveSimpleSelection(arg){
	return simpleexec.getProjectField(arg['projfile'], arg['keys']['keyoritype'])
	.then(projinfo => {
		 return simpleexec.valsToArray(projinfo)
	 })
	.then(projinfovals => {
		var selectiontxt = ""
		var count = 0
		for(var cluster of arg['selection']){
			selectiontxt += cluster['state'] + "\n"
			count++
		}
		while(count < projinfovals.length){
			selectiontxt += 0 + "\n"
			count++
		}
		arg['selection'] = []
		return fs.writeFile(arg['projfile'].replace(".simple", "_selected.txt"), selectiontxt)
	})
	.then(() => {
		return this.task.start(arg)
	})
	.then(() => {
		return ({})
	})
  }
  
  prepareMRCVolume(arg){
	var id = Math.floor(Math.random() * 10000)
    var config = {
        input: [ { name: 'em', filename: arg['volume']}],
        isPeriodic: false,
        outputFilename: "/tmp/" + id + ".mdb",
        blockSize: 96
      }

    return(pack.default(config.input, config.blockSize, config.isPeriodic, config.outputFilename))
    .then(() => {
        return({mdb : id})
    })  
  }
  
  replaceRelativePath(dir, relpath){
	var fullpath = relpath
	var parentdir = path.dirname(dir)
	if(relpath.includes('../')){
		fullpath = fullpath.replace('..', parentdir)
	}
	else if(relpath.includes('./')){
		fullpath = fullpath.replace('.', dir)
	}
	else if(fullpath.charAt(0) != "/"){
		fullpath = dir + "/" + fullpath
	}
	return fullpath
  }
  
  view2D(arg){
	  return this.readMRCHeader(arg['filename'])
	  .then(json => {
		return({html : this.view2d({filename:arg['filename']}), header:json['header']})
	  })
  }
}

module.exports = View
