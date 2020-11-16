const fs = require(global.base + '/node_modules/fs-extra')
const spawn = require(global.base + '/node_modules/child-process-promise').spawn
const path = require('path')
const grepit = require(global.base + '/node_modules/grepit')

const sqlite = require('./sqlite')

class SimpleExec {

  constructor(){}
  
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
  
  getProjectInfo(file){
    return fs.stat(file)
    .then(filestat => {
      if(filestat.isFile()){
        var command = "simple_exec"
        var commandargs = ["prg=print_project_info", "projfile=" + file]
        return spawn(command, commandargs, {capture: ['stdout']})
      }else{
        return
      }
    })
  }

  getProjectField(file, oritype){
    return fs.stat(file)
    .then(filestat => {
      if(filestat.isFile()){
        var command = "simple_exec"
        var commandargs = ["prg=print_project_field", "projfile=" + file, "oritype=" + oritype]
        return spawn(command, commandargs, {capture: ['stdout']})
      }else{
        return
      }
    })
  }

  getProjectVals(file, oritype, keys){  
    return fs.stat(file)
    .then(filestat => {
      if(filestat.isFile()){
        var command = "simple_private_exec"
        var commandargs = [
          "prg=print_project_vals",
          "projfile=" + file,
          "oritype=" + oritype,
          "keys=" + keys,
        ]
        return spawn(command, commandargs, {capture: ['stdout']})
      }else{
        return
      }
    })
  }
  
  valsToArray(text){
	  return new Promise((resolve, reject) => {
		  var returnarray = []
		  for(var line of text.split('\n')){
			  var arrayelement = {}
			  for(var element of line.split(" ")){
				  if(element.includes("=")){
					  var keyval = element.split("=")
					  arrayelement[keyval[0]] = keyval[1]
				  }
			  }
			  if(Object.keys(arrayelement).length > 0){
				returnarray.push(arrayelement)
			  }
		  }
		  resolve(returnarray)
	  })
  }
  
  createProject(name, folder){
	  return(spawn("simple_exec", ["prg=new_project", "projname=" + name], {cwd: folder, capture: ['stdout']}))
  }
  
  createCavgs(projfile){
	  console.log(path.dirname(projfile) + '/selected_clusters.mrc')
	  return(spawn("simple_private_exec", ["prg=export_cavgs", "projfile=" + projfile, "outstk=" + path.dirname(projfile) + '/selected_clusters.mrc'], {cwd: path.dirname(projfile), capture: ['stdout']}))
  }
  
  relionExport(projfile){
	var dir = path.dirname(projfile) + '_RELION'
	var micrographspresent = false
	var particlespresent = false
	var mics = []
	var stks = []
	var ptcls2D = []
	
	return fs.ensureDir(dir)
	.then(() => {
	  return this.getProjectInfo(projfile)
	})
	.then(projinfo => {
		micrographspresent = (projinfo.stdout.includes('mic')) ? true : false
		particlespresent = (projinfo.stdout.includes('ptcl2D')) ? true : false
		return
	})
	.then(() => {
		if(micrographspresent){
			return fs.ensureDir(dir + '/micrographs')
			.then(() => {
				return this.getProjectField(projfile, 'mic')
			})
			.then(projvals => {
				return this.valsToArray(projvals.stdout)
			})
			.then(projvalarray => {
				mics = projvalarray
				var promises = []
				for(var mic of mics){
					promises.push(fs.symlink(this.replaceRelativePath(path.dirname(projfile), mic['intg']), dir + '/micrographs/' + path.basename(mic['intg'])))
				}
				return Promise.all(promises)
			})
		}else{
			return
		}
	})
	.then(() => {
		if(particlespresent){
			return fs.ensureDir(dir + '/particles')
			.then(() => {
				return this.getProjectField(projfile, 'stk')
			})
			.then(projvals => {
				return this.valsToArray(projvals.stdout)
			})
			.then(projvalarray => {
				stks = projvalarray
				var promises = []
				for(var stk of stks){
					promises.push(fs.symlink(this.replaceRelativePath(path.dirname(projfile), stk['stk']), dir + '/particles/' + path.basename(stk['stk']) + 's'))
				}
				return Promise.all(promises)
			})
			.then(() => {
				return this.getProjectField(projfile, 'ptcl2D')
			})
			.then(projvals => {
				return this.valsToArray(projvals.stdout)
			})
			.then(projvalarray => {
				ptcls2D = projvalarray
				return
			})
			.then(() => {
				for(var i in stks){
					stks[i]['stackcount'] = 0
				}
				return
			})
		}else{
			return
		}
	})
	.then(() => {
		var attributes = []
		var starfile = ""
        starfile += "data_\n"
        starfile += "loop_\n"
        
        if(micrographspresent){
			starfile += "_rlnMicrographName\n"
			attributes.push(['mic', 'intg'])
		}
		
		if(particlespresent){
			starfile += "_rlnImageName\n"
			attributes.push(['stk', 'stk'])
		}
		
		if(mics[0] != undefined && mics[0]['dfx'] != undefined){
			starfile += "_rlnDefocusU\n"
			attributes.push(['mic', 'dfx'])
		}else if(stks[0] != undefined && stks[0]['dfx'] != undefined){
			starfile += "_rlnDefocusU\n"
			attributes.push(['stk', 'dfx'])
		}else if(ptcls2D[0] != undefined && ptcls2D[0]['dfx'] != undefined){
			starfile += "_rlnDefocusU\n"
			attributes.push(['ptcl2D', 'dfx'])
		}
		
		if(mics[0] != undefined && mics[0]['dfy'] != undefined){
			starfile += "_rlnDefocusV\n"
			attributes.push(['mic', 'dfy'])
		}else if(stks[0] != undefined && stks[0]['dfy'] != undefined){
			starfile += "_rlnDefocusV\n"
			attributes.push(['stk', 'dfy'])
		}else if(ptcls2D[0] != undefined && ptcls2D[0]['dfx'] != undefined){
			starfile += "_rlnDefocusV\n"
			attributes.push(['ptcl2D', 'dfy'])
		}
	
		if(mics[0] != undefined && mics[0]['angast'] != undefined){
			starfile += "_rlnDefocusAngle\n"
			attributes.push(['mic', 'angast'])
		}else if(stks[0] != undefined && stks[0]['angast'] != undefined){
			starfile += "_rlnDefocusAngle\n"
			attributes.push(['stk', 'angast'])
		}else if(ptcls2D[0] != undefined && ptcls2D[0]['angast'] != undefined){
			starfile += "_rlnDefocusAngle\n"
			attributes.push(['ptcl2D', 'angast'])
		}
		
		if(mics[0] != undefined && mics[0]['kv'] != undefined){
			starfile += "_rlnVoltage\n"
			attributes.push(['mic', 'kv'])
		}else if(stks[0] != undefined && stks[0]['kv'] != undefined){
			starfile += "_rlnVoltage\n"
			attributes.push(['stk', 'kv'])
		}else if(ptcls2D[0] != undefined && ptcls2D[0]['kv'] != undefined){
			starfile += "_rlnVoltage\n"
			attributes.push(['ptcl2D', 'kv'])
		}
		
		if(mics[0] != undefined && mics[0]['cs'] != undefined){
			starfile += "_rlnSphericalAberration\n"
			attributes.push(['mic', 'cs'])
		}else if(stks[0] != undefined && stks[0]['cs'] != undefined){
			starfile += "_rlnSphericalAberration\n"
			attributes.push(['stk', 'cs'])
		}else if(ptcls2D[0] != undefined && ptcls2D[0]['cs'] != undefined){
			starfile += "_rlnSphericalAberration\n"
			attributes.push(['ptcl2D', 'cs'])
		}
		
		if(mics[0] != undefined && mics[0]['fraca'] != undefined){
			starfile += "_rlnAmplitudeContrast\n"
			attributes.push(['mic', 'fraca'])
		}else if(stks[0] != undefined && stks[0]['fraca'] != undefined){
			starfile += "_rlnAmplitudeContrast\n"
			attributes.push(['stk', 'fraca'])
		}else if(ptcls2D[0] != undefined && ptcls2D[0]['fraca'] != undefined){
			starfile += "_rlnAmplitudeContrast\n"
			attributes.push(['ptcl2D', 'fraca'])
		}
		
		if(mics[0] != undefined && mics[0]['smpd'] != undefined){
			starfile += "_rlnMagnification\n"
			starfile += "_rlnDetectorPixelSize\n"
			attributes.push(['mic', 'smpd'])
			attributes.push(['mic', 'mag'])
		}else if(stks[0] != undefined && stks[0]['smpd'] != undefined){
			starfile += "_rlnMagnification\n"
			starfile += "_rlnDetectorPixelSize\n"
			attributes.push(['stk', 'smpd'])
			attributes.push(['stk', 'mag'])
		}else if(ptcls2D[0] != undefined && ptcls2D[0]['smpd'] != undefined){
			starfile += "_rlnMagnification\n"
			starfile += "_rlnDetectorPixelSize\n"
			attributes.push(['ptcl2D', 'mag'])
			attributes.push(['ptcl2D', 'smpd'])
		}
		
		if(ptcls2D[0] != undefined && ptcls2D[0]['xpos'] != undefined){
			starfile += "_rlnCoordinateX\n"
			attributes.push(['ptcl2D', 'xpos'])
		}
		
		if(ptcls2D[0] != undefined && ptcls2D[0]['ypos'] != undefined){
			starfile += "_rlnCoordinateY\n"
			attributes.push(['ptcl2D', 'ypos'])
		}
		
		var filename
		if(micrographspresent == true && particlespresent == false){
			filename = 'micrographs.star'
			for(var mic of mics){
				if(Number(mic['state']) > 0){
					for(var attribute of attributes){
						var rawvalue
						if(attribute[0] == 'mic'){
							rawvalue = 'micrographs/' + path.basename(mic[attribute[1]])
						}
						if(attribute[1] == 'dfx' || attribute[1] == 'dfy'){
							rawvalue = Number(rawvalue) * 10000
						}else if(attribute[1] == 'mag'){
							rawvalue = 10000
						}
						starfile += rawvalue + ' '
					}
					starfile += '\n'
				}
			}
		}else if (particlespresent == true && micrographspresent == true) {
			filename = 'particles.star'
			for (var ptclind in ptcls2D){
				var stkind = Number(ptcls2D[ptclind]['stkind']) - 1
				stks[stkind]['stackcount']++
				var micrograph 
				var base = path.basename(stks[stkind]['stk']).replace("ptcls_from_", "")
				
				for(var mic of mics){
					if(mic['intg'].includes(base)){
						micrograph = mic
						break
					}
				}
				
				if(Number(ptcls2D[ptclind]['state']) > 0){
					for(var attribute of attributes){
						var rawvalue
						if(attribute[0] == 'stk'){
							rawvalue = stks[stkind][attribute[1]]
						}else if (attribute[0] == 'mic'){
							rawvalue = mic[attribute[1]]
						}else if (attribute[0] == 'ptcl2D'){
							rawvalue = ptcls2D[ptclind][attribute[1]]
						}else if (attribute[0] == 'mic'){
							rawvalue = micrograph[attribute[1]]
						}
						if(attribute[1] == 'stk'){
							rawvalue = stks[stkind]['stackcount'] +'@particles/' + path.basename(rawvalue) + 's'
						}else if(attribute[1] == 'intg'){
							rawvalue = 'micrographs/' + path.basename(rawvalue)
						}else if(attribute[1] == 'dfx' || attribute[1] == 'dfy'){
							rawvalue = Number(rawvalue) * 10000
						}else if(attribute[1] == 'mag'){
							rawvalue = 10000
						}else if (attribute[1] == 'xpos' || attribute[1] == 'ypos'){
							rawvalue = Number(rawvalue) + (Number(stks[stkind]['box']) / 2)
						}
						
						starfile += rawvalue + ' '
					}
					starfile += '\n'
				}
			}
		}else if (particlespresent == true && micrographspresent == false) {
			filename = 'particles.star'
			for (var ptclind in ptcls2D){
				var stkind = Number(ptcls2D[ptclind]['stkind']) - 1
				stks[stkind]['stackcount']++
				if(Number(ptcls2D[ptclind]['state']) > 0){
					for(var attribute of attributes){
						var rawvalue
						if(attribute[0] == 'stk'){
							rawvalue = stks[stkind][attribute[1]]
						}else if (attribute[0] == 'ptcl2D'){
							rawvalue = ptcls2D[ptclind][attribute[1]]
						}
						if(attribute[1] == 'stk'){
							rawvalue = stks[stkind]['stackcount'] +'@particles/' + path.basename(rawvalue) + 's'
							
						}else if(attribute[1] == 'dfx' || attribute[1] == 'dfy'){
							rawvalue = Number(rawvalue) * 10000
						}else if(attribute[1] == 'mag'){
							rawvalue = 10000
						}
						
						starfile += rawvalue + ' '
					}
					starfile += '\n'
				}
			}
		}
		return fs.writeFile(dir + '/' + filename, starfile, 'utf8')
	})
	.catch(err => {
	  console.log(err)
	})
  }

  getCommandArgs(arg){
		return new Promise ((resolve, reject) => {
			var commandargs = ["prg=" + arg['type']]
			var environmentargs = ['prg=update_project']

			for(var key of Object.keys(arg['keys'])){
			  if(arg['keys'][key]!= "" && !key.includes('keyenv')){
				if(arg['keys'][key] == "true"){
				  commandargs.push(key.replace('key', '') + "=yes")
				}else if(arg['keys'][key] == "false"){
				  commandargs.push(key.replace('key', '') + "=no")
				}else{
				  commandargs.push(key.replace('key', '') + "=" + arg['keys'][key])
				}
			  } else if (arg['keys'][key] != "" && key.includes('keyenv')){
				environmentargs.push(key.replace('keyenv', '') + "=" + arg['keys'][key])
			  }
			}
			
			if(arg['projfile']){
				commandargs.push("projfile=" + arg['projfile'])
				environmentargs.push("projfile=" + arg['projfile'])
			}
			resolve([commandargs, environmentargs])
		})
	}
	
	createDir(arg){
		var commandargs
		this.jobid = false
		return this.getCommandArgs(arg)
		.then(commandarguments => {
			return sqlite.sqlQuery("INSERT into " + arg['projecttable'] + " (name, description, arguments, status, view, type, parent, folder) VALUES ('" + arg['name'] + "','" + arg['description'] + "','" + JSON.stringify(arg) + "','running', '" + JSON.stringify(arg['view']) + "', '" + arg['type'] + "', '" + arg['projfile'] + "', 'null')")
		})
		.then(rows => {
			return sqlite.sqlQuery("SELECT seq FROM sqlite_sequence WHERE name='" + arg['projecttable'] + "'")
		})
		.then(rows => {
			this.jobid = rows[0]['seq']
			console.log('JOBID', this.jobid)
			console.log(arg['projectfolder'] + '/' + this.jobid + '_' + arg['type'])
			return fs.mkdir(arg['projectfolder'] + '/' + this.jobid + '_' + arg['type'])
		})
		.then(() => {
			return sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET folder='" + arg['projectfolder'] + "/" + this.jobid + '_' + arg['type'] + "' WHERE id=" + this.jobid)
		})
		.then(() => {	
			return sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET status='Finished' WHERE id=" + this.jobid)
		})
		.then(() => {
			return fs.copyFile(arg['projfile'], arg['projectfolder'] + '/' + this.jobid + '_' + arg['type'] + '/' + arg['projfile'].split('/').pop())
		})
	}

	
	exec(arg){
		this.execdir = false
		this.jobid = false
		this.buffer = false
		this.progress = 0
		this.iteration = 0
		var commandargs
		
		return sqlite.sqlQuery("INSERT into " + arg['projecttable'] + " (name, description, arguments, status, view, type, parent) VALUES ('" + arg['name'] + "','" + arg['description'] + "','" + JSON.stringify(arg) + "','pending', '" + JSON.stringify(arg['view']) + "', '" + arg['type'] + "', '" + arg['projfile'] + "')")
		.then(rows => {
			return sqlite.sqlQuery("SELECT seq FROM sqlite_sequence WHERE name='" + arg['projecttable'] + "'")
		})
		.then(rows => {
			this.jobid = rows[0]['seq']
			console.log('JOBID', this.jobid)
		})
		.then(() => {
			return this.getCommandArgs(arg)
		})
		.then(commandarguments => {
			commandargs = commandarguments
			return spawn("simple_exec", commandargs[1], {cwd: arg['projectfolder']})
		})
		.then(output => {
			var promise = spawn("simple_exec", commandargs[0], {cwd: arg['projectfolder']})
			var childprocess = promise.childProcess
			childprocess.stdout.on('data', data => {
				var lines = data.toString().split("\n")
				if(!this.execdir){
					for (var line of lines){
						if(line.includes("EXECUTION DIRECTORY")){
							this.execdir = arg['projectfolder'] + "/" + line.split(" ").pop()
							this.startProgressWatcher(this.execdir, arg['projfile'], arg['projecttable'], this.jobid)
							console.log('PID', childprocess.pid)
							sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET folder='" + this.execdir  + "', pid='" + childprocess.pid + "', status='running' WHERE id=" + this.jobid)
						//	sqlite.sqlQuery("INSERT into " + arg['projecttable'] + " (name, description, arguments, status, view, type, parent, folder, pid) VALUES ('" + arg['name'] + "','" + arg['description'] + "','" + JSON.stringify(arg) + "','running', '" + JSON.stringify(arg['view']) + "', '" + arg['type'] + "', '" + arg['projfile'] + "', '" + this.execdir + "', '" + childprocess.pid + "')")
							break
						}
					}
				}
				this.updateProgress(lines, arg)
				if(this.execdir){
					if(this.buffer != false){
						fs.appendFile(this.execdir + '/simple.log', this.buffer)
						this.buffer = false
					}
					fs.appendFile(this.execdir + '/simple.log', data.toString())
				}else{
					this.buffer = this.buffer + data.toString()
				}
			})
			console.log(`Spawned child pid: ${childprocess.pid}`)
			return promise
		})
		.then(output =>{
			if(output.childProcess.exitCode == 0 && this.execdir && this.jobid){
				if(arg['saveclusters']){
					console.log('clusters')
					return sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET status='Finished' WHERE id=" + this.jobid)
					.then(() => {
						return this.createCavgs(this.execdir + '/' + path.basename(arg['projfile']))
					})
				}else{
					return sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET status='Finished' WHERE id=" + this.jobid)
				}
			}else{
				this.execdir = arg['projectfolder'] + '/' + this.jobid + '_' + arg['type']
				return fs.ensureDir(this.execdir)
				.then(() => {
					return fs.appendFile(this.execdir + '/simple.log', this.buffer)
				})
				.then(() => {	
					return sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET status='Error', folder='" + this.execdir + "' WHERE id=" + this.jobid)
				})
			}
		})
		.catch((err) => {
			this.execdir = arg['projectfolder'] + '/' + this.jobid + '_' + arg['type']
			return fs.ensureDir(this.execdir)
			.then(() => {
				return fs.appendFile(this.execdir + '/simple.log', err)
			})
			.then(() => {
				if(this.buffer != false){
					return fs.appendFile(this.execdir + '/simple.log', this.buffer)
				}else{
					return
				}
			})
			.then(() => {	
				return sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET status='Error', folder='" + this.execdir + "' WHERE id=" + this.jobid)
			})
		})
	}
	
	distrExec(arg){
		this.execdir = false
		this.jobid = false
		this.buffer = false
		this.progress = 0
		this.iteration = 0
		var commandargs
		
		return sqlite.sqlQuery("INSERT into " + arg['projecttable'] + " (name, description, arguments, status, view, type, parent) VALUES ('" + arg['name'] + "','" + arg['description'] + "','" + JSON.stringify(arg) + "','pending', '" + JSON.stringify(arg['view']) + "', '" + arg['type'] + "', '" + arg['projfile'] + "')")
		.then(rows => {
			return sqlite.sqlQuery("SELECT seq FROM sqlite_sequence WHERE name='" + arg['projecttable'] + "'")
		})
		.then(rows => {
			this.jobid = rows[0]['seq']
			console.log('JOBID', this.jobid)
		})
		.then(() => {
			return this.getCommandArgs(arg)
		})
		.then(commandarguments => {
			commandargs = commandarguments
			return spawn("simple_exec", commandargs[1], {cwd: arg['projectfolder']})
		})
		.then(output => {
			var promise = spawn("simple_distr_exec", commandargs[0], {cwd: arg['projectfolder']})
			var childprocess = promise.childProcess
			childprocess.stdout.on('data', data => {
				var lines = data.toString().split("\n")
				if(!this.execdir){
					for (var line of lines){
						if(line.includes("EXECUTION DIRECTORY")){
							this.execdir = arg['projectfolder'] + "/" + line.split(" ").pop()
							this.startProgressWatcher(this.execdir, arg['projfile'], arg['projecttable'], this.jobid)
							console.log('PID', childprocess.pid)
							sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET folder='" + this.execdir  + "', pid='" + childprocess.pid + "', status='running' WHERE id=" + this.jobid)
						//	sqlite.sqlQuery("INSERT into " + arg['projecttable'] + " (name, description, arguments, status, view, type, parent, folder, pid) VALUES ('" + arg['name'] + "','" + arg['description'] + "','" + JSON.stringify(arg) + "','running', '" + JSON.stringify(arg['view']) + "', '" + arg['type'] + "', '" + arg['projfile'] + "', '" + this.execdir + "', '" + childprocess.pid + "')")
							//this.startProgressWatcher(arg['projfile'], )
							break
						}
					}
				}
				this.updateProgress(lines, arg)
				if(this.execdir){
					if(this.buffer != false){
						fs.appendFile(this.execdir + '/simple.log', this.buffer)
						this.buffer = false
					}
					fs.appendFile(this.execdir + '/simple.log', data.toString())
				}else{
					this.buffer = this.buffer + data.toString()
				}
			})
			console.log(`Spawned child pid: ${childprocess.pid}`)
			return promise
		})
		.then(output =>{
			
			if(output.childProcess.exitCode == 0 && this.execdir && this.jobid){
				return sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET status='Finished' WHERE id=" + this.jobid)
			}else{
				this.execdir = arg['projectfolder'] + '/' + this.jobid + '_' + arg['type']
				return fs.ensureDir(this.execdir)
				.then(() => {
					return fs.appendFile(this.execdir + '/simple.log', this.buffer)
				})
				.then(() => {	
					return sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET status='Error', folder='" + this.execdir + "' WHERE id=" + this.jobid)
				})
			}
		})
		.catch((err) => {
			this.execdir = arg['projectfolder'] + '/' + this.jobid + '_' + arg['type']
			return fs.ensureDir(this.execdir)
			.then(() => {
				return fs.appendFile(this.execdir + '/simple.log', err)
			})
			.then(() => {
				if(this.buffer != false){
					return fs.appendFile(this.execdir + '/simple.log', this.buffer)
				}else{
					return
				}
			})
			.then(() => {	
				return sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET status='Error', folder='" + this.execdir + "' WHERE id=" + this.jobid)
			})
		})
	}
	
	startProgressWatcher(execdir, projfile, projecttable, jobid){
		var folder = path.basename(execdir);
		var count = 0
		if(folder.includes('motion_correct') || folder.includes('ctf_estimate') || folder.includes('extract') || folder.includes('preprocess') || folder.includes('pick') ){
			var promise = spawn("simple_exec", ['prg=print_project_field', 'oritype=mic', 'projfile=' + projfile ], {cwd:execdir})
			var childprocess = promise.childProcess
			childprocess.stdout.on('data', data => {
				var lines = data.toString().split("\n")
				console.log(lines)
				for (var line of lines) {
					if ( !line.includes("state=0") && line.length > 0 ){
						count++
					}
				}
			})
			promise.then(() => {
				if ( count < 50 && count > 0 ){
					this.watcher = setInterval(() => {this.progressWatch(execdir, count, projecttable, jobid)}, 10000);
				} else {
					this.watcher = setInterval(() => {this.progressWatch(execdir, count, projecttable, jobid)}, 60000);
				}
			})
		}
	}
	
	progressWatch(execdir, count, projecttable, jobid){
		console.log('watch', execdir, count)
		var files = []
		var donefiles = 0
		
		fs.readdirSync(execdir).forEach(file => {
			files.push(file);
		});
		
		if ( execdir.includes("motion_correct") ){
			for (var i = 0 ; i < files.length ; i++){
				if (files[i].includes(".eps")){
					donefiles++
				}
			}
		} else if ( execdir.includes("preprocess_stream") ){
			count = 0
			for (var i = 0 ; i < files.length ; i++){
				if (files[i].includes("preprocess") && files[i].includes(".simple")){
					count++
				}
			}
			
			files = []
			
			
			if ( fs.existsSync(execdir + "/extract") ) {
				fs.readdirSync(execdir + "/extract").forEach(file => {
					if (file.includes("ptcls_from")){
						donefiles++
					}
				});
			} else if ( fs.existsSync(execdir + "/picker") ) {
				fs.readdirSync(execdir + "/picker").forEach(file => {
					if (file.includes(".box")){
						donefiles++
					}
				});
			} else if ( fs.existsSync(execdir + "/motion_correct") ){
				fs.readdirSync(execdir + "/motion_correct").forEach(file => {
					if (file.includes(".eps")){
						donefiles++
					}
				});
			}
			
			
		} else if ( execdir.includes("ctf_estimate") || execdir.includes("preprocess") ){
			for (var i = 0 ; i < files.length ; i++){
				if (files[i].includes("diag.jpg")){
					donefiles++
				}
			}
		} else if ( execdir.includes("extract") ){
			for (var i = 0 ; i < files.length ; i++){
				if (files[i].includes("ptcls_from")){
					donefiles++
				}
			}
		} else if ( execdir.includes("pick") ){
			for (var i = 0 ; i < files.length ; i++){
				if (files[i].includes(".box")){
					donefiles++
				}
			}
		}
		
		this.progress = Math.round(100 * donefiles / count)
		
		console.log(this.progress)
		
		sqlite.sqlQuery("UPDATE " + projecttable + " SET view='" + this.progress  + "' WHERE id=" + jobid)
		
	}
	
	
	updateProgress(lines, arg){
		if(arg['type'] == 'cluster2D' || arg['type'] == 'refine3D' ){
			for (var line of lines){
				if(line.includes("SEARCH SPACE SCANNED")){
					line = line.replace(/ +(?= )/g,'');
					var scan = line.split(" ")
					var progress = Math.floor((Number(scan[6]) * Number(scan[6])) / 100)
					
					this.iteration += 1
					if(this.iteration == 1){
						this.progress = 5
					}else if(this.iteration == 2){
						this.progress = 10
					}else if(this.iteration == 3){
						this.progress = 15
					}else if(progress > this.progress){
						this.progress = progress
					}
					sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET view='" + this.progress  + "' WHERE id=" + this.jobid)
					console.log('progress', this.progress)
				}
			}
		}else if(arg['type'] == 'cleanup2D'){
			for (var line of lines){
				if(line.includes("CLASS OVERLAP")){
					line = line.replace(/ +(?= )/g,'');
					var scan = line.split(" ")
					var progress = Math.floor((Number(scan[3]) * Number(scan[3])) * 100)
	
					if(progress > this.progress){
						this.progress = progress
					}
					sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET view='" + this.progress  + "' WHERE id=" + this.jobid)
					console.log('progress', this.progress)
				}
			}
		}else if(arg['type'] == 'initial_3Dmodel'){
			for (var line of lines){
				if(line.includes("SEARCH SPACE SCANNED")){
					line = line.replace(/ +(?= )/g,'');
					var scan = line.split(" ")
					var progress = Math.floor((Number(scan[6]) * Number(scan[6])) / 100)
					this.iteration += 1
					if(this.iteration <= 30){
						this.progress += 1
					}else if(!this.probabilistic){
						progress = progress * 0.8
						if(progress > this.progress){
							this.progress = progress
						}
					}else if(this.probabilistic){
						if(progress > this.progress){
							this.progress = progress
						}
					}
					sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET view='" + this.progress  + "' WHERE id=" + this.jobid)
					console.log('progress', this.progress)
				}else if(line.includes("PROBABILISTIC REFINEMENT")){
					this.probabilistic = true
				}
			}
		}
	}
}

module.exports = new SimpleExec()
