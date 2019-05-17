declare var modules
import * as fs from 'fs'

class Module {

	private view2dview
	private	viewini3dview
	private cavgsview
	private ini3diterationview
	private ptclsview2d
	private ptclsview
	private pack
	private importmoviesview
	private ctfview
	private guinierwidget
	private autopicktrainingwidget
	private makepickrefswidget
	private preprocessstreamview
	private filetabgeneratorwidget
	
	public metadata = {
			"moduletitle" :"Simple",
			"modulename" :"simple",
			"tasks" : {
			}
	}
		
	constructor(){
		process.stdout.write("Loading module - Simple ... ")
		var interfacejson = require("./simple_user_interface.json")
		for(var command of interfacejson){
			var commandkeys = Object.keys(command)
			var task = command[commandkeys[0]]
			task['pages'] = []
			for(var i = 1; i < commandkeys.length; i++){
				var page = {}
				page['title'] = commandkeys[i]
				page['keys'] = []
				for(var key of command[commandkeys[i]]){
					page['keys'].push(key)
				}
				if(page['keys'].length > 0){
					task['pages'].push(page)
				}
			}
			this.metadata['tasks'][commandkeys[0]] = task
		}
		
		// Attach compenv
		for(var simpletask of Object.keys(this.metadata.tasks)){
			for (var taskpage of this.metadata['tasks'][simpletask]['pages']){
				if(taskpage['title'] == "computer controls"){
					taskpage['keys'].push({"key": "envuser_email", "keytype": "str", "descr_short": "Your e-mail address", "descr_long": "Your e-mail address; e.g. myname@uni.edu", "descr_placeholder": "", "required": false})
					taskpage['keys'].push({"key": "envjob_memory_per_task", "keytype": "num", "descr_short": "Memory per part", "descr_long": "Memory per part; MB per part{1600}", "descr_placeholder": "", "required": false})
					taskpage['keys'].push({"key": "envqsys_partition", "keytype": "str", "descr_short": " Name of SLURM/PBS partition", "descr_long": " Name of SLURM/PBS partition; give part name", "descr_placeholder": "", "required": false})
					taskpage['keys'].push({"key": "envqsys_qos", "keytype": "str", "descr_short": "Schedule priority", "descr_long": "Schedule priority; give priority", "descr_placeholder": "", "required": false})
					taskpage['keys'].push({"key": "envqsys_reservation", "keytype": "str", "descr_short": "Name of reserved partition", "descr_long": "Name of reserved partition; give yourpart", "descr_placeholder": "", "required": false})
					taskpage['keys'].push({"key": "envtime_per_image", "keytype": "num", "descr_short": "Time per image", "descr_long": "Time per image; in seconds{100}", "descr_placeholder": "", "required": false})
					taskpage['keys'].push({"key": "envuser_account", "keytype": "str", "descr_short": "User account name in SLURM/PBS", "descr_long": "User account name in SLURM/PBS; e.g. Account084", "descr_placeholder": "", "required": false})
					taskpage['keys'].push({"key": "envuser_project", "keytype": "str", "descr_short": "User project name in SLURM/PBS", "descr_long": "User project name in SLURM/PBS; e.g. Project001", "descr_placeholder": "", "required": false})
				}
			}
		}
	
		// Attach import
		var importtask = {}
		importtask["name"]  =  "import",
		importtask["descr_short"] = "Import an external simple job into the gui",
		importtask["descr_long"] = "Import an external simple job into the gui to view its output",
		importtask["executable"] = "simple_import"
		importtask['pages'] = [] 
		var page = {}
		page['title'] = "parameter input/output"
		page['keys'] = []
		var options = []
		for(var simpletask of Object.keys(this.metadata.tasks)){
			options.push(this.metadata['tasks'][simpletask]['name'])
		}
		page['keys'].push({"key": "job_type", "keytype": "multi", "descr_short": "Job Type", "descr_long": "", "descr_placeholder": "", "required": true, "options": options})
		importtask['pages'].push(page)
		this.metadata['tasks']['import'] = importtask
		
		// Attach widgets
		for(var widgetpage of this.metadata['tasks']['postprocess']['pages']){
			for(var widgetkeys of widgetpage['keys']){
				if(widgetkeys['key'] == "bfac"){
					widgetkeys['keytype'] = "widget"
					widgetkeys['widget'] = "guinierplotwidget.view(this)"
				}
			}
			
		}
		
		for(var widgetpage of this.metadata['tasks']['pick']['pages']){
			for(var widgetkeys of widgetpage['keys']){
				if(widgetkeys['key'] == "thres"){
					widgetkeys['keytype'] = "widget"
					widgetkeys['widget'] = "autopicktrainingwidget.view(this)"
				}else if(widgetkeys['key'] == "ndev"){
					widgetkeys['keytype'] = "widget"
					widgetkeys['widget'] = "autopicktrainingwidget.view(this)"
				}else if(widgetkeys['key'] == "refs"){
					widgetkeys['keytype'] = "widget"
					widgetkeys['widget'] = "pickrefswidget.view(this)"
				}
			}
			
		}
		
		for(var widgetpage of this.metadata['tasks']['import_movies']['pages']){
			for(var widgetkeys of widgetpage['keys']){
				if(widgetkeys['key'] == "filetab"){
					widgetkeys['keytype'] = "widget"
					widgetkeys['widget'] = "filetabgeneratorwidget.view('keyfiletab')"
				}
			}
			
		}
		
		const pug = require('pug')
		this.pack = require("../../../external/DensityServer/pack/main.js")
		this.view2dview = pug.compileFile('views/simple-view2d.pug')
		this.viewini3dview = pug.compileFile('views/simple-viewini3d.pug')
		this.cavgsview = pug.compileFile('views/simple-cavgs.pug')
		this.ini3diterationview = pug.compileFile('views/simple-ini3diteration.pug')
		this.ptclsview2d = pug.compileFile('views/simple-viewptcls2d.pug')
		this.ptclsview = pug.compileFile('views/simple-viewptcls.pug')
		this.importmoviesview = pug.compileFile('views/simple-viewimportmics.pug')
		this.ctfview = pug.compileFile('views/simple-viewctf.pug')
		this.guinierwidget = pug.compileFile('views/simple-guinierplotwidget.pug')
		this.autopicktrainingwidget = pug.compileFile('views/simple-autopicktrainingwidget.pug')
		this.makepickrefswidget = pug.compileFile('views/simple-makepickrefswidget.pug')
		this.preprocessstreamview = pug.compileFile('views/simple-viewpreprocessstream.pug')
		this.filetabgeneratorwidget = pug.compileFile('views/simple-filetabgeneratorwidget.pug')
		process.stdout.write("Done\n")
	}
	
	public refineSelection(file){
		var spawn = require('child-process-promise').spawn
		var stat = fs.statSync(file)
		if(stat.isFile()){
			if(file.includes(".simple")){
				var command = "simple_exec"
				var commandargs = []
				commandargs.push("prg=print_project_info")
				commandargs.push("projfile=" + file)
				return spawn(command, commandargs, {capture: ['stdout']})
					.then((result) => {
						var taskarray = []
						if(result.stdout.includes("micrographs")){
							taskarray.push("ctf_estimate")
							taskarray.push("motion_correct")
							taskarray.push("extract")
						}else{
							taskarray.push("import_movies")
						}
						if(result.stdout.includes("per-micrograph stack")){
							taskarray.push("cluster2D")
						}else{
							taskarray.push("import_particles")
						}
						if(result.stdout.includes("per-particle 2D")){
							taskarray.push("cluster3D")
						}else{
							taskarray.push("import_particles")
						}
						if(result.stdout.includes("per-particle 3D" )){
							taskarray.push("refine3D")
						}else{
							taskarray.push("")
						}
						if(result.stdout.includes("per-cluster 2D")){
							taskarray.push("initial_3Dmodel")
						}else{
							taskarray.push("")
						}
						return taskarray
					})
			}
		}else if(stat.isDirectory()){
			return new Promise((resolve, reject) => {
				var taskarray = []
				taskarray.push("preprocess_stream")
				resolve(taskarray)
			})
		}
	}
	
	public getCompEnv(modules, arg){
		var spawn = require('child-process-promise').spawn
		var stat = fs.statSync(arg['projfile'])
		if(stat.isFile() && arg['projfile'].includes(".simple")){
			var command = "simple_exec"
			var commandargs = []
			commandargs.push("prg=print_project_field")
			commandargs.push("projfile=" + arg['projfile'])
		}
	}
	
	public execute(modules, arg){
		var spawn = require('child_process').spawn
		var spawnPromise = require('child-process-promise').spawn
		var path = require('path')
		var command = arg['executable']
		
		var inputpath = arg['inputpath']

		var keys = arg['keys']
		var keynames = Object.keys(keys);
		var commandargs = []
		var compenvargs = ['prg=update_project', 'projname=project']
		commandargs.push("prg=" + type)
		commandargs.push("mkdir=no")
		console.log(arg)
		for(var key of keynames){
			if(keys[key] != "" && ! key.includes('keyenv')){
				if(keys[key] == "true"){
					commandargs.push(key.replace('key', '') + "=yes")
				}else if(keys[key] == "false"){
					commandargs.push(key.replace('key', '') + "=no")
				}else{
					commandargs.push(key.replace('key', '') + "=" + keys[key])
				}
			} else if (keys[key] != "" && key.includes('keyenv')){
				compenvargs.push(key.replace('keyenv', '') + "=" + keys[key])
			}
		}
		
		arg['view'] = null
		var type = arg['type']
		if(type == "import"){
			type = arg['keys']['keyjob_type']
		}
		
		if(type == "cluster2D") {
			arg['view'] = {mod : "simple", fnc : "view2d"}
		}else if(type == "initial_3Dmodel") {
			arg['view'] = {mod : "simple", fnc : "viewini3d"}
		}else if(type == "import_movies" || type == "motion_correct") {
			arg['view'] = { mod : "simple", fnc : "viewImportMovies" }
		}else if(type == "ctf_estimate") {
			arg['view'] = { mod : "simple", fnc : "viewCTFEstimate" }
		}else if(type == "import_particles" || type == "extract") {
			arg['view'] = { mod : "simple", fnc : "viewParticles" }
		}else if(type == "preprocess_stream") {
			arg['view'] = { mod : "simple", fnc : "viewPreprocessStream" }
		}else if(type == "preprocess") {
			arg['view'] = { mod : "simple", fnc : "viewPreprocessStream" }
		}
		
		console.log(command, commandargs)
		
		var JSON 
		return modules['available']['core']['taskCreate'](modules, arg)
		.then((json) => {
			JSON = json
			if(command == "simple_import"){
				var rootname = path.basename(inputpath, ".simple")
				fs.rmdirSync(json['jobfolder'])
				fs.symlinkSync(path.dirname(inputpath), json['jobfolder'])
				fs.mkdirSync(json['jobfolder'] + "/" + rootname)
				fs.copyFileSync(inputpath, json['jobfolder'] + "/" + rootname + "/" + rootname + ".simple")
				compenvargs.push('projfile=' + rootname + "/" + rootname + ".simple")
				return spawnPromise("simple_exec", compenvargs, {cwd: json['jobfolder'], capture: ['stdout']})
			} else {
				fs.mkdirSync(json['jobfolder'] + "/project")
				fs.copyFileSync(inputpath, json['jobfolder'] + "/project/project.simple")
				
				console.log(compenvargs)
				return spawnPromise("simple_exec", compenvargs, {cwd: json['jobfolder'], capture: ['stdout']})
			}
		})
		.then((result) => {
			if(command == "simple_import"){
				var rootname = path.basename(inputpath, ".simple")
				fs.renameSync(JSON['jobfolder'] + "/" + rootname + "/project.simple", JSON['jobfolder'] + "/project.simple")
				modules['available']['core']['updateStatus'](arg['projecttable'], JSON['jobid'], "Finished")
			}else{
				fs.renameSync(JSON['jobfolder'] + "/project/project.simple", JSON['jobfolder'] + "/project.simple")
				const out = fs.openSync(JSON['jobfolder'] + '/task.log', 'a')
				const err = fs.openSync(JSON['jobfolder'] + '/task.log', 'a')
				var executeprocess = spawn(command, commandargs, {detached: true, cwd: JSON['jobfolder'], stdio: [ 'ignore', out, err ]})
				executeprocess.on('exit', function(code){
				console.log(`child process exited with code ${code}`);
					if(code !== null){
						modules['available']['core']['updateStatus'](arg['projecttable'], JSON['jobid'], "Finished")
					}
				})
				executeprocess.on('error', function(error){
					console.log(`child process exited with error ${error}`)
					modules['available']['core']['updateStatus'](arg['projecttable'], JSON['jobid'], "Error")
				})
				
				console.log(`Spawned child pid: ${executeprocess.pid}`)
				modules['available']['core']['updatePid'](arg['projecttable'], JSON['jobid'], executeprocess.pid)
			}
			return({folder:JSON['jobfolder']})
		})

	}
	
	public view2d (modules, arg) {
		var contents = fs.readdirSync(arg['folder'])
		var iterations = []
		for(var file of contents){
			if(file.includes("cavgs") && file.includes(".mrc") && ! file.includes("even") && ! file.includes("odd") && ! file.includes("ranked")){
				iterations.push(file)
			}
		}
		iterations.sort()
		var view = {}
		view['iterations'] = iterations
		view['status'] = arg['status']
		view['folder'] = arg['folder']
		return new Promise((resolve, reject) => {
			resolve({html : this.view2dview(view), func : "simple.selectRecent2D()"})
		})
	}
	
	public getCavgs (mods, arg) {
		var spawn = require('child-process-promise').spawn
		var path = require('path')
		var view = {}
		if(arg['status']) {
			var projfile = path.dirname(arg['file']) + "/project.simple"
			var command = "simple_exec"
			var commandargs = []
			commandargs.push("prg=print_project_field")
			commandargs.push("projfile=" + projfile)
			commandargs.push("oritype=cls2D")
			var cavgarray = []
			return spawn(command, commandargs, {capture: ['stdout']})
				.then((result) => {
					var lines = result.stdout.split("\n")
					for(var line of lines){
						var cavg = {}
						var elements = line.split(" ")
						for(var element of elements){
							var key = element.split("=")[0]
							var value = element.split("=")[1]
							cavg[key] = value
						}
						cavgarray.push(cavg)
					}
					return 
				})
				.then (() => {
					return modules['available']['core']['readMRCHeader'](arg['file'])
				})
				.then((header) => {
					view['thumbcount'] = header['nz'] 
					view['path'] = arg['file']
					view['status'] = arg['status']
					view['cavgarray'] = cavgarray
					view['folder'] = path.dirname(arg['file'])
					view['projectfile'] = projfile
					return({html : this.cavgsview(view)})
				})
		} else {
			return modules['available']['core']['readMRCHeader'](arg['file'])
				.then((header) => {
					view['thumbcount'] = header['nz']
					view['path'] = arg['file']
					view['status'] = arg['status']
					return({html : this.cavgsview(view)})
				})
		}
	}
	
	public viewini3d (modules, arg) {
		var contents = fs.readdirSync(arg['folder'])
		var iterations = []
		for(var file of contents){
			if(file.includes("recvol_state01_iter") && file.includes("pproc")){
				iterations.push(file)
			}
		}
		iterations.sort()
		var view = {}
		view['iterations'] = iterations
		if(fs.existsSync(arg['folder'] + "/rec_final_pproc.mrc")){
			view['final'] = true
		}

		view['status'] = arg['status']
		view['folder'] = arg['folder']
		
		return new Promise((resolve, reject) => {
			resolve({html : this.viewini3dview(view)})
		})
	}
	
	public viewIni3dIteration(modules, arg) {
		var path = require("path")
		var grepit = require('grepit')
		var id = Math.floor(Math.random() * 10000)
		var config = {
				input: [ { name: 'em', filename: arg['file']}],
				isPeriodic: false,
				outputFilename: "/tmp/" + id + ".mdb",
				blockSize: 96
			}
		var view = {}
		var iterations = 0
		arg['folder'] = path.dirname(arg['file'])
		var contents = fs.readdirSync(arg['folder'])
		
		for(var file of contents){
			if(file.includes("recvol_state01_iter") && file.includes("pproc")){
				iterations++
			}
		}
		
		if(arg['file'].includes("final")){
			view['final'] = true
			if(fs.existsSync(arg['folder'] + "/RESOLUTION_STATE01_ITER0" + iterations)){
				var fsc50 = grepit(/RESOLUTION AT FSC=0.500 DETERMINED TO/, arg['folder'] + "/RESOLUTION_STATE01_ITER0" + iterations)
				var fsc143 = grepit(/RESOLUTION AT FSC=0.143 DETERMINED TO/, arg['folder'] + "/RESOLUTION_STATE01_ITER0" + iterations)
				view['fsc50'] = fsc50[0].replace('>>> RESOLUTION AT ','')
				view['fsc143'] = fsc143[0].replace('>>> RESOLUTION AT ','')
			}
			if(fs.existsSync(arg['folder'] + "/cavgs_reprojs.mrc")){
				return (modules.available.core.readMRCHeader(arg['folder'] + "/cavgs_reprojs.mrc"))
					.then((header) => {
						view['reprojections'] = arg['folder'] + "/cavgs_reprojs.mrc"
						view['reprojectionscount'] = header['nz']
						return({})
					})
					.then(() => {
						return(this.pack.default(config.input, config.blockSize, config.isPeriodic, config.outputFilename))
					})
					.then(() => {
						return({html : this.ini3diterationview(view), mdb : id})
					})
			}else{
				return(this.pack.default(config.input, config.blockSize, config.isPeriodic, config.outputFilename))
				.then(() => {
						return({html : this.ini3diterationview(view), mdb : id})
				})
			}
		} else {
			return(this.pack.default(config.input, config.blockSize, config.isPeriodic, config.outputFilename))
			.then(() => {
					return({html : this.ini3diterationview(view), mdb : id})
			})
		}
	}
	
	public setupProject(arg){
		var spawn = require('child_process').spawn
		var command = "simple_exec"
		var commandargs = []
		commandargs.push("prg=new_project")
		commandargs.push("projname=project")
		var executeprocess = spawn(command, commandargs, {detached: true, cwd: arg['keys']['keyfolder']})
		executeprocess.on('exit', function(code){
			console.log(`child process exited with code ${code}`)
			fs.renameSync(arg['keys']['keyfolder'] + "/project/project.simple", arg['keys']['keyfolder'] + "/project.simple")
			fs.rmdirSync(arg['keys']['keyfolder'] + "/project")
		})
		executeprocess.stdout.on('data', (data) => {
			console.log(`child stdout:\n${data}`)
		})
		executeprocess.stderr.on('data', (data) => {
			console.log(`child stderr:\n${data}`)
		})
		executeprocess.on('error', function(error){
			console.log(`child process exited with error ${error}`)
		})
		
		console.log(`Spawned child pid: ${executeprocess.pid}`)

	}
	
	public save2DSelection(modules, arg){
		var spawn = require('child-process-promise').spawn
		fs.copyFileSync(arg['projectfile'], arg['file'])
		var selection = ""
		
		for (var state of arg['selection']){
			selection += state + '\n'
		}

		fs.writeFileSync(arg['file'] + ".txt", selection)
		
		var command = "simple_private_exec"
		var commandargs = []
		commandargs.push("prg=update_project_stateflags")
		commandargs.push("projfile=" + arg['file'])
		commandargs.push("oritype=cls2D")
		commandargs.push("infile=" + arg['file'] + ".txt")

		return spawn(command, commandargs, {capture: ['stdout']})
			.then(() => {
				return({status : "ok"})
			})
			.catch((error) => {
				return({status : "error"})
			})
	}
	
	public saveMicsSelection(modules, arg){
		var spawn = require('child-process-promise').spawn
		fs.copyFileSync(arg['projectfile'], arg['file'])
		var selection = ""
		
		for (var state of arg['selection']){
			selection += state + '\n'
		}

		fs.writeFileSync(arg['file'] + ".txt", selection)
		
		var command = "simple_private_exec"
		var commandargs = []
		commandargs.push("prg=update_project_stateflags")
		commandargs.push("projfile=" + arg['file'])
		commandargs.push("oritype=mic")
		commandargs.push("infile=" + arg['file'] + ".txt")

		return spawn(command, commandargs, {capture: ['stdout']})
			.then(() => {
				return({status : "ok"})
			})
			.catch((error) => {
				return({status : "error"})
			})
	}
	
	private writeMicsStar(projfile){
		var spawn = require('child-process-promise').spawn
		var command = "simple_private_exec"
		var commandargs = []
		var starcontents = ""
		commandargs.push("prg=print_project_vals")
		commandargs.push("projfile=" + projfile)
		commandargs.push("oritype=mic")
		commandargs.push("keys=movie,intg,dfx,dfy,angast,ctfscore")
		
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line of lines){
					var elements = line.split((/[ ]+/))
					if(elements[2] == 1){
						starcontents += elements[3] + " "
						starcontents += elements[4] + " "
						starcontents += elements[5] + " "
						starcontents += elements[6] + " "
						starcontents += elements[7] + " "
						starcontents += elements[8] + "\n"
					}
				}
				console.log(starcontents)
				return 
			})
	}
	
	public viewParticles2D(modules, arg){
		var spawn = require('child-process-promise').spawn
		var command = "simple_private_exec"
		var commandargs = []
		var classcontents = []
		var particles = []
		var stks = []
		commandargs.push("prg=print_project_vals")
		commandargs.push("projfile=" + arg['projfile'])
		commandargs.push("oritype=ptcl2D")
	//	commandargs.push("keys=class,x,y,dfx,dfy,angast,frameid,stkind,bfac")
		commandargs.push("keys=class,x,y,dfx,dfy,angast,stkind,bfac")
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line of lines){
					var elements = line.split((/[ ]+/))
					elements.shift()
					//if(Number(elements[2]) == arg['class'] && elements[1] != "0"){
					//	classcontents.push(elements)
					//}
					particles.push(elements)
				}
				return
			})
			.then(() => {
				commandargs = []
				commandargs.push("prg=print_project_vals")
				commandargs.push("projfile=" + arg['projfile'])
				commandargs.push("oritype=stk")
				commandargs.push("keys=stk,fromp")
				return spawn(command, commandargs, {capture: ['stdout']})
			})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line of lines){
					var elements = line.split((/[ ]+/))
					stks.push([elements[3], elements[4]])
				}
				return
			})
			.then(() => {
				for(var i = 0; i < particles.length; i++){
					if(Number(particles[i][2]) == arg['class'] && particles[i][1] != "0"){
						var stk = stks[Number(particles[i][8]) - 1][0]
						var frameid = i - (Number(stks[Number(particles[i][8]) - 1][1]) - 1);
						particles[i][8] = stk
						particles[i].push(frameid)
						classcontents.push(particles[i])
					}
					
				}
				return({html : this.ptclsview2d({ptcls : classcontents})})
			})
	}
	
	public viewImportMovies(modules, arg){
		var spawn = require('child-process-promise').spawn
		var command = "simple_private_exec"
		var commandargs = []
		var mics = []
		var projectfile = arg['folder'] + "/project.simple"
		commandargs.push("prg=print_project_vals")
		commandargs.push("projfile=" + projectfile)
		commandargs.push("oritype=mic")
		commandargs.push("keys=intg")
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line of lines){
					if (line.includes("key: intg is missing in segment")){
						break
					}
					var elements = line.split((/[ ]+/))
					if(elements[3]){
						mics.push([ elements[3], elements[3].split('\\').pop().split('/').pop()])
					}
				}
				return
			})
			.then(() => {
				return({html : this.importmoviesview({movies : mics, folder : arg['folder'], projectfile : projectfile}), func : "viewer.loadImages('mic')"})
			})
	}
	
	public viewCTFEstimate(modules, arg){
		var spawn = require('child-process-promise').spawn
		var command = "simple_private_exec"
		var commandargs = []
		var fits = []
		var projectfile = arg['folder'] + "/project.simple"
		commandargs.push("prg=print_project_vals")
		commandargs.push("projfile=" + projectfile)
		commandargs.push("oritype=mic")
		commandargs.push("keys=intg,dfx,dfy,angast")
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line of lines){
					var elements = line.split((/[ ]+/))
					if(elements[3]){

						fits.push([arg['folder'] + "/" + elements[3].split('\\').pop().split('/').pop().replace(".mrc", ".jpg"), elements[4], elements[5], elements[6]])
					}
				}
				return
			})
			.then(() => {
				return({html : this.ctfview({fits : fits, folder : arg['folder'], projectfile : projectfile}), func : "viewer.loadImages('fits')"})
			})
	}
	
	public viewPreprocessStream(modules, arg){
		var spawn = require('child-process-promise').spawn
		var command = "simple_private_exec"
		var commandargs = []
		var preproc = []
		var boxfiles
		var projectfile = arg['folder'] + "/project.simple"
		commandargs.push("prg=print_project_vals")
		commandargs.push("projfile=" + projectfile)
		commandargs.push("oritype=mic")
		commandargs.push("keys=thumb,dfx,dfy,angast,ctf_estimatecc,dferr,ctfscore,cc90")
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line of lines){
					var elements = line.split((/[ ]+/))
					if(elements[3]){
						preproc.push([elements[3], elements[3].replace("motion_correct", "ctf_estimate").replace("thumb.jpg", "ctf_estimate_diag.jpg"), elements[4], elements[5], elements[6], elements[7], elements[8], elements[9], elements[10]])
					}
				}
				commandargs = []
				commandargs.push("prg=print_project_vals")
				commandargs.push("projfile=" + projectfile)
				commandargs.push("oritype=mic")
				commandargs.push("keys=boxfile,intg,xdim,ydim,thumb")
				return spawn(command, commandargs, {capture: ['stdout']})
			})
			.then((result) => {
				var lines = result.stdout.split("\n")
				if(! lines[0].includes("missing in segment")){
					boxfiles = []
					for(var line of lines){
						var elements = line.split((/[ ]+/))
						if(elements[3]){
							boxfiles.push([elements[3], elements[4], elements[5], elements[6], elements[7]])
						}
					}
				}
			})
			.then(() => {
				return({html : this.preprocessstreamview({preproc : preproc, folder : arg['folder'], projectfile : projectfile, boxfiles : boxfiles}), func : "simple.viewPreprocessStream()"})
			})
	}
	
	public viewParticles(modules, arg){
		var spawn = require('child-process-promise').spawn
		var command = "simple_private_exec"
		var commandargs = []
		var stks = []
		commandargs.push("prg=print_project_vals")
		commandargs.push("projfile=" + arg['folder'] + "/project.simple")
		commandargs.push("oritype=stk")
		commandargs.push("keys=stk,nptcls")
		
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line in lines){
					var elements = lines[line].split((/[ ]+/))
					if(elements[3]){
						stks.push([elements[3], elements[3].split('\\').pop().split('/').pop().replace(".mrc",""), elements[4]])
					}
				}
				return({html : this.ptclsview({stks : stks})})
			})
	}
	
	public viewStkParticles(modules, arg){
		var spawn = require('child-process-promise').spawn
		var command = "simple_private_exec"
		var commandargs = []
		var ptcls = []
		var stks = []
		commandargs.push("prg=print_project_vals")
		commandargs.push("projfile=" + arg['folder'] + "/project.simple")
		commandargs.push("oritype=ptcl2D")
		commandargs.push("keys=class,x,y,dfx,dfy,angast,frameid,stkind")
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line of lines){
					var elements = line.split((/[ ]+/))
					elements.shift()
					if(Number(elements[2]) == arg['class'] && elements[1] != "0"){
						ptcls.push(elements)
					}
				}
				return
			})
			.then(() => {
				commandargs = []
				commandargs.push("prg=print_project_vals")
				commandargs.push("projfile=" + arg['folder'] + "/project.simple")
				commandargs.push("oritype=stk")
				commandargs.push("keys=stk")
				return spawn(command, commandargs, {capture: ['stdout']})
			})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line of lines){
					var elements = line.split((/[ ]+/))
					if(elements[3]){
						stks.push([elements[3], elements[3].split('\\').pop().split('/').pop().replace(".mrc","")])
					}
				}
				return
			})
			.then(() => {
				for(var i = 0; i < ptcls.length; i++){
					var stk = stks[Number(ptcls[i][9]) - 1]
					ptcls[i][9] = stk
					
				}
				return({html : this.ptclsview({ptcls : ptcls, stks : stks})})
			})
	}
	
	public getGuinierWidget(modules, arg){
		return new Promise((resolve, reject) => {
			resolve({html : this.guinierwidget({})})
		})
	}
	
	public getFileTabGeneratorWidget(modules, arg){
		return new Promise((resolve, reject) => {
			var files = []
			if(arg['folder'] != undefined && arg['folder'] != ""){
				var contents = fs.readdirSync(arg['folder'])
				for (var object of contents){
					if(object.charAt(0) != "."){
						var stat = fs.statSync(arg['folder'] + "/" + object)
						if(stat.isFile()){
							if(arg['filter'] != undefined && arg['filter'] != ""){
								if(object.includes(arg['filter'])){
									files.push(object)
								}
							}else{
								files.push(object)
							}
						}
					}
				}
				files.sort()
			}
			var view = {}
			view['files'] = files
			if(arg['folder'] != undefined ){
				view['folder'] = arg['folder']
			}else{
				view['folder'] = ""
			}
			if(arg['filter'] != undefined ){
				view['filter'] = arg['filter']
			}else{
				view['filter'] = ""
			}
			resolve({html : this.filetabgeneratorwidget(view)})
		})
	}
	
	public saveFileTabGeneratorWidget(modules, arg){
		return new Promise((resolve, reject) => {
			var selection = ""
			for (var file of arg['files']){
				selection += arg['folder'] + "/" + file + '\n'
			}

			fs.writeFileSync(arg['filename'], selection)
			resolve({status : "success"})
		})
	}

	public calculateGuinier(modules, arg){
		var spawn = require('child-process-promise').spawn
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
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
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
	
	public getAutopickTrainingWidget(modules, arg){
		var spawn = require('child-process-promise').spawn
		var path = require("path")
		var command = "simple_private_exec"
		var commandargs = []
		var mics = []
		commandargs.push("prg=print_project_vals")
		commandargs.push("projfile=" + arg['projfile'])
		commandargs.push("oritype=mic")
		commandargs.push("keys=intg,xdim,ydim")
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line of lines){
					var elements = line.split((/[ ]+/))
					elements.shift()
					if(elements.length > 0){
						mics.push({path : elements[2], name : path.basename(elements[2]).replace(".mrc", ""), xdim : elements[3], ydim : elements[4]})
					}
				}
				return({html : this.autopicktrainingwidget({mics : mics})})
			})
	}
	
	public trainAutopick(modules, arg){
		var spawn = require('child-process-promise').spawn
		var command = "simple_private_exec"
		var commandargs = []
		var selection = []
		
		var trainingfolder = arg['projectfolder'] + "/.training"
		
		if(!fs.existsSync(trainingfolder)){
			fs.mkdirSync(trainingfolder)
		}
		
		var files = fs.readdirSync(trainingfolder)
		
		for (var file of files){
			if(file.includes(".box")){
				fs.unlinkSync(trainingfolder + "/" + file)
			}
		}
		
		
		fs.copyFileSync(arg['projfile'], trainingfolder + "/project.simple")
		
		commandargs.push("prg=print_project_vals")
		commandargs.push("projfile=" + trainingfolder + "/project.simple")
		commandargs.push("oritype=mic")
		commandargs.push("keys=intg")
		
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line of lines){
					var elements = line.split((/[ ]+/))
					elements.shift()
					if(elements.length > 0){
						selection.push(0)
					}
				}
				
				selection[arg['micrographid']] = 1
				
				var selectionfile = ""
				
				for (var state of selection){
					selectionfile += state + '\n'
				}
				
				fs.writeFileSync(trainingfolder + "/selection.txt", selectionfile)
				
				return
			})
			.then(() => {
				commandargs = []
				commandargs.push("prg=update_project_stateflags")
				commandargs.push("projfile=" + trainingfolder + "/project.simple" )
				commandargs.push("oritype=mic")
				commandargs.push("infile=" + trainingfolder + "/selection.txt")
				return spawn(command, commandargs, {capture: ['stdout']})
					.then(() => {
						return({status : "ok"})
					})
					.catch((error) => {
						return({status : "error"})
					})
			})
			.then(() => {
				command = "simple_distr_exec"
				commandargs = []
				commandargs.push("prg=pick")
				commandargs.push("refs=" + arg['refs'])
				commandargs.push("ndev=" + arg['ndev'])
				commandargs.push("thres=" + arg['thres'])
				commandargs.push("nparts=1")
				commandargs.push("mkdir=no")
				return spawn(command, commandargs, {capture: ['stdout'], cwd: trainingfolder})
					.then(() => {
						return({status : "ok"})
					})
					.catch((error) => {
						return({status : "error"})
					})
			})
			.then((result) => {
				console.log(result.stderr)
				console.log(result.stdout)
				var files = fs.readdirSync(trainingfolder)
				var boxfile
				var boxes 
				for (var file of files){
					if(file.includes(".box")){
						boxfile = trainingfolder + "/" + file
						boxes = fs.readFileSync(boxfile, {encoding : 'utf8'})
						break
					}
				}
				
				var lines = boxes.split("\n")
				var coordinates = []
				var boxsize
				for(var line of lines){
					var elements = line.split((/[ , \t]+/))
					if(elements.length > 2){
						coordinates.push([elements[0], elements[1]])
						boxsize = elements[3]
					}
				}
				return({coordinates : coordinates, boxsize : boxsize})
			})
			
	}
	
	public getMakePickRefsWidget(modules, arg){
		return new Promise((resolve, reject) => {
			resolve({html : this.makepickrefswidget({})})
		})
	}
	
	public createAutopickRefs(modules, arg){
		var folder
		return this.execute(modules, arg)
			.then((response) => {
				folder = response['folder']
				return({pickrefs : folder + "/pickrefs.mrc"})
			})
	}
	
	public readProjFile(projfile, section, keys){
		var spawn = require('child-process-promise').spawn
		var command = "simple_private_exec"
		var commandargs = []
		var queryresult = []
		commandargs.push("prg=print_project_vals")
		commandargs.push("projfile=" + projfile)
		commandargs.push("oritype=" + section)
		commandargs.push("keys=" + keys)
		return spawn(command, commandargs, {capture: ['stdout']})
			.then((result) => {
				var lines = result.stdout.split("\n")
				for(var line in lines){
					var elements = lines[line].split((/[ ]+/))
					elements.shift()
					queryresult.push(elements)
				}
				return(queryresult)
			})
	}
	
	public getBoxes(modules, arg){
		return new Promise((resolve, reject) => {
			var boxes = fs.readFileSync(arg['boxfile'], {encoding : 'utf8'})
			var lines = boxes.split("\n")
			var coordinates = []
			var boxsize
			for(var line of lines){
				var elements = line.split((/[ , \t]+/))
				if(elements.length > 2){
					coordinates.push([elements[1], elements[2]])
					boxsize = elements[3]
				}
			}
			resolve({coordinates : coordinates, boxsize : boxsize})
		})
	}
	
}

module.exports = new Module()
