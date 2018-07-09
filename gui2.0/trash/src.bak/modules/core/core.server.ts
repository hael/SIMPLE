import * as fs from 'fs'
import * as binaryfile from 'binary-file'
import * as sharp from 'sharp'

class Module {
	
	constructor(){
		console.log("Loaded module Core")
	}
	
	folderListing (modules, arg) {
		return new Promise(function(resolve, reject){
			var path = arg.pth
			var view = arg.viw
			var filter = arg.flt
			var contents = fs.readdirSync(path)
			var files = []
			var folders = []
			for (var object of contents){
				if(object.charAt(0) != "."){
					var stat = fs.lstatSync(path + "/" + object)
					if(stat.isFile()){
						if(filter != undefined){
							if(object.includes(filter)){
								files.push(object)
							}
						}else{
							files.push(object)
						}
					} else if (stat.isDirectory()){
						folders.push(object)
					}
				}
			}
			files.sort()
			folders.sort()
			resolve({ "path" : path, "files" : files, "folders" : folders, "view" : view, filter : filter})
		})
	}

	taskListing(modules, arg) {
		return new Promise(function(resolve, reject){
			var tasks = []
			for (var module in modules.available){
				var metadata = modules.available[module].metadata
				if(metadata != undefined){
					tasks.push(metadata)
				}
			} 
			resolve({view : "taskselector", tasks : tasks})
		})
	}

	taskInput(modules, arg){
		var module = arg.module
		var task = arg.task
		return new Promise(function(resolve, reject){
			var taskarray = modules.available[module].metadata['tasks'][task]
			taskarray['modulename'] = modules.available[module].metadata['modulename']
			taskarray['moduletitle'] = modules.available[module].metadata['moduletitle']
			console.log(taskarray)
			resolve({view : "taskinput", taskarray : taskarray})
		})
	}

	createMDB(modules, arg) {
		return new Promise(function(resolve, reject){
			var id = Math.floor(Math.random() * 10000);
			var packmain = require("../../../external/DensityServer/build/pack/main");
			var config = {
				input: [ { name: 'em', filename: arg.pth }],
				isPeriodic: false,
				outputFilename: "/tmp/" + id + ".mdb",
				blockSize: 96
			}
			packmain.default(config.input, config.blockSize, config.isPeriodic, config.outputFilename)
			
			resolve({ "pth" : id, "fil" : arg.pth })
		})
	}

	view2DInit(modules, arg){
		var nz
		var path = arg.pth
		const datafile = new binaryfile(path, 'r', true)
		return datafile.open()
			.then(function(){
				console.log("Opened " + path)
				return datafile.readInt32(8)
			})
			.then(function(val){
				nz = val
				return datafile.close()
			})
			.then(function () {
				console.log(nz)
				console.log('File closed')
				var thumbnails = []
				for (var i = 0; i < nz; i++){
					var thumbnail = [{file : path, frame : i}]
					thumbnails.push(thumbnail)
				}
				return({thumbnails : thumbnails, view : "thumbnailview" })
			})
	}

	mrc2JPEG(modules, path, frame, width){
		var nx
		var ny
		var mode
		var nsymbt
		var data
		
		const datafile = new binaryfile(path, 'r', true)
		return datafile.open()
			.then(function(){
				return datafile.readInt32(0)
			})
			.then(function(val){
				nx = val
				return datafile.readInt32(4)
			})
			.then(function(val){
				ny = val
				return datafile.readInt32(12)
			})
			.then(function(val){
				mode = val
				return datafile.readInt32(92)
			})
			.then(function(val){
				nsymbt = val
				if(mode == 2){
					return datafile.read(nx * ny * 4, 1024 + nsymbt + Number(frame) * nx * ny * 4)
				}// need to resturn on error!
			})
			.then(function(val){
				data = val
				return datafile.close()
			})
			.then(function(val){
				var mean = 0
				var variance = 0
				var tmp
				
				for(var i = 0; i < nx * ny; i++){
					tmp = mean
					mean = tmp + data.readFloatLE(i * 4)
				}
				mean /= nx * ny
				
				for(var i = 0; i < nx * ny; i++){
					tmp = variance
					variance = tmp + Math.pow((data.readFloatLE(i * 4) - mean), 2)
				}
				variance /= nx * ny
				var sd = Math.sqrt(variance)
				var pixels = []
				var norm
				
				for(var i = 0; i < nx * ny; i++){
					norm = (data.readFloatLE(i * 4) - mean) / sd
					if(norm > 10){
						pixels.push(255)
					} else if (norm < -10){
						pixels.push(0)
					} else {
						pixels.push(Math.round((norm * 12.8) + 128))
					}
				}
				
				var pixbuf =  Buffer.from(pixels)
	 
				return sharp(pixbuf, { raw : {
					width : nx,
					height : ny,
					channels : 1
					}
				}).resize(Number(width))
				  .jpeg()
				  .toBuffer()
			})
			.then(function (image) {
				return({image : image})
			})
	}
	
	save2D(modules, arg) {
		var source = arg.src
		var destination = arg.pth
		var selection = arg.selection
		var nx
		var ny
		var mode
		var nsymbt
		var data
		
		const datain = new binaryfile(source, 'r', true)
		const dataout = new binaryfile(destination, 'w', true)
		
		return datain.open()
			.then(function(){
				return datain.readInt32(0)
			})
			.then(function(val){
				nx = val
				return datain.readInt32(4)
			})
			.then(function(val){
				ny = val
				return datain.readInt32(12)
			})
			.then(function(val){
				mode = val
				return datain.readInt32(92)
			})
			.then(function(val){
				nsymbt = val
				return dataout.open()
			})
			.then(function(){
				return datain.read(1024, 0)
			})
			.then(function(val){
				return dataout.write(val)
			})
			.then(async function(val){
				for(var selected of selection){ 
					await datain.read(nx * ny * 4, 1024 + nsymbt + Number(selected) * nx * ny * 4).then(function(val){return dataout.write(val)})
				}
				return datain.close()
			})
			.then(function(){
				return dataout.writeInt32(selection.length, 8)
			})
			.then(function(){
				return dataout.writeInt32(0, 92)
			})
			.then(function(){
				return dataout.close()
			})
			.then(function(){
				return ({joe :"joe"})
			}).catch(function (err) {
				console.log(`There was an error: ${err}`);
			})
	}
	
	createProject(modules, arg){
		var projecttable = "t" + Math.floor(Math.random() * Math.floor(1000000))
		
		var query = "CREATE TABLE IF NOT EXISTS projects (id INTEGER PRIMARY KEY AUTOINCREMENT,	name TEXT, user TEXT, description TEXT, history TEXT UNIQUE, folder TEXT)"
		this.sqlQuery(query)
		
		var query = "CREATE TABLE IF NOT EXISTS " + projecttable + " (id INTEGER PRIMARY KEY AUTOINCREMENT,	name TEXT, description TEXT, arguments TEXT, folder TEXT, pid TEXT, status TEXT, type TEXT, view TEXT)"
		this.sqlQuery(query)
		
		var query = "INSERT INTO projects (name, user, description, history, folder) VALUES ('" + arg['name'] + "','" + arg['user'] + "','" + arg['description'] + "','" + projecttable + "','" + arg['folder'] + "')"
		this.sqlQuery(query)
		
		return new Promise(function(resolve, reject){
			resolve({"status" : "success"})
		})
	}
	
	listProjects(modules, arg){
		var query = "SELECT * FROM projects"
		return this.sqlQuery(query)
			.then(function(rows){
				return ({view :"projectslist", projects: rows})
			})
	}
	
	projectHistory(modules, arg){
		var query = "SELECT * FROM projects WHERE id=" + arg['projectid']
		return this.sqlQuery(query)
			.then(function(rows){
				query = "SELECT * FROM " + rows[0]['history']
				return(modules['available']['core'].sqlQuery(query))
			})
			.then(function(rows){
				return ({view :"projecthistory", history: rows})
			})
			
	}
	
	createTask(modules, arg){
		console.log("rtask")
		var query = "INSERT into " + arg['projecttable'] + " (name, description, arguments, status) VALUES ('joe','joe', '"+ arg + "','running')"
		var jobfolder
		return this.sqlQuery(query)
			.then(function(rows){
				query = "SELECT seq FROM sqlite_sequence WHERE name='" + arg['projecttable'] + "'"
				return modules['available']['core'].sqlQuery(query)
			})
			.then(function(rows){
				var jobid = rows[0]['seq']
				jobfolder = arg['projectfolder'] + "/" + jobid + "_" + arg['type']
				fs.mkdirSync(jobfolder)
				query = "UPDATE " + arg['projecttable'] + " SET folder='" + jobfolder + "' WHERE id=" + jobid 
				return modules['available']['core'].sqlQuery(query)
			})
			.then(function(rows){
				return({jobfolder: jobfolder})
			})
	}	

	private sqlQuery(query){
		console.log(query)
		return new Promise(function(resolve,reject){
			var sqlite3 = require('sqlite3').verbose()
			var db = new sqlite3.Database('simple.sqlite')
			db.all(query, function(err, rows){
				if(err == null){
					resolve(rows)
				} else {
					reject(err)
				}
			})
			db.close()
		})
		
	}
	

}

module.exports = new Module()


