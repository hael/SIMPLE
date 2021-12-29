const fs = require(global.base + '/node_modules/fs-extra');
const pug = require(global.base + '/node_modules/pug')
const simpleexec = require('./simpleExec')
const path = require('path')
//const spawn = require(global.base + '/node_modules/child-process-promise').spawn

class Relion{
	
	constructor(){
		this.jobfolders = []
	}
		
	analyse(dir){
		this.jobfolders = []
		if (fs.existsSync(dir)){
			try {
				this.findJobs(dir)
			} catch (err){
				return({jobsfound:false, dir:dir})
			}
			
			return(this.getJobData())
		} else {
			return({exists:false})
		}
	}
	
	findJobs(dir){
		for (var content1 of fs.readdirSync(dir)){
			if (fs.lstatSync(dir + "/" + content1).isDirectory()){
				for (var content2 of fs.readdirSync(dir + "/" + content1)){
					if (content2.includes("job")){
						if (fs.existsSync(dir + "/" + content1 + "/" + content2 + "/" + "job_pipeline.star")){
							this.jobfolders.push(dir + "/" + content1 + "/" + content2)
						}
					}
				}
			}
		}
	}
	
	getJobData(){
		var alljobs = []
		for (var dir of this.jobfolders){
			var object = {}
			
			object.id = this.getJobId(dir)
			object.type = this.getJobType(dir)
			object.directory = dir
			object.parent = this.getJobParent(dir)
			object.params = this.getJobParams(dir)
			
			alljobs.push(object)
		}
		
		return alljobs
	}
	
	getJobId(dir){
		var dirsplit = dir.split("/")
		var id = dirsplit[dirsplit.length - 1]/*.replace("job","")*/
		return id
	}
	
	getJobType(dir){
		var dirsplit = dir.split("/")
		var type = dirsplit[dirsplit.length - 2]
		return type
	}
	
	getJobParent(dir){
		// job_pipeline.star
		var pipeline = fs.readFileSync(dir + "/job_pipeline.star","UTF8")
		
		var pipesplit = pipeline.split("\n")
		
		for (var i = 0 ; i < pipesplit.length ; i++){
			if (pipesplit[i].includes("_rlnPipeLineEdgeProcess #2")){
				return pipesplit[i+1].split("/")[1]
				break
			}
		}
		
		return null
		
	}
	
	getJobParams(dir){
		// run.job
		var params = {}
		var paramstring = fs.readFileSync(dir + "/run.job","UTF8")
		
		paramstring = paramstring.substring(0, paramstring.length - 1)	// cut newline off end
		
		for (var param of paramstring.split("\n")){
			var paramsplit = param.split(" == ")
			params[paramsplit[0]] = paramsplit[1]
		}

		return params
	}
  
}

module.exports = Relion
