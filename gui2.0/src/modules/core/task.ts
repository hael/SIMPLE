declare var modules

import * as fs from 'fs'

export default class Task{

	private taskselector
	private taskinput
	private sqlite3
	
	constructor() {
		const pug = require('pug')
		this.taskselector = pug.compileFile('views/core-taskselector.pug')
		this.taskinput = pug.compileFile('views/core-taskinput.pug')
		this.sqlite3 = require('sqlite3').verbose()
	}
	
	getSelector(arg) {
		var tasks = {}
		tasks['folder'] = arg['folder']
		tasks['tasks'] = []
		for(var module in modules['available']){
			if(typeof(modules['available'][module]['metadata']) == "object"){
				console.log(module)
				tasks['tasks'].push(modules['available'][module]['metadata'])
			}
		}
		return new Promise((resolve, reject) => {
			resolve({html : this.taskselector(tasks)})
		})
	}
	
	refineSelector(arg){
		return new Promise( async (resolve, reject) => {
			var refinedSelection = {}
				if(fs.existsSync(arg['path'])){
				for(var module in modules['available']){
					if(typeof(modules['available'][module]['refineSelection']) == "function"){
						await modules['available'][module]['refineSelection'](arg['path'])
							.then((val) => {refinedSelection[modules['available'][module]['metadata']['moduletitle']] = val})
					}
				}
			}
			resolve({refinedselection : refinedSelection})
		})
	}
	
	setup(arg){
		var module = arg.module
		var task = arg.task
		var inputpath = arg['inputpath']
		return new Promise((resolve, reject) => {
			var taskarray = {}
			taskarray['taskarray'] = modules.available[module].metadata['tasks'][task]
			taskarray['taskarray']['modulename'] = modules.available[module].metadata['modulename']
			taskarray['taskarray']['moduletitle'] = modules.available[module].metadata['moduletitle']
			taskarray['taskarray']['inputpath'] = inputpath
			resolve({html : this.taskinput(taskarray)})
		})
	}
	
	createNew(arg){
		var query = "INSERT into " + arg['projecttable'] + " (name, description, arguments, status, view, type) VALUES ('" + arg['name'] + "','" + arg['description'] + "','" + JSON.stringify(arg) + "','running', '" + JSON.stringify(arg['view']) + "', '" + arg['type'] + "')"

		var jobfolder
		var jobid 
		return this.sqlQuery(query)
			.then((rows) => {
				var query = "SELECT seq FROM sqlite_sequence WHERE name='" + arg['projecttable'] + "'"
				return this.sqlQuery(query)
			})
			.then((rows) => {
				jobid = rows[0]['seq']
				jobfolder = arg['projectfolder'] + "/" + jobid + "_" + arg['type']
				fs.mkdirSync(jobfolder)
				query = "UPDATE " + arg['projecttable'] + " SET folder='" + jobfolder + "' WHERE id=" + jobid 
				return this.sqlQuery(query)
			})
			.then((rows) => {
				return({jobfolder: jobfolder, jobid : jobid})
			})
	}
	
	updatePid(table, jobid, jobpid){
		var query = "UPDATE " + table + " SET pid='" + jobpid + "' WHERE id=" + jobid
		return this.sqlQuery(query)
	}
	
	updateStatus(table, jobid, status){
		var query = "UPDATE " + table + " SET status='" + status + "' WHERE id=" + jobid
		return this.sqlQuery(query)
	}
	
	delete(arg){
		var query = "DELETE FROM " + arg['history'] + " WHERE id=" + arg['jobid']
		return this.sqlQuery(query)
	}
	
	kill(arg){
		process.kill(arg['jobpid'])
		return this.updateStatus(arg['history'], arg['jobid'], "Killed")
	}
	
	private sqlQuery(query){
		return new Promise((resolve,reject) => {
			var db = new this.sqlite3.Database('simple.sqlite')
			db.all(query, (err, rows) => {
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
