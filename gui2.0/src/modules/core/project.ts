declare var modules

import * as fs from 'fs'

export default class Project{

	private projectselector
	private projecthistory
	private logview
	private createproject
	private sqlite3
	
	constructor(){
		const pug = require('pug')
		this.projectselector = pug.compileFile('views/core-projectselector.pug')
		this.projecthistory = pug.compileFile('views/core-projecthistory.pug')
		this.logview = pug.compileFile('views/core-logview.pug')
		this.createproject = pug.compileFile('views/core-createproject.pug')
		this.sqlite3 = require('sqlite3').verbose()
		var query = "CREATE TABLE IF NOT EXISTS projects (id INTEGER PRIMARY KEY AUTOINCREMENT,	name TEXT, user TEXT, description TEXT, history TEXT UNIQUE, folder TEXT)"
		this.sqlQuery(query)
	}
	
	get(user){
		var query = "SELECT * FROM projects WHERE user='" + user + "'"
		return this.sqlQuery(query)
			.then((rows) => {
				var projects = {}
				projects['projects'] = rows
				return ({html : this.projectselector(projects)})
			})
	}
	
	getHistory(arg){
		var query = "SELECT * FROM " + arg['history']
		return this.sqlQuery(query)
			.then((rows) => {
				var history = {}
				history['history'] = rows
				return ({html : this.projecthistory(history)})
			})
	}
	
	getNew(){
		return new Promise((resolve, reject) => {
			resolve({html : this.createproject({})})
		})
	}
	
	createNew(arg){
		var projecttable = "t" + Math.floor(Math.random() * Math.floor(1000000))
		var query = "CREATE TABLE IF NOT EXISTS " + projecttable + " (id INTEGER PRIMARY KEY AUTOINCREMENT,	name TEXT, description TEXT, arguments TEXT, folder TEXT, pid TEXT, status TEXT, type TEXT, view TEXT)"
		for(var module in modules['available']){
			if(typeof(modules['available'][module]['setupProject']) == "function"){
				modules['available'][module].setupProject(arg)
			}
		}
		return this.sqlQuery(query)
			.then((rows) => {
				var query = "INSERT INTO projects (name, user, description, history, folder) VALUES ('" + arg['keys']['keyname'] + "','" + arg['user'] + "','" + arg['keys']['keydescription'] + "','" + projecttable + "','" + arg['keys']['keyfolder'] + "')"
				return this.sqlQuery(query)
			})
			.then((rows) => {
				return ({})
			})
	}
	
	getLog(arg){
		return new Promise((resolve, reject) => {
			var log
			if(fs.existsSync(arg['folder'] + "/task.log")) {
				log = fs.readFileSync(arg['folder'] + "/task.log", 'utf8')
			}
			resolve(log)
		})
		.then((data) => {
			return ({html : this.logview({data : data})})
		})
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
