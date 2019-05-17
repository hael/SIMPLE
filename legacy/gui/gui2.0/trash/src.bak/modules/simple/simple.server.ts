import * as fs from 'fs'

class Module {
	
	public metadata = {
			"moduletitle" :"Simple",
			"modulename" :"simple",
			"tasks" : {
			}
	}
		
	constructor(){
		console.log("Loaded module Simple")
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
				task['pages'].push(page)
			}
			this.metadata['tasks'][commandkeys[0]] = task
		}
	}
	
	public execute = function(modules, arg){
		var spawn = require('child_process').spawn
		var command = arg['executable']
		var type = arg['type']
		var keys = arg['keys']
		var keynames = Object.keys(keys);
		var commandargs = []
		commandargs.push("prg=" + type)
		for(var key of keynames){
			if(keys[key] != ""){
				if(keys[key] == "true"){
					commandargs.push(key.replace('key', '') + "=yes")
				}else if(keys[key] == "false"){
					commandargs.push(key.replace('key', '') + "=no")
				}else{
					commandargs.push(key.replace('key', '') + "=" + keys[key])
				}
			}
		}
		
		return modules['available']['core']['createTask'](modules, arg)
			.then(function(json){
				console.log(command)
				console.log(commandargs)
				var executeprocess = spawn(command, commandargs, {detached: true, cwd: json['jobfolder']})
				executeprocess.on('exit', function(code){
					console.log(`child process exited with code ${code}`);
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
				return({})
			})
	}
}

module.exports = new Module()
