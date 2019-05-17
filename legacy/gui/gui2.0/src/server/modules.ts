import * as fs from 'fs'

export default class Modules {
	
	public available
		
	constructor(){
		this.available = {}
		var modulefolders = fs.readdirSync("modules")
		for (var modulename of modulefolders){
			this.available[modulename] = {}
			if(fs.existsSync("modules" + "/" + modulename + "/" + modulename + ".server.js")){
				var stat = fs.lstatSync("modules" + "/" + modulename + "/" + modulename + ".server.js")
				if(stat.isFile()){
					this.available[modulename] = require("../modules" + "/" + modulename + "/" + modulename + ".server.js")
				}
			}
		}
	}
	
}
