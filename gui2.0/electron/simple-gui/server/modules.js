const fs = require('fs')

class Modules {
	
	constructor(){
		this.available = {}
		var modulefolders = fs.readdirSync("modules")
		for (var modulename of modulefolders){
			if(fs.existsSync("modules" + "/" + modulename + "/" + modulename + ".server.js")){
				var stat = fs.statSync("modules" + "/" + modulename + "/" + modulename + ".server.js")
				if(stat.isFile()){
					this.available[modulename] = require("../modules" + "/" + modulename + "/" + modulename + ".server.js")
				}
			}
		}
	}
	
}

module.exports = Modules
