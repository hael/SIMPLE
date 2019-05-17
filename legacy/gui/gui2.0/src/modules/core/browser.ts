import * as fs from 'fs'



export default class Browser{

	private browser
	
	constructor() {
		const pug = require('pug')
		this.browser = pug.compileFile('views/core-browser.pug')
	}
	
	getAll(arg) {
		return new Promise((resolve, reject) => {
			var files = []
			var folders = []
			var contents = fs.readdirSync(arg['path'])
			for (var object of contents){
				if(object.charAt(0) != "."){
				//	var stat = fs.lstatSync(arg['path'] + "/" + object)
					var stat = fs.statSync(arg['path'] + "/" + object)
					if(stat.isFile()){
						if(arg['filter'] != undefined){
							if(object.includes(arg['filter'])){
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
			var view = {}
			view['files'] = files
			view['folders'] = folders
			view['path'] = arg['path']
			view['filter'] = arg['filter']
			view['buttons'] = arg['buttons']
			resolve({html : this.browser(view)})
		})
	}
	
	
}
