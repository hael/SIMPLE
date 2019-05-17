
export default class HeaderButtons{

	private headerbuttons
	
	constructor(){
		const pug = require('pug')
		this.headerbuttons = pug.compileFile('views/core-headerbuttons.pug')
	}

	get(){
		var buttons = {}
		buttons['buttons'] = []
		buttons['buttons'].push({value : "View", function :"browser.show({buttons : [{ name : 'view2d', action : 'viewer.view2d()'}, { name : 'view3d', action : 'viewer.view3d()'}], path : '/tmp', gauze : true})"})
		buttons['buttons'].push({value : "Task", function :"taskselector.show()"})
		return new Promise((resolve, reject) => {
			resolve(
				{html : this.headerbuttons(buttons)}
			)
		})
		
	}
	
}
