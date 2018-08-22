declare var modules

export default class View{

	private view2dview
	private view3dview
	private pack
	
	constructor() {
		const pug = require('pug')
		this.view2dview = pug.compileFile('views/core-view2d.pug')
		this.view3dview = pug.compileFile('views/core-view3d.pug')
		this.pack = require("../../../external/DensityServer/pack/main.js")
	}
	
	public view2d (arg) {
		var view = {}
		if(arg['file'].includes(".mrc")){
			return modules['available']['core']['readMRCHeader'](arg['file'])
				.then((header) => {
					view['thumbcount'] = header['nz']
					view['path'] = arg['file']
					return({html : this.view2dview(view)})
				})
		}
	}
	
	public view3d(arg) {
		var id = Math.floor(Math.random() * 10000)
		var config = {
				input: [ { name: 'em', filename: arg['file']}],
				isPeriodic: false,
				outputFilename: "/tmp/" + id + ".mdb",
				blockSize: 96
			}
			
		return(this.pack.default(config.input, config.blockSize, config.isPeriodic, config.outputFilename))
		.then(() => {
				return({html : this.view3dview(), mdb : id})
		})
	}
}
