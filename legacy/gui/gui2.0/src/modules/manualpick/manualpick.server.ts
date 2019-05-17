declare var modules
import * as fs from 'fs'

class Module {

	private manualpickwidget

	public metadata = {
			"moduletitle" :"Manual Picking",
			"modulename" :"manualpick",
			"tasks" : {
				manualpick : { 
					name: 'manualpick',
					descr_short: 'Manual Pick',
					descr_long: 'is a program for centering a volume and mapping the shift parameters back to the particle images',
					executable: 'simple_exec',
					pages: [ 
						 { 
							title: 'picking', 
							keys: [
								{
									key: "pick",
									keytype: "widget",
									widget: "manualpickingwidget.view(this)",
									descr_short: "Pick",
									descr_long: "Pick",
									descr_placeholder: "json string of coordinates",
									required: true
								}
							]
						}
					] 
				}
			}
	}
		
	constructor(){
		process.stdout.write("Loading module - ManualPick ... ")
		const pug = require('pug')
		this.manualpickwidget = pug.compileFile('views/manualpick-manualpickingwidget.pug')
		process.stdout.write("Done\n")
	}
	
	public getManualPickWidget(modules, args){
		var path = require('path')
		var micrographs = []
		if(args['micfile'].includes('.simple')){
			return modules['available']['simple'].readProjFile(args['micfile'], "mic", "intg,xdim,ydim")
				.then((mics) => {
					for (var micrograph of mics){
						if(micrograph.length != 0){
							var micname = path.basename(micrograph[2])
							micrographs.push({name: micname, path: micrograph[2], xdim: micrograph[3], ydim: micrograph[4]})
							
						}
					}

					return({html : this.manualpickwidget({micrographs: micrographs})})
				})
		}
	}
	
	public execute(modules, args){
		return modules['available']['core']['taskCreate'](modules, args)
			.then((json) => {
				var coords = JSON.parse(args['keys']['keypick'])
				fs.appendFileSync(json['jobfolder'] + "/task.log", "Writing coordinates ...")
				for(var coord of coords){
					fs.appendFileSync(json['jobfolder'] + "/" + coord['name'].replace(".mrc", ".box"), coord['xcoord'] + "\t" + coord['ycoord'] + "\t" + coord['boxsize'] + "\t" + coord['boxsize'] + "\n")
				}
				fs.appendFileSync(json['jobfolder'] + "/task.log", "Wrote " + coords.length + " coordinates")
				modules['available']['core']['updateStatus'](args['projecttable'], json['jobid'], "Finished")
				return({})
			})
	}
	
}

module.exports = new Module()
