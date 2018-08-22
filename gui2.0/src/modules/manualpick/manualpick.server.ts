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
					//testing
					 micrographs = [
						{name: "test", path: "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337193_Data_23336346_23336347_20180523_1726-163966_intg.mrc", xdim: 3860, ydim: 3860},
						{name: "test", path: "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337193_Data_23336361_23336362_20180523_1725-163963_intg.mrc", xdim: 3860, ydim: 3860},
						{name: "test", path:  "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337193_Data_23336368_23336369_20180523_1726-163964_intg.mrc", xdim: 3860, ydim: 3860},
						{name:  "test", path: "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337193_Data_23336374_23336375_20180523_1726-163965_intg.mrc", xdim: 3860, ydim: 3860},
						{name: "test", path:  "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337193_Data_23336381_23336382_20180523_1727-163967_intg.mrc", xdim: 3860, ydim: 3860},
						{name: "test", path:  "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337194_Data_23336346_23336347_20180523_1730-163971_intg.mrc", xdim: 3860, ydim: 3860},
						{name:  "test", path: "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337194_Data_23336361_23336362_20180523_1728-163968_intg.mrc", xdim: 3860, ydim: 3860},
						{name: "test", path:  "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337194_Data_23336368_23336369_20180523_1729-163969_intg.mrc", xdim: 3860, ydim: 3860},
						{name: "test", path:  "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337194_Data_23336374_23336375_20180523_1729-163970_intg.mrc", xdim: 3860, ydim: 3860},
						{name: "test", path:  "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337194_Data_23336381_23336382_20180523_1730-163972_intg.mrc", xdim: 3860, ydim: 3860},
						{name: "test", path:  "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337195_Data_23336346_23336347_20180523_1733-163976_intg.mrc", xdim: 3860, ydim: 3860},
						{name:  "test", path: "/beegfs/lea/ApoFer/2018-05-23_KriosK2/1_pipeline/pipeline/micrographs/FoilHole_23337195_Data_23336361_23336362_20180523_1732-163973_intg.mrc", xdim: 3860, ydim: 3860}
					 ]
					
					return({html : this.manualpickwidget({micrographs: micrographs})})
				})
		}
	}
	
}

module.exports = new Module()
