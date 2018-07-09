import fetch from 'node-fetch'

class Module {
	
	constructor(){
		console.log("Loaded module Arctica Sync")
	}
	
	metadata2 = {
		"name" : "Arctica Sync",
  		"modulename" : "arctica_sync",
  		"tasks" : [
  			{
  				"name": "Arctica Sync",
  				"descr_short": "Sync micrographs from the arctica",
 				"descr_long": "Sync micrographs from the arctica",
 				"executable": "newSync",
 				"pages" : [
  					{
  						"title": "",
  						"keys" : [
  							{
 	 							"key":"",
 								"keytype": "",
 								"descr_short": "",
 								"descr_long":"",
 								"descr_placeholder": "",
 								"required": ""
 							}
  						]
  					}
  				]
  
  			}
  		]
	}

	newSync(modules, arg){
		arg['fnc'] = "new"
		return fetch("http://localhost:8089", {
			method: 'POST',
			body : JSON.stringify(arg),
			headers: {'content-type': 'application/json'}
		}).then(function(){
			return({viw : "taskinput", taskarray : []})
		}).catch(function(err){
			console.log(err)
		})
	}
	
	status(modules, arg){
		return fetch("http://localhost:8089", {
			method: 'POST',
			body : JSON.stringify({fnc : "status"}),
			headers: {'content-type': 'application/json'}
		}).then(function(result){
			return(result.json())
		}).catch(function(err){
			console.log(err)
			return({html : "Cannot communicate with sync worker. Is it running?"})
		})
	}
	
}



/*	metadata : function() {
		
	}
		header : "Arctica Sync",
		name : "arctica_sync",
		tasks : {
			Sync : {
				help : "Sync micrographs from the arctica",
				arguments : [ 
					{
						name	: "src",
						label	: "Source",
						type	: "folder",
						help	: "Source folder"
					},
					{
						name	: "dst",
						label	: "Destination",
						type	: "folder",
						help	: "Destination folder"
					},
					{
						name	: "idt",
						label	: "Identifier",
						type	: "text",
						help	: "Filename identifier"
					},
					{
						name	: "rem",
						label	: "Remove",
						type	: "check",
						help	: "Remove"
					}
				],
				widget : {
					execute : "status"
					
				},
				execute : "newSync"
			}
		}
	},
*/	


module.exports = new Module()
