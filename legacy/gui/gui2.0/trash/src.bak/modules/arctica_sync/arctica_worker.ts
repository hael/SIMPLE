var port = 8089

var fs = require('fs');

var currentsyncs = []

var express = require('express')

var server = express()

server.use(express.json())
server.use(express.urlencoded({ extended: false }))

server.post('/', function (req, res) {
	console.log(req.body)
	var fnc = req.body.fnc
	if(fnc == "status"){
		var status = getStatus()
		res.send({html : status})
	} else if (fnc == "new"){
		var sync = {}
		sync['source'] = req.body.src
		sync['destination'] = req.body.dst
		
		sync['remove'] = req.body.rem
		sync['identifier'] = req.body.idt
		sync['lastnew'] = 0
		sync['count'] = 0
		sync['queue'] = []
		var datasetname =  req.body.src.split('\\').pop().split('/').pop()
		if (!fs.existsSync(req.body.dst + "/" + datasetname)){
			fs.mkdirSync(req.body.dst + "/" + datasetname)
		}
		if (!fs.existsSync(req.body.dst + "/" + datasetname + "/movies")){
			fs.mkdirSync(req.body.dst + "/" + datasetname + "/movies")
		}
		sync['realdestination'] = req.body.dst + "/" + datasetname + "/movies"
		currentsyncs.push(sync)
		res.send({})
	}
})

function getStatus(){
	var status = "<table>"
	for(var sync of currentsyncs){
		if(sync['source'] != "undefined"){
			status += "<tr>"
			status += "<td>" + sync['source'] + "</td>"
			status += "<td>" + sync['destination'] + "</td>"
			status += "<td>" + sync['identifier'] + "</td>"
			status += "<td>" + sync['remove'] + "</td>"
			status += "<td>" + sync['lastnew'] + "</td>"
			status += "<td>" + sync['count'] + "</td>"
			status += "</tr>"
		}
	}
	status += "</table>"
	return(status)
}

function findFiles(){
	console.log("finding files")
	var find = require('find')
	for(var sync of currentsyncs){
		find.eachfile(/\sync['identifier']/, sync['source'], function(file) {
			sync['queue'].push(file)
		}).error(function(err) {})
	}
}

function syncFiles(){
	console.log("syncing files")
	for(var sync of currentsyncs){
		console.log(sync)
		var date
		while(sync['queue'].length > 0){
			if(sync['source'] != "undefined"){
				var filename =  sync['queue'][0].split('\\').pop().split('/').pop()
				console.log("copy " + sync['queue'][0] + " to " + sync['realdestination'] + "/" + filename)
				//fs.createReadStream(sync['queue'][0]).pipe(fs.createWriteStream(sync['realdestination'] + "/" + filename))
				// remove file if ok
				sync['queue'].shift()
				date = new Date()
				sync['lastnew'] = date.getTime()
			}
		}
	}
	date = new Date()
	var now = date.getTime()
	for(var sync of currentsyncs){
		if((now - sync['lastnew']) > 1440000){
			sync = {}
		}
	}
	syncloop(3000)
}

function syncloop(timeout){
	setTimeout(function(){
		findFiles()
	}, timeout)
	setTimeout(function(){
		syncFiles()
	}, timeout + 6000)
}

server.listen(port)
syncloop(3000)
