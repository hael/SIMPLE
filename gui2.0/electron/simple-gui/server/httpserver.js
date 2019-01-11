const fs = require('fs')
const express = require('express')
const morgan = require('morgan')
const compression = require('compression')
const debug = require('debug')('inferno:HTTPserver')
const auth = require('http-auth')
const pug = require('pug')
//const densityserver = require("../../external/DensityServer/server/web-api")

class HTTPServer {

	constructor (args) {
		
		var multiuser = args['multiuser'] ? true : false
		
		var indexpug = pug.compileFile('server/views/index.pug')
		
		if(multiuser){
			this.basic = auth.basic({
				realm: "simple",
				file: __dirname + "/users.htpasswd"
			})
		
			this.basic.on('success', (result, req) => {
				console.log(`User authenticated: ${result.user}`);
			});
			 
			this.basic.on('fail', (result, req) => {
				console.log(`User authentication failed: ${result.user}`);
			});
			 
			this.basic.on('error', (error, req) => {
				console.log(`Authentication error: ${error.code + " - " + error.message}`);
			})
		}
		
		this.server = express()
		
		this.server.use(morgan('dev'))
		this.server.use(express.json())
		this.server.use(express.urlencoded({ extended: false }))
		this.server.use(express.static('public'))
		this.server.use(compression({ level: 6, memLevel: 9, chunkSize: 16 * 16384, filter: () => { return true; } }))
		if(multiuser){
			this.server.use(auth.connect(this.basic))
		}
		
		this.server.set('view engine', 'pug')
		
		this.server.get('/', (req, res) => {
			res.send(indexpug())
		})
		
		this.server.get('/image', (req, res) => {
		//	modules['available']['core']['getImage'](req.query).then((result) => {
		//		res.contentType('jpg')
		//		res.send(result.image)
		//	})
		})
		
		this.server.post('/', (req, res) => {
			var fnc = req.body.fnc
			var mod = req.body.mod
			var arg = req.body.arg
			if(multiuser){
				arg['user'] = req.user
			} else {
				arg['user'] = "simple"
			}
			if(fnc != "undefined" && mod != "undefined"){
				//modules.available[mod][fnc](modules, arg).then(function(result){
				//	if(result.view){
				//		res.render(result.view, result)
				//	}else{
				//		res.send(result)
				//	}
				//}) 
			}
		})
		

		this.server.use(this.errorHandler)

	//	this.densityserver.default(this.server)
	}
	
	start(port) {
		this.server.listen(port)
		console.log("HTTP server is running on port " + port);
	}
	
	stop(){
		this.server.close()
	}
	
	errorHandler(err, req, res, next) {
		debug(err)
		if (res.headersSent) {
			return next(err)
		}
		res.status(500)
		res.render('error', { error: err })
	}

}

module.exports = HTTPServer
