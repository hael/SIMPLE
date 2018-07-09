import * as fs from 'fs'

/*	module : {
 * 		name : <>,
 * 		available : <>,
 * 		tasks : [
 * 			{
 * 				"name": <>,
 * 				"descr_short": <>,
 *				"descr_long": <>,
 *				"executable": <>,
 *				"pages" : [
 * 					{
 * 						"title": <>,
 * 						"keys" : [
 * 							{
 *	 							"key": <>,
 *								"keytype": <file|folder|check|text>,
 *								"descr_short": <>,
 *								"descr_long": <>,
 *								"descr_placeholder": <>,
 *								"required": <true/false>
 *							}
 * 						]
 * 					}
 * 				]
 * 
 * 			}
 * 		]
 * }
*/



class Modules {
	
	public available
		
	constructor(){
		this.available = {}
		var modulefolders = fs.readdirSync("build/modules")
		for (var modulename of modulefolders){
			this.available[modulename] = {}
			if(fs.existsSync("build/modules" + "/" + modulename + "/" + modulename + ".server.js")){
				var stat = fs.lstatSync("build/modules" + "/" + modulename + "/" + modulename + ".server.js")
				if(stat.isFile()){
					this.available[modulename] = require("../modules" + "/" + modulename + "/" + modulename + ".server.js")
				}
			}
		}
	}
}

class HTTPServer {

	private		server
	private		express
	private		logger
	private		debug
	private		compression
	private		densityserver
	private		auth
	private		basic
	
	constructor () {
		/* Import dependencies */
		this.express = require('express')
		this.logger = require('morgan')
		this.compression = require('compression')
		this.debug = require('debug')('inferno:HTTPserver')
		this.densityserver = require("../../external/DensityServer/build/server/web-api")
		this.auth = require('http-auth')
		this.basic = this.auth.basic({
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
		
		/* Create server */
		this.server = this.express()
		
		/* Setup server */
		this.server.use(this.logger('dev'))
		this.server.use(this.express.json())
		this.server.use(this.express.urlencoded({ extended: false }))
		this.server.use(this.express.static('public'))
		this.server.use(this.compression({ level: 6, memLevel: 9, chunkSize: 16 * 16384, filter: function () { return true; } }))
		this.server.use(this.auth.connect(this.basic))
		
		/* Setup PUG view engine */
		this.server.set('view engine', 'pug')
		
		/* Setup routes */
		this.server.get('/', function(req, res) {
			res.render('index', {})
		})
		
		this.server.get('/image', function(req, res) {
			var fnc = req.query.fnc
			var mod = req.query.mod
			var pth = req.query.pth
			var frm = req.query.frm
			var wth = req.query.wth
			
			if(fnc != "undefined" && mod != "undefined"){
			
				modules.available[mod][fnc](modules, pth, frm, wth).then(function(result){
					res.contentType('jpeg')
					res.send(result.image)
				}) // add catch send error
			}
		})
		
		this.server.get('/landing', function(req, res) {
			res.render('landing', {ui : ui})
		})
		
		this.server.post('/', function (req, res) {
			var fnc = req.body.fnc
			var mod = req.body.mod
			var arg = req.body.arg
			arg['user'] = req.user
			if(fnc != "undefined" && mod != "undefined"){
				modules.available[mod][fnc](modules, arg).then(function(result){
					if(result.view){
						res.render(result.view, result)
					}else{
						res.send(result)
					}
				}) // add catch send error
			}
		})
		
		/* Error handler*/
		this.server.use(this.errorHandler)
		
		/* Density server*/
		this.densityserver.default(this.server)
	}
	
	public start(port) {
	
		/* Start listening */
		this.server.listen(port)
		
	}
	
	private errorHandler(err, req, res, next) {
		this.debug(err)
		if (res.headersSent) {
			return next(err)
		}
		res.status(500)
		res.render('error', { error: err })
	}

}

var modules = new Modules()
var httpserver = new HTTPServer()
var ui = "inferno"

if(process.argv[2] == "simple"){
	ui = "simple"
}else if(process.argv[2] == "cosmic"){
	ui = "cosmic"
}
	
httpserver.start(8090)
