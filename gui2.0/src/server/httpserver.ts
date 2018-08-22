import * as fs from 'fs'
import Modules from "./modules"

export default class HTTPServer {

	private		server
	private		express
	private		logger
	private		debug
	private		compression
	private		densityserver
	private		auth
	private		basic
	private		modules
	
	constructor (modules, ui, multiuser) {
		
		/* Import dependencies */
		this.express = require('express')
		this.logger = require('morgan')
		this.compression = require('compression')
		this.debug = require('debug')('inferno:HTTPserver')
		this.densityserver = require("../../external/DensityServer/server/web-api")
		this.auth = require('http-auth')
		
		if(multiuser){
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
		}
		
		/* Create server */
		this.server = this.express()
		
		/* Setup server */
		this.server.use(this.logger('dev'))
		this.server.use(this.express.json())
		this.server.use(this.express.urlencoded({ extended: false }))
		this.server.use(this.express.static('public'))
		this.server.use(this.compression({ level: 6, memLevel: 9, chunkSize: 16 * 16384, filter: function () { return true; } }))
		if(multiuser){
			this.server.use(this.auth.connect(this.basic))
		}
		
		/* Setup PUG view engine */
		this.server.set('view engine', 'pug')
		
		/* Setup routes */
		this.server.get('/', function(req, res) {
			res.render('./index.pug', {ui : ui})
		})
		
		this.server.get('/image', (req, res) => {
			modules['available']['core']['getImage'](req.query).then((result) => {
				res.contentType('jpeg')
				res.send(result.image)
				}) // add catch send error
		})
		
		this.server.post('/', function (req, res) {
			var fnc = req.body.fnc
			var mod = req.body.mod
			var arg = req.body.arg
			if(multiuser){
				arg['user'] = req.user
			} else {
				arg['user'] = "simple"
			}
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
