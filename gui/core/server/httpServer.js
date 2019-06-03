const express = require('express')
const morgan = require('morgan')
const compression = require('compression')
const debug = require('debug')('inferno:HTTPserver')
const auth = require('http-auth')
const pug = require('pug')
const bodyParser = require('body-parser');

const densityserver = require("../../external/DensityServer/server/web-api")

const Project = require('./project')
const Task = require('./task')
const View = require('./view')
const Browser = require('./browser')
const Tutorial = require('./tutorial')
const Widget = require('./widget')

class HTTPServer {

  constructor (args) {
	  
	this.project = new Project()
	this.task = new Task()
	this.view = new View()
	this.browser = new Browser()
	this.widget = new Widget()
	this.tutorial = new Tutorial()
	
    var multiuser = args['multiuser'] ? true : false

    var indexpug = pug.compileFile(global.appPath + '/core/client/index.pug')

    if(multiuser){
      this.basic = auth.basic({
        realm: "simple",
        file: global.userdata + "/users.htpasswd"
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

   
    
    var bodyParser = require('body-parser');
	this.server.use(express.json({limit: "100mb"}));
	this.server.use(express.urlencoded({limit: "100mb", extended: false, parameterLimit:1000000}))
	
	this.server.use(morgan('dev'))
    this.server.use(express.json())

    this.server.use(express.static(global.appPath + '/public'))
    this.server.use(compression({ level: 6, memLevel: 9, chunkSize: 16 * 16384, filter: () => { return true; } }))
    
    this.server.get('/ping', (req, res) => {
      res.send("pong")
    })
    
    this.server.use("/tutorialdata", express.static(global.simplepath +'/tutorials/public'))
    this.server.get('/tutorial', (req, res) => {
		this.tutorial[req.query.fnc](req.query.arg).then(result=>{
			res.send(result)
		})
    })
    
    if(multiuser){
      this.server.use(auth.connect(this.basic))
    }

    this.server.set('view engine', 'pug')

    this.server.get('/', (req, res) => {
      res.send(indexpug())
    })

    this.server.get('/image', (req, res) => {
	 if(req.query.stackfile.includes(".mrc")){
      return this.view.getMRCJPEG(req.query)
      .then(result => {
        res.contentType('jpg')
        res.send(result.image)
      })
      .catch(err => {
		  console.log(err)
		  res.send({error : err})
	  })
     }else if(req.query.stackfile.includes(".jpg")){
		return this.view.getJPEG(req.query)
		.then(result => {
			res.contentType('jpg')
			res.send(result.image)
		})  
		.catch(err => {
		  console.log(err)
		  res.send({error : err})
		})
	 }
    })

    this.server.post('/', (req, res) => {
		console.log('here')
      var fnc = req.body.fnc
      var cls = req.body.cls
      var arg = req.body.arg
      if(multiuser){
        arg['user'] = req.user
      } else {
        arg['user'] = "simple"
      }
      console.log(cls, fnc, arg)
      return this[cls][fnc](arg)
		.then(result => {
          res.send(result)
        })
        .catch(err => {
			console.log(err)
          res.send({error : err})
        })
     
    })

    this.server.use(this.errorHandler)

    densityserver.default(this.server)
  }

  start(port) {
    this.server.listen(port)
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
