const child = require('child_process')
const fs = require('fs-extra')
const commandExistsSync = require('command-exists').sync
const process = require('process')
const path = require('path')

// Required for pkg packaging
const express = require('express')
const morgan = require('morgan')
const compression = require('compression')
const debug = require('debug')('inferno:HTTPserver')
const auth = require('http-auth')
const pug = require('pug')
const bodyParser = require('body-parser')
const sqlite3 = require('sqlite3').verbose()
const spawn = require('child-process-promise').spawn
const grepit = require('grepit')
const running = require('is-running')
const sharp = require('sharp');
const {getHeader, toPixels} = require('mrchandler')
//

//const HTTPServer = require("./server/httpServer")
//const simpleexec = require('./server/simpleExec')

var multiuserport = 8095

let httpServer
let url

global.appPath = path.dirname(__dirname)
if(!process.env.SIMPLE_PATH){
	console.log("Error: SIMPLE_PATH not set")
	process.exit(1)
}else{
	global.userdata = process.env.SIMPLE_PATH + '/gui_data'
	global.simplepath = process.env.SIMPLE_PATH
}

console.log(global.appPath)
console.log(__dirname)

global.base = __dirname
global.exe =  __filename
global.slurm = false
global.pbs = false

const HTTPServer = require(global.simplepath + '/gui_data/server/src/httpServer')
const simpleexec = require(global.simplepath + '/gui_data/server/src/simpleExec')

fs.ensureDirSync(global.userdata)
fs.ensureFileSync(global.userdata + '/users.htpasswd')
fs.ensureDirSync(global.simplepath + '/tutorials')
fs.ensureDirSync(global.simplepath + '/tutorials/public')

console.log("Application data location ", global.userdata)

if (!commandExistsSync('simple_exec')){
	console.log("Error: simple_exec not found")
	process.exit(1)
}

if (!commandExistsSync('simple_private_exec')){
	console.log("Error: simple_private_exec not found")
	process.exit(1)
}

if (!commandExistsSync('gs')){
	console.log("Error: ghostscript command gs not found")
	process.exit(1)
}

if(commandExistsSync('sbatch')){
	console.log("Found Slurm submission environment")
	global.slurm = true
}

if(commandExistsSync('qsub')){
	console.log("Found PBS/Torque submission environment")
	global.pbs = true
}


function startMultiuserServer() {
  console.log("Starting multiuser server")
  
  var jsonui = child.spawnSync("simple_private_exec",  ["prg=write_ui_json"], {cwd: global.userdata})

  if (jsonui.status == 0){
		console.log("Generated json ui")
  }else{
		console.log("Error: failed to generate json ui")
		process.exit(1)
  }
	
  console.log('Server pid ' + process.pid)
  fs.writeFileSync(global.userdata + '/server.pid', process.pid.toString())
  httpServer = new HTTPServer({multiuser:true})
  httpServer.start(multiuserport)
  url = "http://localhost:" + multiuserport
}


switch(process.argv[2]){
  case "start":
	startMultiuserServer()
    break
  case 'execute':
	console.log('execute')
	var arg = JSON.parse(process.argv[3])
	console.log(arg)

	if(arg['executable'] == "all"){
		console.log("Executable all modified to simple_exec")
		arg['executable'] = "simple_exec"
	}

	if(arg['executable'] == "simple_exec"){
		return simpleexec.exec(arg, arg['jobid'])
		.then(() => {
			console.log('RUNNER Closing')
			process.exit(0)
		})
		.catch((err) => {
			console.log('RUNNER ERROR: ' + err)
			setTimeout(() => {
				process.exit(1)
			}, 3000)

		})
	}else if(arg['executable'] == "simple_distr_exec"){
		return simpleexec.distrExec(arg, arg['jobid'])
		.then(() => {
			console.log('RUNNER Closing')
			process.exit(0)
		})
		.catch((err) => {
			console.log('RUNNER ERROR: ' + err)
			setTimeout(() => {
				process.exit(1)
			}, 3000)
		})
	}else{
		return simpleexec.createDir(arg, arg['jobid'])
		.then(() => {
			console.log('RUNNER Closing')
			process.exit(0)
		})
		.catch((err) => {
			console.log('RUNNER ERROR: ' + err)
			setTimeout(() => {
				process.exit(1)
			}, 3000)
		})
	}
	break
  default:
    console.log('No argument given')
}
