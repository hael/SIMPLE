const {app, BrowserWindow, dialog, clipboard, globalShortcut} = require('electron')
const findPort = require('find-open-port')
const fetch = require('node-fetch')
const path = require('path')
const fs = require('fs-extra')
const prompt = require('electron-prompt')
const commandExistsSync = require('command-exists').sync
const child = require('child_process')
const process = require('process')

var multiuserport = 8095

let httpServer
let mainWindow
let url

app.commandLine.appendSwitch('ignore-gpu-blacklist', 'true')

global.base = __dirname
global.appPath = app.getAppPath()
global.userdata = app.getPath('userData')
global.exe = app.getPath('exe')
global.slurm = false
global.pbs = false

if(!process.env.SIMPLE_PATH){
	console.log("Error: SIMPLE_PATH not set")
	process.exit(1)
}else{
	global.simplepath = process.env.SIMPLE_PATH
}

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
  fs.writeFileSync(global.userdata + '/server.pid', process.pid)
  httpServer = new HTTPServer({multiuser:true})
  httpServer.start(multiuserport)
  url = "http://localhost:" + multiuserport
}

function startSingleUser() {
	console.log("Application data location ", global.userdata)
	var jsonui = child.spawnSync("simple_private_exec",  ["prg=write_ui_json"], {cwd: global.userdata})

	if (jsonui.status == 0){
		console.log("Generated json ui")
	}else{
		console.log("Error: failed to generate json ui")
		process.exit(1)
	}
	
	app.on('ready', createWindow)
    app.on('ready', () => {
      globalShortcut.register('CommandOrControl+X', () => {
        mainWindow.webContents.executeJavaScript('tutorial.getHTML()',  result => {
          var filename = Math.floor((Math.random() * 1000000) + 1) + ".html"
          fs.writeFile(global.simplepath + '/tutorials/public/' + filename,result)
          clipboard.writeText(filename)
        })
      })
    })
    app.on('window-all-closed', () => {
      if (process.platform !== 'darwin') {
        app.quit()
      }
    })
    app.on('activate', () => {
      if (mainWindow === null) {
        createWindow()
      }
    })
    app.on('login', function(event, webContents, request, authInfo, callback) {
      event.preventDefault()
      callback(global.username, global.password)
    })
    app.on('ready', () => {
		fetch("http://localhost:" + multiuserport + "/ping")
		.then(() => {
		  console.log("Found multiuser server instance")
		  url = "http://localhost:" + multiuserport
		  return prompt({
			  title : "Login",
			  label : "Username",
			  inputAttrs:{
				type:"text"
			  }
			  
		  })
		  .then(result => {
			if(result === null){
				console.log("User cancelled")
				process.exit(1)
			} else {
				global.username = result
			}
			return prompt({
				  title : "Login",
				  label : "Password",
				  inputAttrs:{
					type:"password"
				  }	  
				})
			})
			.then(result => {
				if(result === null){
					console.log("User cancelled")
					process.exit(1)
				} else {
					global.password = result
				}
				mainWindow.loadURL(url)
			})
			.catch(console.error)
		})
		.catch(err => {
		  console.log("Starting personal server instance")
		  httpServer = new HTTPServer({})
		  findPort()
		  .then(port => {
			console.log("Personal server running on port", port)
			httpServer.start(port)
			url = "http://localhost:" + port
			mainWindow.loadURL(url)
		  })
		})
	})
}

const createWindow = () => {
  mainWindow = new BrowserWindow({
  minWidth: 800,
  minHeight: 800,
  width: 1000,
  height: 1000,
  show: false,
  icon: path.join(global.appPath, 'icons/png/32x32.png'),
  title: "SIMPLE"
  })

  mainWindow.on('closed', () => {
    mainWindow = null
  })

  mainWindow.once('ready-to-show', () => {
    mainWindow.show()
  })

}

switch(process.argv[1]){
  case "multiuser":
	startMultiuserServer()
    break
  case 'execute':
	console.log('execute')
	var arg = JSON.parse(process.argv[2])
	console.log(arg)
	
	if(arg['executable'] == "all"){
		console.log("Executable all modified to simple_exec")
		arg['executable'] = "simple_exec"
	}

	if(arg['executable'] == "simple_exec"){
		return simpleexec.exec(arg, arg['jobid'])
		.then(() => {
			console.log('RUNNER Closing')
			app.quit()
		})
		.catch((err) => {
			console.log('RUNNER ERROR: ' + err)
			setTimeout(() => {
				app.exit(1)
			}, 3000)
		})
	}else if(arg['executable'] == "single_exec"){
		return simpleexec.singleExec(arg, arg['jobid'])
		.then(() => {
			console.log('RUNNER Closing')
			app.quit()
		})
		.catch((err) => {
			console.log('RUNNER ERROR: ' + err)
			setTimeout(() => {
				app.exit(1)
			}, 3000)
		})
	}else{
		return simpleexec.createDir(arg, arg['jobid'])
		.then(() => {
			console.log('RUNNER Closing')
			app.quit()
		})
		.catch((err) => {
			console.log('RUNNER ERROR: ' + err)
			setTimeout(() => {
				app.exit(1)
			}, 3000)
		})
	}
	break
  default:
    startSingleUser()

}
