const {app, BrowserWindow} = require('electron')
const findPort = require('find-open-port')
const HTTPServer = require("./server/httpserver")
const Modules = require("./server/modules")

let modules
let httpServer
let mainWindow
let url

const createServer = () => {
	httpServer = new HTTPServer({
		
	})
	findPort()
	.then((port) => {
      httpServer.start(port)
      url = "http://localhost:" + port
	})	
}

const createWindow = () => {
  mainWindow = new BrowserWindow({
    width: 800,
    height: 600,
    titleBarStyle: 'hidden',
	show: false
  })

  mainWindow.setMenu(null)

  mainWindow.loadURL(url)

  mainWindow.webContents.openDevTools()
  
  mainWindow.on('closed', () => {
    mainWindow = null
  })
  
  mainWindow.once('ready-to-show', () => {
	mainWindow.show()
  })
  
}

console.log("Loading modules")
modules = new Modules()


console.log("Starting server")

createServer()

console.log("Opening browser window")

app.on('ready', createWindow)

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

