declare var global:any

import * as fs from 'fs'
import Modules from "./server/modules"
import HTTPServer from "./server/httpserver"

var ui = "inferno"
var multiuser = false

if(process.argv[2] == "simple"){
	ui = "simple"
}else if(process.argv[2] == "cosmic"){
	ui = "cosmic"
}

if(process.argv[3] == "multiuser"){
	multiuser = true
}

global.modules = new Modules()

var httpserver = new HTTPServer(global.modules, ui, multiuser)

httpserver.start(8090)
