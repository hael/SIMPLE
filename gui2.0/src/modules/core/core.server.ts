import * as fs from 'fs'
import * as binaryfile from 'binary-file'
//import * as sharp from 'sharp'
import HeaderButtons from "./headerbuttons"
import Project from "./project"
import Browser from "./browser"
import Task from "./task"
import MRC from "./mrc"
import JPG from "./jpg"
import View from "./view"

class Module {
	
	private headerbuttons
	private project
	private	browser
	private	task
	private mrc
	private view
	private jpg
	
	constructor(){
		process.stdout.write("Loading module - Core ... ")
		this.headerbuttons = new HeaderButtons()
		this.project = new Project()
		this.browser = new Browser()
		this.task = new Task()
		this.mrc = new MRC()
		this.jpg = new JPG()
		this.view = new View()
		process.stdout.write("Done\n")
	}
	
	headerButtons (modules, arg) {
		return this.headerbuttons.get()
	}
	
	projectSelector (modules, arg) {
		return this.project.get(arg['user'])
	}
	
	projectHistory (modules, arg) {
		return this.project.getHistory(arg)
	}
	
	folderListing (modules, arg) {
		return this.browser.getAll(arg)
	}

	taskSelector (modules, arg) {
		return this.task.getSelector(arg)
	}
	
	refineTasks(modules, arg) {
		return this.task.refineSelector(arg)
	}
	
	taskSetup(modules, arg) {
		return this.task.setup(arg)
	}
	
	deleteTask(modules, arg) {
		return this.task.delete(arg)
	}
	
	killTask(modules, arg) {
		return this.task.kill(arg)
	}
	
	showLog(modules, arg) {
		return this.project.getLog(arg)
	}
	
	projectCreator(modules, arg) {
		return this.project.getNew()
	}
	
	createProject(modules, arg) {
		return this.project.createNew(arg)
	}
	
	taskCreate(modules, arg) {
		return this.task.createNew(arg)
	}
	
	updatePid(table, jobid, jobpid) {
		return this.task.updatePid(table, jobid, jobpid)
	}
	
	updateStatus(table, jobid, status) {
		return this.task.updateStatus(table, jobid, status)
	}
	
	readMRCHeader(file) {
		return this.mrc.readHeader(file)
	}
	
	getImage(arg) {
		if(arg['stackfile'].includes(".mrc")){
			return(this.mrc.toJPEG(arg))
		} else if(arg['stackfile'].includes(".jpg")){
			return(this.jpg.toJPEG(arg))
		}
	}
	
	view2D(modules, arg) {
		return this.view.view2d(arg)
	}
	
	view3D(modules, arg) {
		return this.view.view3d(arg)
	}
	
}

module.exports = new Module()


