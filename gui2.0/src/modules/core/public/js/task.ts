class TaskSelector {
	
	private popup
	private	gauze
	public selectedtask
	public selectedmodule

	constructor() {
		this.popup = document.getElementById('popup')
		this.gauze = document.getElementById('gauze')
	}
	
	refresh() {
		var request = {
			mod : "core",
			fnc : "taskSelector",
			arg : {}
		}
		if(projectselector.selectedfolder != undefined){
			request['arg']['folder'] = projectselector.selectedfolder
		}else{
			request['arg']['folder'] = "/tmp"
		}
		postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				this.gauze.style.display = "block"
				this.popup.innerHTML = json.html
				this.popup.className = "taskselector"
			})
	}
	
	show(){
		if(projectselector.selectedtable != undefined){
			this.refresh()
		} else {
			alert("Please select a project first!")
		}
	}

	showHelp(text){
		var taskhelp = document.getElementById('taskhelp') as HTMLElement
		taskhelp.style.display = "block"
		var tasktext = document.getElementById('text') as HTMLElement
		tasktext.innerHTML = text
	}
	
	hide() {
		this.popup.innerHTML = ""
		this.popup.className = "popup"
		this.gauze.style.display = "none"
	}
	
	hideHelp() {
		var taskhelp = document.getElementById('taskhelp') as HTMLElement
		taskhelp.style.display = "none"
		var tasktext = document.getElementById('text') as HTMLElement
		tasktext.innerHTML = ""
	}
	
	refineTasks(){
		(<HTMLInputElement>document.getElementById('inputpath')).value = browser.selection
		var request = {
			mod : "core",
			fnc : "refineTasks",
			arg : {path : browser.selection}
		}
		postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				for(var module in json['refinedselection']){
					for(var task of json['refinedselection'][module]){
						var search = "[data-module=" + module + "][data-task=" + task + "]"
						for (var element of document.querySelectorAll(search)){
							element.className = "taskdisplay"
						}
					}
				}
			})
	}
	
	showAllTasks(module){
		for (var task of document.querySelectorAll('[data-module]')){
			if(task.className == "advanced" || task.className == "advanceddisplay" || task.className == "advancedselected"){
				task.className = "task"
			}
		}
		//for(var task of document.querySelectorAll('[data-module='+ module + ']')){
		for(var task of document.querySelectorAll('[data-module]')){
			task.className = "taskdisplay"
		}
	}
	
	select(element){
		for (var task of document.querySelectorAll('[data-task]')){
			if(task.className == "taskselected"){
				task.className = "taskdisplay"
			}
		}
		element.className = "taskselected"
		document.getElementById('taskdescription').innerHTML =  element.getAttribute('data-description')
		this.selectedtask = element.getAttribute('data-task')
		this.selectedmodule = element.getAttribute('data-modulename')
	}
}

var taskselector
window.addEventListener('load', () => {taskselector = new TaskSelector()})

class TaskSetup {

	private popup
	private	gauze
	private advanced

	constructor() {
		this.popup = document.getElementById('popup')
		this.gauze = document.getElementById('gauze')
	}
	
	refresh(input) {
		var request = {
			mod : "core",
			fnc : "taskSetup",
			arg : {}
		}
		request['arg']['module'] = taskselector.selectedmodule
		request['arg']['task'] = taskselector.selectedtask
		request['arg']['inputpath'] = input
		return postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				this.gauze.style.display = "block"
				this.popup.innerHTML = json.html
				this.popup.className = "taskinput"
				this.advanced = false
			})
	}
	
	show(){
		this.refresh((<HTMLInputElement>document.getElementById('inputpath')).value)
	}
	
	changepage(thispage, newpage){
		var inputpages = document.getElementsByClassName('inputpage') as HTMLCollectionOf<HTMLDivElement>
		var requiredelements = inputpages[thispage].getElementsByClassName('required')
		var change = true
		for(var required of requiredelements){
			if((<HTMLInputElement>required.getElementsByClassName('argument')[0]).value == ""){
				change = false
			}
		}
		
		if(change){
			for(var inputpage of inputpages){
				inputpage.style.display = "none"
			}
			inputpages[newpage].style.display = "block"
		}
	}
	
	startTask(){
		var request = {}
		var fnc = document.getElementById('fnc') as HTMLInputElement
		var mod = document.getElementById('mod') as HTMLInputElement
		var type = document.getElementById('type') as HTMLInputElement
		var inputpath = document.getElementById('inputpath') as HTMLInputElement
		var jobname = document.getElementById('jobname') as HTMLInputElement
		var jobdescription = document.getElementById('jobdescription') as HTMLInputElement
		request['fnc'] = "execute"
		request['mod'] = mod.value
		request['arg'] = {}
		request['arg']['keys'] = {}
		request['arg']['projecttable'] = projectselector.selectedtable
		request['arg']['executable'] = fnc.value
		request['arg']['mod'] = mod.value
		request['arg']['type'] = type.value
		request['arg']['projectfolder'] = projectselector.selectedfolder
		request['arg']['inputpath'] = inputpath.value
		request['arg']['name'] = jobname.value
		request['arg']['description'] = jobdescription.value
		var args = document.getElementsByClassName('argument') as HTMLCollectionOf<HTMLInputElement>
		for(var argument of args){
			if(argument.parentElement.parentElement.className == "required" && argument.value == ""){
				argument.style.borderColor = "red"
				setTimeout(() => { 
					argument.style.borderColor = null
				}, 1000)
				return
			}
			if(argument.checked){
				request['arg']['keys'][argument.id] = "true"
			} else {
				request['arg']['keys'][argument.id] = argument.value
			}
		}
		postAjaxPromise(request).then(function(response){
			return response.text()
		}).then((html) => {
			alert("Job Started")
		}).then(() =>{
			taskselector.hide()
			projectselector.refreshHistory()
		})
	}
	
	setValue(elementid){
		var element = document.getElementById(elementid) as HTMLInputElement
		element.value = browser.selection
	}
	
	toggleAdvanced(element){
		var advanced = element.parentElement.parentElement.getElementsByClassName("advanced")
		if(this.advanced){
			for(var input of advanced){
				input.style.display = "none"
			}
			this.advanced = false
		} else {
			for(var input of advanced){
				input.style.display = "table-row"
			}
			this.advanced = true
		}
	}

}

var tasksetup
window.addEventListener('load', () => {tasksetup = new TaskSetup()})
