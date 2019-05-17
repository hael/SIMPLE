class ProjectSelector {
	
	private pselement
	private viewwindowbody
	private popup
	private	gauze
	public selectedtable
	public selectedfolder
	
	constructor(){
		this.pselement = document.getElementById('projectselector')
		this.viewwindowbody = document.getElementById('viewwindowbody')
		this.popup = document.getElementById('popup')
		this.gauze = document.getElementById('gauze')
		
		var request = {
			mod : "core",
			fnc : "projectSelector",
			arg : {}
		}
		postAjaxPromise(request)
			.then(response => response.json())
			.then (json => {
				this.pselement.innerHTML = json.html
			})
		setInterval(() => {
			if(this.selectedtable != undefined){
				this.refreshHistory()
			}
		}, 20000);
		
	}
	
	showHistory(history, name, folder, description){
		this.selectedtable = history
		this.selectedfolder = folder
		var request = {
			mod : "core",
			fnc : "projectHistory",
			arg : {history : history}
		}
		postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				document.getElementById('viewwindowtitle').innerHTML = "Project - " + name
				document.getElementById('viewwindowtitle').setAttribute('title', description)
				this.viewwindowbody.innerHTML = json.html
			})
	}
	
	refreshHistory(){
		var request = {
			mod : "core",
			fnc : "projectHistory",
			arg : {history : this.selectedtable}
		}
		postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				this.viewwindowbody.innerHTML = json.html
			})
	}
	
	showCreate(){
		var request = {
			mod : "core",
			fnc : "projectCreator",
			arg : {}
		}
		postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				this.gauze.style.display = "block"
				this.popup.innerHTML = json.html
				this.popup.className = "taskinput"
			})
	}
	
	createProject(){
		var request = {}
		request['fnc'] = "createProject"
		request['mod'] = "core"
		request['arg'] = {}
		request['arg']['keys'] = {}
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
		}).then(function(html) {
			alert("Project Created")
			taskselector.hide()
			var request = {
				mod : "core",
				fnc : "projectSelector",
				arg : {}
			}
			return postAjaxPromise(request)
		})
		.then(response => response.json())
		.then (json => {
			this.pselement.innerHTML = json.html
		})
	}
	
}

var projectselector

window.addEventListener('load', function(){projectselector = new ProjectSelector()})

