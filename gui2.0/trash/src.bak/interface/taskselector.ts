class TaskSelector {
	
	public task
	public module

	private selector
	private gauze
	private viewerbody
	private viewertitle

	constructor(){
		this.selector = document.getElementById('taskselector')
		this.gauze = document.getElementById('gauze')
		this.viewerbody = document.getElementById('viewerbody')
		this.viewertitle = document.getElementById('viewertitle')
	}
	
	public show(){
		this.gauze.style.display = "block"
		this.selector.style.display = "block"
		this.requestTasks()
	}
	
	public hide(){
		this.selector.style.display = "none"
		this.selector.innerHTML = ""
		this.gauze.style.display = "none"
	}
	
	public requestTasks () {
		var request = {}
		request['mod'] = "core"
		request['fnc'] = "taskListing"
		request['arg'] = {}
		postAjax(request, this.draw.bind(this));
	}

	public draw (data) {
		this.selector.innerHTML = data
	}
	
	public showModule(element, module){
		this.module = module
		document.getElementById('modulehelp').innerHTML = ""
		var modules = document.getElementById('moduleselect').getElementsByTagName('tr')
		for(var mod of modules){
			mod.style.backgroundColor = ""
		}
		element.style.backgroundColor = "rgb(134, 215, 113)"
		var tasks = document.getElementById('moduletasks').getElementsByTagName('tr')
		for(var task of tasks){
			task.style.display = "none"
		}
		for(var task of tasks){
			if(task.dataset.module == module){
				task.style.display = "block"
			}
		}
	}
	
	public showTask(element, selectedtask, help){
		var tasks = document.getElementById('moduletasks').getElementsByTagName('tr')
		for(var task of tasks){
			task.style.backgroundColor = ""
		}
		element.style.backgroundColor = "rgb(134, 215, 113)"
		document.getElementById('modulehelp').innerHTML = help
		this.task = selectedtask
	}
	
	public changepage(pageid){
		var inputpages = document.getElementsByClassName('inputpage') as HTMLCollectionOf<HTMLDivElement>
		for(var inputpage of inputpages){
			inputpage.style.display = "none"
		}
		inputpages[pageid].style.display = "block"
	}
	
	public select(){
		this.hide()
		var request = {}
		request['mod'] = "core"
		request['fnc'] = "taskInput"
		request['arg'] = {}
		request['arg']['module'] = this.module
		request['arg']['task'] = this.task
		postAjaxPromise(request).then(function(response){
			return response.text()
		})
		.then(function(html) {
			taskselector.viewertitle.innerHTML = ""
			taskselector.viewerbody.innerHTML = html
			var widget = document.getElementById('widget')
			if(typeof widget != "undefined"){
			var widgetrequest = {}
				widgetrequest['mod'] = taskselector.module
				widgetrequest['fnc'] = widget.dataset.execute
				postAjaxPromise(widgetrequest).then(function(response){
					return response.json()
				}).then(function(data){
					widget.innerHTML = data['html']
				})
			}
		})
		.catch(function(err){
			alert(err)
		})
	}
	
	public startTask(){
		var request = {}
		var fnc = document.getElementById('fnc') as HTMLInputElement
		var mod = document.getElementById('mod') as HTMLInputElement
		var type = document.getElementById('type') as HTMLInputElement
		request['fnc'] = "execute"
		request['mod'] = mod.value
		request['arg'] = {}
		request['arg']['keys'] = {}
		request['arg']['projecttable'] = project.currentselectedtable
		request['arg']['executable'] = fnc.value
		request['arg']['type'] = type.value
		request['arg']['projectfolder'] = project.currentselectedfolder
		var args = document.getElementsByClassName('argument') as HTMLCollectionOf<HTMLInputElement>
		for(var argument of args){
			if(argument.checked){
				request['arg']['keys'][argument.id] = "true"
			} else {
				request['arg']['keys'][argument.id] = argument.value
			}
		}
		postAjaxPromise(request).then(function(response){
			return response.text()
		}).then(function(html) {
			alert("Job Started")
			taskselector.viewertitle.innerHTML = ""
			taskselector.viewerbody.innerHTML = ""
		})
	}
}

var taskselector = new TaskSelector()
