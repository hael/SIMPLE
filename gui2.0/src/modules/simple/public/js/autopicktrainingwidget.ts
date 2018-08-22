class AutoPickTrainingWidget {

	private widgetpopup
	
	constructor(){
		this.widgetpopup = document.getElementById('widgetpopup')
	}
	
	refresh(){
		var request = {
			mod : "simple",
			fnc : "getAutopickTrainingWidget",
			arg : {}
		}
		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				this.widgetpopup.innerHTML = json.html
				this.widgetpopup.className = "guinierwidgetpopup"
			})
	}
	
	view(element){
		this.refresh()
	}
	
	hide(){
		this.widgetpopup.innerHTML = ""
		this.widgetpopup.className = ""
	}
	
	generate(){
		var request = {}
		request['fnc'] = "createAutopickRefs"
		request['mod'] = "simple"
		request['arg'] = {}
		request['arg']['keys'] = {}
		request['arg']['projecttable'] = projectselector.selectedtable
		request['arg']['executable'] = "simple_exec"
		request['arg']['mod'] = "simple"
		request['arg']['type'] = "make_pickrefs"
		request['arg']['projectfolder'] = projectselector.selectedfolder
		request['arg']['name'] = "Generated picking references"
		request['arg']['description'] = "Created from the autopick training widget"
		
		request['arg']['keys']['keyvol1'] = (<HTMLInputElement>document.getElementById('autopickvol')).value
		request['arg']['keys']['keystk'] = (<HTMLInputElement>document.getElementById('autopickclasses')).value
		request['arg']['keys']['keysmpd'] = (<HTMLInputElement>document.getElementById('autopicksmpd')).value
		request['arg']['keys']['keycontrast'] = (<HTMLInputElement>document.getElementById('autopickcontrast')).value
		request['arg']['keys']['keypgrp'] = (<HTMLInputElement>document.getElementById('autopickpgrp')).value

		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				(<HTMLInputElement>document.getElementById('autopickrefs')).value = json.pickrefs
			})
	}

}

var autopicktrainingwidget
window.addEventListener('load', () => {autopicktrainingwidget = new AutoPickTrainingWidget()})
