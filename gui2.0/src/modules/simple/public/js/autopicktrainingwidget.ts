class AutoPickTrainingWidget {

	private widgetpopup
	private img
	
	constructor(){
		this.widgetpopup = document.getElementById('widgetpopup')
		this.img = new Image()
		this.img.addEventListener('load',() => {
			this.ctx.drawImage(this.img,0,0)
		}, false)
	}
	
	
	
	refresh(){
		var request = {
			mod : "simple",
			fnc : "getAutopickTrainingWidget",
			arg : {}
		}
		request['arg']['projfile'] = (<HTMLInputElement>document.getElementById('inputpath')).value
		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				this.widgetpopup.innerHTML = json.html
				this.widgetpopup.className = "guinierwidgetpopup"; 
				(<HTMLInputElement>document.getElementById('autopickrefs')).value = (<HTMLInputElement>document.getElementById('keyrefs')).value
			})
	}
	
	view(element){
		this.refresh()
	}
	
	hide(){
		this.widgetpopup.innerHTML = ""
		this.widgetpopup.className = ""
	}
	
	pick(){
		var request = {}
		request['fnc'] = "trainAutopick"
		request['mod'] = "simple"
		request['arg'] = {}
		request['arg']['projectfolder'] = projectselector.selectedfolder
		request['arg']['projfile'] = (<HTMLInputElement>document.getElementById('inputpath')).value
		request['arg']['micrographid'] = (<HTMLSelectElement>document.getElementById('micrographselector')).value
		request['arg']['refs'] = (<HTMLInputElement>document.getElementById('autopickrefs')).value
		request['arg']['ndev'] = (<HTMLInputElement>document.getElementById('autopickquality')).value
		request['arg']['thres'] = (<HTMLInputElement>document.getElementById('autopickthreshold')).value

		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				(<HTMLInputElement>document.getElementById('autopickrefs')).value = json.pickrefs
			})
	}

	done(){
		
	}
}

var autopicktrainingwidget
window.addEventListener('load', () => {autopicktrainingwidget = new AutoPickTrainingWidget()})
