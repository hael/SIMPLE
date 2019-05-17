class FileTabGeneratorWidget {

	private	plugin
	private widgetpopup
	private callingelement
	
	constructor(){
		this.widgetpopup = document.getElementById('widgetpopup')
	}
	
	refresh(){
		var request = {
			mod : "simple",
			fnc : "getFileTabGeneratorWidget",
			arg : {}
		}
		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				this.widgetpopup.innerHTML = json.html
				this.widgetpopup.className = "filetabgeneratorwidgetpopup"
			})
	}
	
	view(element){
		this.callingelement = element
		this.refresh()
	}
	
	hide(){
		this.widgetpopup.innerHTML = ""
		this.widgetpopup.className = ""
	}
	
	search(){
		var request = {
			mod : "simple",
			fnc : "getFileTabGeneratorWidget",
			arg : {}
		}
		request['arg']['folder'] = (<HTMLInputElement>document.getElementById('filetabfolder')).value
		request['arg']['filter'] = (<HTMLInputElement>document.getElementById('filetabfilter')).value
		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				this.widgetpopup.innerHTML = json.html
				this.widgetpopup.className = "filetabgeneratorwidgetpopup"
			})
	}
	
	selectAll(){
		var elements = document.getElementsByClassName('filecheck')
		for(var element of elements){
			(<HTMLInputElement>element).checked = true
		}
	}
	
	save(){
		var request = {
			mod : "simple",
			fnc : "saveFileTabGeneratorWidget",
			arg : {}
		}
		request['arg']['filename'] = (<HTMLInputElement>document.getElementById('filetaboutputfolder')).value + "/" + (<HTMLInputElement>document.getElementById('filetaboutputfilename')).value
		request['arg']['folder'] = (<HTMLInputElement>document.getElementById('filetabfolder')).value
		request['arg']['files'] = []
		var elements = document.getElementsByClassName('filecheck')
		for(var element of elements){
			if ((<HTMLInputElement>element).checked){
				request['arg']['files'].push((<HTMLInputElement>element).value)
			}
		}
		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				if(json['status'] == "success"){
					alert("Selection Saved");
					(<HTMLInputElement>document.getElementById('keyfiletab')).value = request['arg']['filename']
					this.hide()
				} else {
					alert("Error saving selection");
				}
			})
	}

}

var filetabgeneratorwidget
window.addEventListener('load', () => {filetabgeneratorwidget = new FileTabGeneratorWidget()})
