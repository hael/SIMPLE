class History2 {

	private popup
	private	gauze

	constructor() {
		this.popup = document.getElementById('popup')
		this.gauze = document.getElementById('gauze')
	}
	
	showLog(folder){
		var request = {
			mod : "core",
			fnc : "showLog",
			arg : {}
		}
		request['arg']['folder'] = folder
		postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				this.gauze.style.display = "block"
				this.popup.innerHTML = json.html
				this.popup.className = "viewlogfile"
			})
	}
	
	hideLog(){
		this.popup.innerHTML = ""
		this.popup.className = "popup"
		this.gauze.style.display = "none"
	}
	
	viewOutput(arg, folder, status){
		var json = JSON.parse(arg)
		var request = {
			mod : json['mod'],
			fnc : json['fnc'],
			arg : {}
		}
		request['arg']['folder'] = folder
		request['arg']['status'] = status
		return postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				viewer.show(json.html)
				if(json.func){
					var tmpFunc = new Function(json.func)
					tmpFunc()
				}
			})
		
	}
	
	resetActions(element){
		element.parentElement.selectedIndex = 0
	}
	
	rerunJob(argstring){
		var args = JSON.parse(argstring)
		taskselector.selectedmodule = args.mod
		taskselector.selectedtask = args.type
		return tasksetup.refresh(args.inputpath)
			.then(() => {
				var keys = Object.keys(args.keys)
				for(var key of keys){
					if(args.keys['keyname'] != ""){
						var input = document.getElementById(key) as HTMLInputElement
						input.value = args.keys[key]
					}
				}
			})
	}
}

var history2
window.addEventListener('load', () => {history2 = new History2()})
