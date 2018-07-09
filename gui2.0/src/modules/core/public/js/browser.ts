class Browser {
	
	private popup
	private	gauze
	private args
	public selection
	
	constructor() {
		this.popup = document.getElementById('browserpopup')
		this.gauze = document.getElementById('gauze')
	}
	
	show(args){
		this.args = args
		this.refresh()
	}
	
	refresh() {
		var request = {
			mod : "core",
			fnc : "folderListing",
			arg : this.args
		}
		if(document.getElementById('filter')){
			request['arg']['filter'] = (<HTMLInputElement>document.getElementById('filter')).value
		}
		postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				if(this.args['gauze']){
					this.gauze.style.display = "block"
				}
				this.popup.innerHTML = json.html
				this.popup.className = "browser"
				document.getElementById('sourcefolder').addEventListener("keyup", (event) => {
					event.preventDefault()
					if (event.keyCode === 13) {
						this.stepDown((<HTMLInputElement>document.getElementById('sourcefolder')).value)
					}
				})
				document.getElementById('filter').addEventListener("keyup", (event) => {
					event.preventDefault()
					if (event.keyCode === 13) {
						this.refresh()
					}
				})
			})
	}
	
	stepUp(path) {
		var newpath = path.split('/')
		newpath.pop()
		var newpath = newpath.join('/')
		if (newpath == ""){
			newpath = "/"
		}
		this.args['path'] = newpath
		this.selection = newpath
		this.refresh()
	}
	
	stepDown(path) {
		this.selection = path
		this.args['path'] = path
		this.refresh()
	}
	
	hide() {
		this.popup.innerHTML = ""
		this.popup.className = "popup"
		if(this.args['gauze']){
			this.gauze.style.display = "none"
		}
	}
	
	select(filename, element) {
		this.selection = filename
		for (var row of document.getElementsByClassName('filerowselected')){
			row.className = "filerow"
		}
		element.className = "filerowselected"
	}
	
	setValue(id){
		(<HTMLInputElement>document.getElementById(id)).value = this.selection
	}
}

var browser
window.addEventListener('load', () => {browser = new Browser()})
