class FileSelector {
	
	public selection
	public target
	public action

	private selector
	private gauze
	private updategauze
	
	constructor(){
		this.selector = document.getElementById('fileselector')
		this.gauze = document.getElementById('gauze')
	}
	
	public callback (){
		console.log(this.selection)
	}
	
	public show(update, folder="/tmp"){
		if(update){
			this.updategauze = true
		}else{
			this.updategauze = false
		}
		this.gauze.style.display = "block"
		this.selector.style.display = "block"
		this.request(folder)
	}
	
	public showFolder(update){
		if(update){
			this.updategauze = true
		}else{
			this.updategauze = false
		}
		this.gauze.style.display = "block"
		this.selector.style.display = "block"
		this.folderRequest("/tmp")
	}
	
	public hide(){
		this.selector.style.display = "none"
		this.selector.innerHTML = ""
		if(this.updategauze){
			this.gauze.style.display = "none"
		}
	}
	
	public request (path) {
		var filter = document.getElementById('filter') as HTMLInputElement
		var request
		if(filter != null && filter.value != ""){
			request = {mod : "core", fnc : "folderListing", arg : { pth : path, flt : filter.value, viw : "fileselector" }}
		}else{
			request = {mod : "core", fnc : "folderListing", arg : { pth : path, viw : "fileselector" }}
		}
		postAjax(request, this.draw.bind(this));
	}
	
	public folderRequest (path) {
		var request
		request = {mod : "core", fnc : "folderListing", arg : { pth : path, viw : "folderselector" }}
		postAjax(request, this.draw.bind(this));
	}
	
	public requestUp (path) {
		var newpath = path.split('/')
		newpath.pop()
		var newpath = newpath.join('/')
		if (newpath == ""){
			newpath = "/"
		}
		this.request(newpath)
	}
	
	public requestFolderUp (path) {
		var newpath = path.split('/')
		newpath.pop()
		var newpath = newpath.join('/')
		if (newpath == ""){
			newpath = "/"
		}
		this.folderRequest(newpath)
	}
	
	public draw (data) {
		this.selector.innerHTML = data
		var sourcefolder = document.getElementById('sourcefolder') as HTMLInputElement
		var filter = document.getElementById('filter') as HTMLInputElement
		sourcefolder.addEventListener("keyup", function(event) {
			event.preventDefault()
			if (event.keyCode === 13) {
				fileselector.request(sourcefolder.value)
			}
		})
		
		filter.addEventListener("keyup", function(event) {
			event.preventDefault()
			if (event.keyCode === 13) {
				fileselector.request(sourcefolder.value)
			}
		})
	}
	
	public select (filename, element) {
		this.selection = filename
		for (var row of document.getElementsByClassName('filerow')){
			var rowHTML = row as HTMLElement
			rowHTML.style.backgroundColor = ""
		}
		element.style.backgroundColor = "rgb(134, 215, 113)"
	}
	
	public setValue(){
		var element = document.getElementById(this.target) as HTMLInputElement
		element.value = this.selection
	}
	
	public performAction(){
		this.action()
		this.hide()
	}
	
}

var fileselector = new FileSelector()
