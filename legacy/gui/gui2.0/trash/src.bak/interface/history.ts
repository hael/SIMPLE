class Project {

	public currentselected
	public currentselectedpid
	public currentselectedtable
	public currentselectedfolder

	constructor(){
	}
	
	public show(){
		document.getElementById('projectdiv').style.display = "inline-block"
		document.getElementById('newprojectbutton').style.display = "inline-block"
		this.updateProjects()
	}
	
	public hide(){

	}
	
	public showNew(){
		document.getElementById('gauze').style.display = "block"
		document.getElementById('newproject').style.display = "inline-block"
	}
	
	public hideNew(){
		document.getElementById('newproject').style.display = "none"
		document.getElementById('gauze').style.display = "none"
	}
	
	public createNew(){
		var request = {}
		request['mod'] = "core"
		request['fnc'] = "createProject"
		request['arg'] = {}
		request['arg']['name'] = (<HTMLInputElement>document.getElementById('newprojectname')).value
		request['arg']['description'] = (<HTMLInputElement>document.getElementById('newprojectdescription')).value
		request['arg']['folder'] = (<HTMLInputElement>document.getElementById('newprojectfolder')).value
		postAjaxPromise(request).then(function(response){
			return response.json()
		})
		.then(function(json) {
			if(json['status'] == "success"){
				alert("Project successfully created")
				project.hideNew()
				project.updateProjects()
			}
			
		})
	}
	
	public updateProjects(){
		var request = {}
		request['mod'] = "core"
		request['fnc'] = "listProjects"
		request['arg'] = {}
		postAjaxPromise(request).then(function(response){
			return response.text()
		})
		.then(function(html) {
			document.getElementById('project').innerHTML=html
		})
	}
	
	public selectProject(pid, name, projecttable, projectfolder){
		project.currentselected = name
		project.currentselectedpid = pid
		project.currentselectedtable = projecttable
		project.currentselectedfolder = projectfolder
		var request = {}
		request['mod'] = "core"
		request['fnc'] = "projectHistory"
		request['arg'] = {}
		request['arg']['projectid'] = pid
		postAjaxPromise(request).then(function(response){
			return response.text()
		})
		.then(function(html) {
			document.getElementById('viewertitle').innerHTML = ""
			document.getElementById('viewerbody').innerHTML = html
		})
	}
	
}

var project = new Project()

window.onload = function(){
	project.show()
}
