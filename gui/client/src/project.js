class Project {

  constructor(){
    this.updateProjectSelector()
    setInterval(() => {
      if(this.selectedtable != undefined){
        this.refreshHistory()
      }
    }, 10000);

  }
  
  updateProjectSelector(){
	var request = {
      cls : "project",
      fnc : "getProjectSelector",
      arg : {}
    }
    fetcher.fetchJSON(request)
    .then(response => response.json())
    .then (json => {
		document.getElementById('projectselector').innerHTML = json.html
    })  
  }

  selectProject(element){
    if (element.dataset.history == "createproject"){
      this.showCreate()
    }else{
      this.showHistory(element.dataset.history, element.dataset.name, element.dataset.folder, element.dataset.description)
    }
  }

  showHistory(history, name, folder, description){
    this.selectedtable = history
    this.selectedfolder = folder
    this.selectedname = name
    this.selecteddescription = description
    document.getElementById('projectnametitle').innerHTML = "Project - " + name
    document.getElementById('projectnametitle').title = description
    this.refreshHistory()
  }

  refreshHistory(){
    var request = {
      cls : "project",
      fnc : "getProjectHistory",
      arg : {
        history : this.selectedtable,
        name : this.selectedname,
        description : this.selecteddescription
      }
    }
    fetcher.fetchJSON(request)
      .then(response => response.json())
      .then (json => {
		this.drawHistoryTree(json.data)
      })
  }
  
  drawHistoryTree(data){
	var viewwindow = document.getElementById('viewwindow')
	viewwindow.innerHTML = ""
	for(var element of data){
		if(element['parent'] == this.selectedfolder + "/" + this.selectedname + ".simple"){
			var div = document.createElement('div')
			div.className = "tree"
			div.id = element['folder'] + "/" + this.selectedname + ".simple"
			div.appendChild(createBubble(element))
			var branch = document.createElement('div')
			if(localStorage.getItem('branch_' + element.id) != null){
				branch.style.display = 'none'
			}
			branch.className = "branch"
			div.appendChild(branch)
			viewwindow.appendChild(div)
			element['included'] = true
		}
	}
	var updated = true
	while(updated){
		var updated = false;
		var trees = document.getElementsByClassName('tree')
		for(var tree of trees){
			for(var element of data){
				if(element['parent'] == tree.id && element['included'] != true){
					updated = true
					var div = document.createElement('div')
					div.className = "tree"
					div.id = element['folder'] + "/" + this.selectedname + ".simple"
					var stalk = document.createElement('div')//
					stalk.className = "stalk"//
					stalk.appendChild(document.createElement('div'))//
					div.appendChild(stalk)//
					div.appendChild(createBubble(element))
					var branch = document.createElement('div')
					if(localStorage.getItem('branch_' + element.id) != null){
						branch.style.display = 'none'
					}
					branch.className = "branch"
					div.appendChild(branch)
					tree.getElementsByClassName('branch')[0].appendChild(div)
					element['included'] = true
				}
			}
		}
	}
	
	var trees = document.getElementsByClassName('tree')
	for(var tree of trees){
		var branch = tree.getElementsByClassName('branch')[0]
		if(branch.getElementsByClassName('tree').length == 0){
			tree.getElementsByClassName('maxminicon')[0].style.display = 'none'
		}
	}
	
	
	function createBubble(element){
		var bubble = document.createElement('div')
		var topline = document.createElement('div')
		var maxminicon = document.createElement('img')
		
		bubble.id = 'branch_' + element.id
		
		topline.className = "topline"
		topline.innerHTML = element.id + " " + element.type
		
		maxminicon.className = 'maxminicon'
		maxminicon.title = 'Minimise branch'
		maxminicon.onclick = ((event) => {
			var tree = event.target.parentElement.parentElement.parentElement
			var branches = tree.getElementsByClassName('branch')
			var bubble = event.target.parentElement.parentElement
			if(event.target.src.includes('img/minus.png')){
				branches[0].style.display = 'none'
				localStorage.setItem('branch_' + element.id, true)
				event.target.src = 'img/plus.png'
				bubble.getElementsByClassName('lowerline')[0].style.visibility = 'hidden'
				bubble.getElementsByClassName('buttonline')[0].style.visibility = 'hidden'
			}else{
				branches[0].style.display = 'flex'
				localStorage.removeItem('branch_' + element.id)
				event.target.src = 'img/minus.png'
				bubble.getElementsByClassName('lowerline')[0].style.visibility = 'unset'
				bubble.getElementsByClassName('buttonline')[0].style.visibility = 'unset'
			}
		})
		topline.appendChild(maxminicon)
		
		var centerline = document.createElement('div')
		centerline.className = "centerline"
		centerline.innerHTML = element.name
		centerline.title = element.description
		
		var lowerline = document.createElement('div')
		lowerline.className = "lowerline"
		lowerline.innerHTML = element.status
		
		if(element.status == "running"){
			var statusbar = document.createElement('div')
			if(element.view == 'undefined'){
				element.view = 0
			}
			statusbar.title = 'Estimated progress ' + element.view + '%'
			statusbar.className = 'statusbar'
			var status = document.createElement('div')
			status.className = 'status'
			status.style.width = element.view + 'px'
			statusbar.appendChild(status)
			lowerline.appendChild(statusbar)
		}
		
		lowerline.title = "PID : " + element.pid
		var buttonline = document.createElement('div')
		buttonline.className = "buttonline"
		
		if(localStorage.getItem('branch_' + element.id) == null){
			maxminicon.src = 'img/minus.png'		
		}else{
			maxminicon.src = 'img/plus.png'
			lowerline.style.visibility = 'hidden'
			buttonline.style.visibility = 'hidden'
		}
		
		var viewimage = document.createElement('img')
		viewimage.src = "img/view-out.png"
		viewimage.title = "View Output"
		viewimage.onclick = (() => {project.viewSimple(element.folder + "/" + project.selectedname + ".simple")})
		buttonline.appendChild(viewimage)
		
		var viewfiles = document.createElement('img')
		viewfiles.src = "img/folder.png"
		viewfiles.title = "View Files"
		viewfiles.onclick = (() => browser.show({buttons : [{ name : 'view2d', action : 'simpleview.view2D()'}, { name : 'view3d', action : 'simpleview.view3D()'}], path : element.folder, gauze : true}))
		buttonline.appendChild(viewfiles)
		
		var rerunjob = document.createElement('img')
		rerunjob.src = "img/rerun.png"
		rerunjob.title = "Rerun Job"
		rerunjob.onclick = (() => {project.rerunJob(element.arguments)})
		buttonline.appendChild(rerunjob)
		
		if(element.status != 'running'){
			var deletejob = document.createElement('img')
			deletejob.src = "img/delete.png"
			deletejob.title = "Delete Job"
			deletejob.onclick = (() => {project.deleteJob(element)})
			buttonline.appendChild(deletejob)
		}else{
			var killjob = document.createElement('img')
			killjob.src = "img/kill.png"
			killjob.title = "Kill Job"
			killjob.onclick = (() => {project.killJob(element)})
			buttonline.appendChild(killjob)
		}
		var viewlog = document.createElement('img')
		viewlog.src = "img/log.png"
		viewlog.title = "View Log"
		viewlog.onclick = (() => {project.viewLog(element.folder + "/simple.log")})
		buttonline.appendChild(viewlog)
		
		//var stalk = document.createElement('div')
		//bubble.appendChild(stalk)
		bubble.appendChild(topline)
		bubble.appendChild(centerline)
		bubble.appendChild(lowerline)
		bubble.appendChild(buttonline)
		
		if(element.status == "Deleted"){
			buttonline.style.visibility = 'hidden'
			bubble.style.color = 'grey'
		}
		
		bubble.className = "bubble"
		return bubble
	}
  }
  
  showCreate(){
    var request = {
      cls : "project",
      fnc : "getCreateProject",
      arg : {}
    }
    fetcher.fetchJSON(request)
    .then(response => response.json())
    .then (json => {
		if(json.html != undefined){
		  document.getElementById('popup').style.display = "block"
		  document.getElementById('popupwindow').className = "taskinput"
		  document.getElementById('popupwindow').innerHTML = json.html
	    }else{
			alert('Error : unable to open create project popup')
		}
    })
  }

  createProject(){
    var request = {
		cls : "project",
		fnc : "createProject",
		arg : {
			keys : {}
		}
	}
	document.getElementById('keyname').value = document.getElementById('keyname').value.replace(" ", "_")
    var args = document.getElementsByClassName('argument')
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
    fetcher.fetchJSON(request)
    .then(response => {
      return response.json()
    })
    .then(json => {
      if(json.status == 0){
		  alert("Project Created")
		}else{
			alert("Error creating project")
		}
      taskselector.hide()
      this.updateProjectSelector()
     })
  }

	rerunJob(argstring){
		var args = JSON.parse(argstring)
		taskselector.selectedtask = args.type
		return tasksetup.show()
		  .then(() => {
			var keys = Object.keys(args.keys)
			for(var key of keys){
			  if(args.keys[key] != ""){
				var input = document.getElementById(key)
				input.value = args.keys[key]
				//input.parentElement.parentElement.className = 'required'
				input.parentElement.parentElement.style = 'display:table-row'
			  }
			}
			var inputpath = document.getElementById('inputpath')
			for(var i = 0; i < inputpath.length; i++){
				  if(args['parent'].includes(inputpath[i].value)){
					  inputpath[i].selected = true
					  break
				  } 
			}
		  })
	}
	
	viewSimple(projfile){
		var request = {
		  cls : "view",
		  fnc : "getViewSimple",
		  arg : {
			 projfile : projfile 
		  }
		}
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			if(json.html != undefined){
				document.getElementById('popup').style.display = "block"
				document.getElementById('popupwindow').className = "simpleview"
				document.getElementById('popupwindow').innerHTML = json.html
				var viewtypes = document.getElementsByClassName('viewtype')
				if(viewtypes.length > 1){
					viewtypes[viewtypes.length - 2].click()
				}else{
					viewtypes[0].click()
				}
			}else{
				alert('Error : unable to open view')
			}
		})  
	}
	
	viewLog(logfile){
		var request = {
		  cls : "view",
		  fnc : "getViewLog",
		  arg : {
			 logfile : logfile 
		  }
		}
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			if(json.html != undefined){
				document.getElementById('popup').style.display = "block"
				document.getElementById('popupwindow').className = "simpleview"
				document.getElementById('popupwindow').innerHTML = json.html
				var mainpane = document.getElementById('mainpane')
				mainpane.scrollTop = mainpane.scrollHeight - mainpane.clientHeight
			}else{
				alert('Error : logfile does not exist')
			}
		})  
	}
	
	toggleLog(element){
		var mainpane = document.getElementById('mainpane')
		var subpane = document.getElementById('subpane')
		if(mainpane.style.display == 'flex'){
			mainpane.style.display = 'none'
			subpane.style.display = 'flex'
			subpane.scrollTop = subpane.scrollHeight - subpane.clientHeight
			element.innerHTML = 'Master Log'
		}else{
			mainpane.style.display = 'flex'
			subpane.style.display = 'none'
			mainpane.scrollTop = mainpane.scrollHeight - mainpane.clientHeight
			element.innerHTML = 'Subprocess Log'
		}	
	}
	
	deleteProject(){
		var request = {
		  cls : "project",
		  fnc : "deleteProject",
		  arg : {
			history : this.selectedtable,
		  }
		}
		if (confirm("Are you sure you wish to delete this project?")) {
			if (confirm("Are you doubly certain that you wish to remove all your hard work? Please note that this does not remove your project folder. If you wish to remove this you will have to do it manually")) {
				fetcher.fetchJSON(request)
					.then(response => response.json())
					.then (json => {
						location.reload()
					})	
			}
		}
	}
	
	
	deleteJob(element){
		var request = {
		  cls : "task",
		  fnc : "delete",
		  arg : {
			 jobid:element.id,
			 folder:element.folder,
			 projectfolder:this.selectedfolder,
			 history : this.selectedtable,
			 removefolder:false
		  }
		}
		if (confirm("Are you sure you wish to delete job " + element.id + "?")) {
			if (confirm("Do you wish to move the job folder to the trash folder?")) {
				request.arg.removefolder = true
			}
			fetcher.fetchJSON(request)
			.then(response => response.json())
			.then (json => {
				this.refreshHistory()
			})  
		} 
	}
	
	killJob(element){
		var request = {
		  cls : "task",
		  fnc : "kill",
		  arg : {
			 pid:element.pid,
			 id:element.id,
			 projecttable:this.selectedtable
		  }
		}
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			this.refreshHistory()
		})  
	}
	
		

}

var project

window.addEventListener('load', () => {project = new Project()});
