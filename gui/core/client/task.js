class TaskSelector {

  constructor() {
  }

  refresh() {
    var request = {
      cls : "task",
      fnc : "getTaskSelector",
      arg : {}
    }
    fetcher.fetchJSON(request)
      .then(response => response.json())
      .then ((json) => {
        document.getElementById('popup').style.display = "block"
        document.getElementById('popupwindow').className = "taskselector"
        document.getElementById('popupwindow').innerHTML = json.html
      })
  }

  show(){
    if(project.selectedtable != undefined){
      this.refresh()
    } else {
      alert("Please select a project first!")
    }
  }
  
  hide() {
	fetcher.killWorker()
    document.getElementById('popup').style.display = "none"
    document.getElementById('popupwindow').className = ""
    document.getElementById('popupwindow').innerHTML = ""
  }

  select(element){
    for (var task of document.querySelectorAll('.taskselected[data-task]')){
      task.className = "taskdisplay"
    }
    element.className = "taskdisplay taskselected"
    document.getElementById('taskdescription').innerHTML =  element.getAttribute('data-description')
    this.selectedtask = element.getAttribute('data-task')
  }  
  
  showHelp(text){
    var taskhelp = document.getElementById('taskhelp')
    taskhelp.style.display = "block"
    var tasktext = document.getElementById('text')
    tasktext.innerHTML = text
  }
  
  hideHelp() {
    var taskhelp = document.getElementById('taskhelp')
    taskhelp.style.display = "none"
    var tasktext = document.getElementById('text')
    tasktext.innerHTML = ""
  }
  
}

var taskselector
window.addEventListener('load', () => {taskselector = new TaskSelector()})


class TaskSetup {

  constructor() {
  }

  show(){
    var request = {
      cls : "task",
      fnc : "getTaskSetup",
      arg : {
		  task : taskselector.selectedtask,
		  projectfolder : project.selectedfolder,
		  projectname : project.selectedname
		  }
    }
    return fetcher.fetchJSON(request)
      .then(response => response.json())
      .then (json => {
        document.getElementById('popup').style.display = "block"
        document.getElementById('popupwindow').className = "taskinput"
        document.getElementById('popupwindow').innerHTML = json.html
      })
  }

  changepage(thispage, newpage){
    var inputpages = document.getElementsByClassName('inputpage')
    var requiredelements = inputpages[thispage].getElementsByClassName('required')
    var change = true
    for(var required of requiredelements){
      if(required.getElementsByClassName('argument')[0].value == ""){
        change = false
      }
    }

    if(change){
      for(var inputpage of inputpages){
        inputpage.style.display = "none"
      }
      inputpages[newpage].style.display = "block"
    }
  }

  startTask(element){
	element.innerHTML = 'Starting'
	element.removeAttribute("onclick")
    var request = {
		cls : "task",
		fnc : "start",
		arg : {
			type : document.getElementById('type').value,
			executable : document.getElementById('fnc').value,
			projectfolder : project.selectedfolder,
			projecttable : project.selectedtable,
			projfile : document.getElementById('inputpath').value,
			parent : document.getElementById('inputpath').value,
			name : document.getElementById('jobname').value,
			description : document.getElementById('jobdescription').value,
			keys : {}
		}
	}
    for(var argument of document.getElementsByClassName('argument')){
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

    fetcher.fetchJSON(request).then(response => {
      return response.json()
    }).then((json) => {
		if(json.status == 'running'){
			alert("Job Started")
		}else{
			alert("There was an error starting your job")
		}
      taskselector.hide()
      project.refreshHistory()
      if(request['arg']['type'] == "manualpick"){
		  var projfile = request['arg']['projectfolder'] + '/' + json['jobid'] + '_manualpick/' + project.selectedname + '.simple'
		  project.viewSimple(projfile)
	  }
    })
  }

  setValue(elementid){
    var element = document.getElementById(elementid)
    element.value = browser.selection
  }

  toggleAdvanced(element){
    var advanced = element.parentElement.getElementsByClassName("advanced")
    var advancedbutton = element.getElementsByClassName("advancedbutton")[0]
    if(advancedbutton.src.includes("minus.png")){
		for(var input of advanced){
			if(input.getElementsByClassName("argument")[0] != undefined && input.getElementsByClassName("argument")[0].value == ''){
				input.style.display = "none"
			}
		}
		advancedbutton.src = "img/plus.png"
	}else{
		for(var input of advanced){
			input.style.display = "table-row"
		}
		advancedbutton.src = "img/minus.png"
	}
  }

}

var tasksetup
window.addEventListener('load', () => {tasksetup = new TaskSetup()})
