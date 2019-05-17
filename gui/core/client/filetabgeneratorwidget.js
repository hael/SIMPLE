class FileTabGeneratorWidget {

  constructor(){
  }

  refresh(){
    var request = {
      cls : "widget",
      fnc : "getFileTabGeneratorWidget",
      arg : {
		  folder : project.selectedfolder,
      projectfolder : project.selectedfolder
		}
    }
    return fetcher.fetchJSON(request)
      .then(response => response.json())
      .then (json => {    
		document.getElementById('widgetpopup').style.display = "block"
		document.getElementById('widgetwindow').className = "filetabgeneratorwidgetpopup"
		document.getElementById('widgetwindow').innerHTML = json.html
      })
  }

  view(element){
    this.callingelement = element
    this.refresh()
  }

  hide(){
	document.getElementById('widgetpopup').style.display = "none"
    document.getElementById('widgetwindow').className = ""
    document.getElementById('widgetwindow').innerHTML = ""
  }

  search(){
	var filetabfolder = document.getElementById('filetabfolder').value
	var filetabfilter = document.getElementById('filetabfilter').value
    var request = {
      cls : "widget",
      fnc : "getFileTabGeneratorWidget",
      arg : {
        folder : filetabfolder,
        filter : filetabfilter,
        projectfolder : project.selectedfolder
      }
    }
    return fetcher.fetchJSON(request)
      .then(response => response.json())
      .then (json => {
		document.getElementById('widgetpopup').style.display = "block"
		document.getElementById('widgetwindow').className = "filetabgeneratorwidgetpopup"
		document.getElementById('widgetwindow').innerHTML = json.html
		document.getElementById('filetabfolder').value = filetabfolder
		document.getElementById('filetabfilter').value = filetabfilter
      })
  }

  selectAll(){
    var elements = document.getElementsByClassName('filecheck')
    for(var element of elements){
      element.checked = true
    }
  }

  save(){
    var request = {
      cls : "widget",
      fnc : "saveFileTabGeneratorWidget",
      arg : {
        filename : document.getElementById('filetaboutputfolder').value + "/" + document.getElementById('filetaboutputfilename').value,
        folder : document.getElementById('filetabfolder').value,
        files : []
      }
    }
    var elements = document.getElementsByClassName('filecheck')
    for(var element of elements){
      if (element.checked){
        request['arg']['files'].push(element.value)
      }
    }
    return fetcher.fetchJSON(request)
      .then(response => response.json())
      .then (json => {
        if(json['status'] == "success"){
          alert("Selection Saved");
          document.getElementById('keyfiletab').value = request['arg']['filename']
          this.hide()
        } else {
          alert("Error saving selection");
        }
      })
  }

}

var filetabgeneratorwidget
window.addEventListener('load', () => {filetabgeneratorwidget = new FileTabGeneratorWidget()})
