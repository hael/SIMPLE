class Browser {

  constructor() {
  }

  show(args){
    args['path'] = (args['path'] == "selectedfolder")  ? project.selectedfolder : args['path']
    if(args['path'] == undefined){
		args['path'] = false
	}
    this.args = args
    this.refresh()
  }

  refresh() {
    var request = {
      cls : "browser",
      fnc : "getAll",
      arg : this.args
    }

    if(document.getElementById('filter')){
      request['arg']['filter'] = document.getElementById('filter').value
    } else if (this.args['filter']){
      request['arg']['filter'] = this.args['filter']
    }
    fetcher.fetchJSON(request)
      .then(response => response.json())
      .then ((json) => {
		if(json.html != undefined){
			document.getElementById('browser').style.display = "block"
			document.getElementById('browserwindow').className = "browser"
			document.getElementById('browserwindow').innerHTML = json.html
			document.getElementById('sourcefolder').addEventListener("keyup", (event) => {
			  event.preventDefault()
			  if (event.keyCode === 13) {
				this.stepDown(document.getElementById('sourcefolder').value)
			  }
			})
			document.getElementById('filter').addEventListener("keyup", (event) => {
			  event.preventDefault()
			  if (event.keyCode === 13) {
				this.refresh()
			  }
			})
		}else{
			alert("Error : failed to browse location. Either location doesn't exist or you have insufficient permissions")
		}
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
    this.selection = path.replace('//','/')
    this.args['path'] = path.replace('//','/')
    this.refresh()
  }

  hide() {
    document.getElementById('browser').style.display = "none"
    document.getElementById('browserwindow').className = ""
    document.getElementById('browserwindow').innerHTML = ""
  }

  select(filename, element) {
    this.selection = filename
    for (var row of document.getElementsByClassName('filerowselected')){
      row.className = "filerow"
    }
    element.className = "filerowselected"
  }

  setValue(id){
    document.getElementById(id).value = this.selection
  }

  createFolder(){
    var newfoldername = document.getElementById('newfoldername').value;
    if(newfoldername != ""){
      var request = {
      cls : "browser",
      fnc : "createFolder",
        arg : {
          path : this.args['path'],
          name : newfoldername
        }
      }
      
      fetcher.fetchJSON(request)
      .then(response => response.json())
      .then (json => {
        this.refresh()
      })
    } 
  }
}

var browser
window.addEventListener('load', () => {browser = new Browser()})
