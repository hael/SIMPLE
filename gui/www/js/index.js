function IncludeHTML() {
  var z, i, elmnt, file, xhttp;
  z = document.getElementsByTagName("*");
  for (i = 0; i < z.length; i++) {
    elmnt = z[i];
    file = elmnt.getAttribute("include-html");
    if (file) {
      xhttp = new XMLHttpRequest();
      xhttp.onreadystatechange = function() {
        if (this.readyState == 4 && this.status == 200) {
          elmnt.innerHTML = this.responseText;
          elmnt.removeAttribute("include-html");
          IncludeHTML();
        }
      }
      xhttp.open("GET", file, true);
      xhttp.send();
      return;
    }
  }
}

function showGauze(){
	var gauze = document.getElementById('gauze');
	gauze.style.display = 'block';
}

function hideGauze(){
	var gauze = document.getElementById('gauze');
	gauze.style.display = 'none';
}

function getArgument(arguments, argumentname){
	var argumentarray = arguments.split(' ');
	var returnvalue = "";
	for (var i = 0; i < argumentarray.length; i++) {
		var keyvalue = argumentarray[i].split('=');
		if (keyvalue[0] == argumentname){
			returnvalue = keyvalue[1];
		}
	}
	return returnvalue;
}

ws = new WebSocket('ws://localhost:8080/');

ws.onerror = function(event){
	console.log(event.data);
}

ws.onopen= function(){
   ws.send('cmd=listprojects usr=local');
}

ws.onmessage = function(event) {
	switch(getArgument(event.data, 'cmd')){
		case 'pipelineview':
			addPipelineView(event.data);
			break;
		case 'pipelineviewpage':
			addPipelineViewData(event.data);
			break;
		case 'pipelinemicrographview':
			addPipelineViewMicrograph(event.data);
			break;
		case 'dir':
			addFolderViewData(event.data);
			break;
		case 'listprojects':
			addProject(event.data);
			break;
		case 'listjobs':
			addJob(event.data);
			break;	
		case 'addjob':
			hideTasks();
			var project = document.getElementById('projectselector').value;
			projectHistory(project);
			break;
	};
	return false;
}

ws.onclose= function(){
   console.log('close');
}
