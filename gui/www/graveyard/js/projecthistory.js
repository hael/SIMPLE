function projectHistory(project){
	var historytable = document.getElementById('historytable');
	historytable.innerHTML = "";
	var commandstring = 'cmd=listjobs usr=local project=' + project; //USER
	ws.send(commandstring);
}

function addJob(data){
	var historytable = document.getElementById('historytable');
	var row = historytable.insertRow();
	var cell1 = row.insertCell();
	cell1.className = "historystatus";
	var cell2 = row.insertCell();
	cell2.className = "historyjobname";
	var cell3 = row.insertCell();
	cell1.innerHTML = "Status";
	cell2.innerHTML = getArgument(data, 'jobname');
	cell3.innerHTML = "<img class=jobdropdownicon src=img/dropdown.png></img>";
	cell3.innerHTML += "<div onclick=pipelineView('"+getArgument(data, 'jobdir')+"')>Pipeline View</div>";

}

function projectHistoryFull(){
	var historyjobnames = document.getElementsByClassName('historyjobname');
	for (var i = 0; i <  historyjobnames.length; i++){
		historyjobnames[i].style.display = 'block';
	}	
}

function projectHistorySmall(){
	var historyjobnames = document.getElementsByClassName('historyjobname');
	for (var i = 0; i <  historyjobnames.length; i++){
		historyjobnames[i].style.display = 'none';
	}
}
