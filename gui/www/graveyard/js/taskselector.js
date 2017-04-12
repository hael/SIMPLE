function selectMode(){
	var modeselector = document.getElementById('modeselector');
	var mode = modeselector.options[modeselector.selectedIndex].value;
	if (mode == 'standard'){
		var advancedprograms = document.getElementsByClassName('advancedprogram');
		for (var i = 0; i <  advancedprograms.length; i++){
			advancedprograms[i].style.display = 'none';
		}
	} else if (mode == 'advanced'){
		var advancedprograms = document.getElementsByClassName('advancedprogram');
		for (var i = 0; i <  advancedprograms.length; i++){
			advancedprograms[i].style.display = 'block';
		}
	}
}

function showTask(taskname){
	hideTasks();
	projectHistorySmall();
	var taskpane = document.getElementById(taskname);
	taskpane.style.display = 'block';
}

function hideTasks(){
	var taskpanes = document.getElementsByClassName('taskwindow');
	projectHistoryFull();
	for (var i = 0; i <  taskpanes.length; i++){
		taskpanes[i].style.display = 'none';
	}
}
