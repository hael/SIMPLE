function selectProject(){
	var projectselector = document.getElementById('projectselector');
	var project = projectselector.options[projectselector.selectedIndex].value;
	if (project == 'createnew'){
		var newprojectpopup = document.getElementById('newprojectpopup');
		showGauze();
		newprojectpopup.style.display = 'block';
	} else {
		projectHistory(project);
	}
}

function addProject(data){
	var projectselector = document.getElementById('projectselector');
	var option = document.createElement("option");
	option.text = getArgument(data, 'projectname');
	option.value = getArgument(data, 'projectname');
	projectselector.add(option);
}

