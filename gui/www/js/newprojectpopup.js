function hideNewProjectPopup(){
	var newprojectpopup = document.getElementById('newprojectpopup');
	newprojectpopup.style.display = 'none';
	clearNewProjectPopupLines();
}

function clearNewProjectPopupLines(){
	var newprojectname = document.getElementById('newprojectname');
	newprojectname.value = "";
	var newprojectfolder = document.getElementById('newprojectfolder');
	newprojectfolder.value = "";
}

function confirmNewProjectPopup(){
	var newprojectname = document.getElementById('newprojectname').value;
	var newprojectfolder = document.getElementById('newprojectfolder').value;
	var commandstring = 'cmd=addproject usr=local project=' + newprojectname + ' projectdir=' + newprojectfolder; 
	ws.send(commandstring);
	hideGauze();
	hideNewProjectPopup();
	location.reload(); 
}

function cancelNewProjectPopup(){
	hideNewProjectPopup();
	hideGauze();
}
