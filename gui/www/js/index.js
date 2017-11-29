function getAjax (url, success) {
    var xhr = new XMLHttpRequest();
    xhr.open('GET', url, true);
    xhr.onreadystatechange = function() {
        if (xhr.readyState>3 && xhr.status==200){
			success(xhr.responseText);
		}
    };
    xhr.send();
}

function getProjects () {
 getAjax('JSONhandler?function=getprojects', function(data){addProjects(data)});
}

function addProjects(data){
	var JSONdata = JSON.parse(data);

	var projectselector = document.getElementById('projectselector');
	
	for (var i = 0; i < JSONdata.projects.length; i++) {
		var option = document.createElement("option");
		option.setAttribute('data-projectname', JSONdata.projects[i].projectname);
		option.setAttribute('data-projectdescription', JSONdata.projects[i].projectdescription);
		option.setAttribute('data-projectfolder', JSONdata.projects[i].projectfolder);
		option.setAttribute('data-projecttable', JSONdata.projects[i].projecttable);
		option.setAttribute('data-projectusage', JSONdata.projects[i].projectusage);
		option.innerHTML = JSONdata.projects[i].projectname;
		projectselector.add(option); 
	}
}
function showFileViewer(){
 var mainpaneiframe = document.getElementById('mainpaneiframe');
 mainpaneiframe.src = "fileviewer.html";
}

function showProjectManager(){
 var mainpaneiframe = document.getElementById('mainpaneiframe');
 mainpaneiframe.src = "projectmanager.html";
}

function showProjecthistory(){
  var mainpaneiframe = document.getElementById('mainpaneiframe');
 mainpaneiframe.src = "projecthistory.html";
}

function showJobSelector (){
 var mainpaneiframe = document.getElementById('mainpaneiframe');
 mainpaneiframe.src = "jobselector.html";
}

function selectProject (){
	var project = document.getElementById('projectselector').value;
	if (project == "projectmanager") {
		showProjectManager();
	} else if (project != ""){
		showProjecthistory();
	}
}

function showHideMainMenu (){
	var menu = document.getElementById('hoverdiv');
	if(menu.style.display == "block"){
		menu.style.display = "none";
	} else {
		menu.style.display = "block";
	}
}


getProjects();
