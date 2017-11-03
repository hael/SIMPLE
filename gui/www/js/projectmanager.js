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
 getAjax('JSONhandler?function=getprojects', function(data){showProjects(data)});
}

function newProjectPopup () {
	showHideHeaderMenu();
	var newprojectpopup = document.getElementById('newprojectpopup');
	var gauze = document.getElementById('gauze');
	newprojectpopup.style.display = "block";
	gauze.style.display = "block";
}


function hideNewProjectPopup () {
	var newprojectpopup = document.getElementById('newprojectpopup');
	var gauze = document.getElementById('gauze');
	newprojectpopup.style.display = "none";
	gauze.style.display = "none";
}


function showProjects(data){
	var JSONdata = JSON.parse(data);

	var projecttable = document.getElementById('projecttable');
	projecttable.innerHTML = "";
	
	for (var i = 0; i < JSONdata.projects.length; i++) {
		var row = projecttable.insertRow(-1);
		var cell1 = row.insertCell(0);
		var cell2 = row.insertCell(1);
		var cell3 = row.insertCell(1);
		cell1.innerHTML = JSONdata.projects[i].projectname;
		cell1.style.width = "100%";
		cell2.innerHTML = "Edit";
		cell3.innerHTML = "Delete";
	}	
}

function showHideHeaderMenu() {
	var headermenu = document.getElementById('headermenu');
	if (headermenu.style.display == "block") {
		headermenu.style.display = "none";
	} else {
		headermenu.style.display = "block";
	}	
}

getProjects();
