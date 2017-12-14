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
		//cell2.innerHTML = "Edit";
		cell2.innerHTML = "";
		cell3.innerHTML = "Delete";
		cell3.setAttribute('data-projectid', JSONdata.projects[i].id);
		cell3.onclick = function(){deleteProject(this)};
	}	
}

function deleteProject(element){
	
	var url = 'JSONhandler?function=deleteproject&projectid=' + element.getAttribute('data-projectid');
	
	if (confirm("Do you wish to delete this project? Please note, you will have to remove the project directory yourself") == true) {
		getAjax(url, function(data){window.location.reload(1)});
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

function folderSelect(element) {
	var selectfilepopup = document.getElementById('selectfilepopup');
	var filetarget = document.getElementById('filetarget');
	filetarget.value = element.getAttribute('data-target');
	selectfilepopup.style.display = "block";
	var gauze = document.getElementById('gauze');
	gauze.style.display = "block";
	getFolderBrowserData();
}

function showFolderBrowserData(data){
	var JSONdata = JSON.parse(data);
	var directories = JSONdata.directories;
	var files = JSONdata.files;
	var selectfiledirectory = document.getElementById('selectfiledirectory');
	var selectfiletable = document.getElementById('selectfiletable');
	selectfiletable.innerHTML = "";
	directories.sort();
	var row = selectfiletable.insertRow(-1);
	var rootdir = JSONdata.rootdirectory.split('/');
	rootdir.pop()
	row.id = rootdir.join('/');
	var cell1 = row.insertCell(0);
	cell1.innerHTML = "<img src=../img/folder.png class=folderimage>";
	var cell2 = row.insertCell(1);
	cell2.innerHTML = "..";
	cell2.style.width = "100%";
	row.ondblclick = function(){getFolderBrowserData(this.id)};
	row.onclick = function(){this.style.background = "#6698ab"; document.getElementById('jobfolder').value=this.id};
	for (var i = 0; i < directories.length; i++) {
		if(directories[i][0] != "."){
			var row = selectfiletable.insertRow(-1);
			row.id = JSONdata.rootdirectory + "/" + directories[i];
			row.setAttribute("data-target",JSONdata.rootdirectory + "/" + directories[i]);
			var cell1 = row.insertCell(0);
			cell1.innerHTML = "<img src=../img/folder.png class=folderimage>";
			var cell2 = row.insertCell(1);
			cell2.innerHTML = directories[i];
			cell2.style.width = "100%";
			row.ondblclick = function(){getFolderBrowserData(this.id)};
			row.onclick = function(){this.style.background = "#6698ab"; document.getElementById(document.getElementById('filetarget').value).value=this.getAttribute('data-target')};
		}
	}
	selectfiledirectory.value = JSONdata.rootdirectory;
}

function getFolderBrowserData(directory, filter){
	var url = '../JSONhandler?function=listdir';
	if (directory){	
		url += "&directoryname=" + directory;
	}
	if (filter){	
		url += "&filefilter=" + filter;
	}
	getAjax(url, function(data){showFolderBrowserData(data)});	
}

function hideFileSelect () {
	var selectfilepopup = document.getElementById('selectfilepopup');
	selectfilepopup.style.display = "none";
}

getProjects();
