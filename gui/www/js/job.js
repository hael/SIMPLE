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

function showJobPage(jobpage){
	var jobiframe = document.getElementById('jobiframe');
	
	jobiframe.src = "jobs/" + jobpage;
}

function showJobs(element, selector){
	var modeimages = document.getElementsByClassName('modeimage');
	for(var i = 0; i < modeimages.length; i++){
		modeimages[i].style.borderColor = "white";
	}
	
	element.style.borderColor = "orange";
	
	var jobs = document.getElementsByClassName('job');
	
	for(var i = 0; i < jobs.length; i++){
		jobs[i].style.display = "none";
	}
	
	var showjobs = document.querySelectorAll('[' + selector + '="yes"]');
	
	for(var i = 0; i < showjobs.length; i++){
		showjobs[i].style.display = "block";
	}
}

function setProjectTable(){
	var projectselector = parent.parent.document.getElementById('projectselector');
	var projecttable = document.getElementById('projecttable');
	projecttable.value = projectselector.options[projectselector.selectedIndex].getAttribute('data-projecttable');
}

function showHistory(){
	var largepopupiframe = parent.parent.document.getElementById('largepopupiframe');
	largepopupiframe.src = "projecthistory.html";
}

function fileSelect(element) {
	var selectfilepopup = document.getElementById('selectfilepopup');
	selectfilepopup.style.display = "block";
	var gauze = document.getElementById('gauze');
	gauze.style.display = "block";
	getFileBrowserData(element);
}

function folderSelect(element) {
	var selectfilepopup = document.getElementById('selectfilepopup');
	selectfilepopup.style.display = "block";
	var gauze = document.getElementById('gauze');
	gauze.style.display = "block";
	getFolderBrowserData(element);
}

function showFileBrowserData(data){
	var JSONdata = JSON.parse(data);
	var directories = JSONdata.directories;
	var files = JSONdata.files;
	var selectfiledirectory = document.getElementById('selectfiledirectory');
	var selectfiletable = document.getElementById('selectfiletable');
	selectfiletable.innerHTML = "";
	directories.sort();
	var row = selectfiletable.insertRow(-1);
	row.id = JSONdata.directory;
	var cell1 = row.insertCell(0);
	cell1.innerHTML = "<img src=../img/folder.png class=folderimage>";
	var cell2 = row.insertCell(1);
	cell2.innerHTML = "..";
	cell2.style.width = "100%";
	row.ondblclick = function(){getFileBrowserData(this.id)};
	row.onclick = function(){this.style.background = "orange"; document.getElementById('jobfolder').value=this.id};
	for (var i = 0; i < directories.length; i++) {
		var row = selectfiletable.insertRow(-1);
		row.id = JSONdata.directory + "/" + directories[i];
		var cell1 = row.insertCell(0);
		cell1.innerHTML = "<img src=../img/folder.png class=folderimage>";
		var cell2 = row.insertCell(1);
		cell2.innerHTML = directories[i];
		cell2.style.width = "100%";
		row.ondblclick = function(){getFileBrowserData(this.id)};
		row.onclick = function(){this.style.background = "orange"; document.getElementById('jobfolder').value=this.id};
	}
	selectfiledirectory.value = JSONdata.directory;
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
	var rootdir = JSONdata.directory.split('/');
	rootdir.pop()
	row.id = rootdir.join('/');
	var cell1 = row.insertCell(0);
	cell1.innerHTML = "<img src=../img/folder.png class=folderimage>";
	var cell2 = row.insertCell(1);
	cell2.innerHTML = "..";
	cell2.style.width = "100%";
	row.ondblclick = function(){getFolderBrowserData(this.id)};
	row.onclick = function(){this.style.background = "orange"; document.getElementById('jobfolder').value=this.id};
	for (var i = 0; i < directories.length; i++) {
		var row = selectfiletable.insertRow(-1);
		row.id = JSONdata.directory + "/" + directories[i];
		var cell1 = row.insertCell(0);
		cell1.innerHTML = "<img src=../img/folder.png class=folderimage>";
		var cell2 = row.insertCell(1);
		cell2.innerHTML = directories[i];
		cell2.style.width = "100%";
		row.ondblclick = function(){getFolderBrowserData(this.id)};
		row.onclick = function(){this.style.background = "orange"; document.getElementById('jobfolder').value=this.id};
	}
	selectfiledirectory.value = JSONdata.directory;
}

function getFileBrowserData(directory, filter){
	var url = '../JSONhandler?function=listdir';
	if (directory){	
		url += "&directoryname=" + directory;
	}
	if (filter){	
		url += "&filefilter=" + filter;
	}
	getAjax(url, function(data){showFileBrowserData(data)});	
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
	var gauze = document.getElementById('gauze');
	selectfilepopup.style.display = "none";
	gauze.style.display = "none";
}

function showHideModeSelect () {
	var modeimages = document.getElementById('modeimages');
	if(modeimages.style.display == "block"){
		modeimages.style.display = "none";	
	} else {
		modeimages.style.display = "block";
	}
}
