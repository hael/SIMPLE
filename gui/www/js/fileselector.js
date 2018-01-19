function fileSelect(element) {
	var selectfilepopup = document.getElementById('selectfilepopup');
	var filetarget = document.getElementById('filetarget');
	var filefilter = document.getElementById('filefilter');
	filetarget.value = element.getAttribute('data-target');
	filefilter.value = element.getAttribute('data-filter');
	var projectselector = parent.parent.document.getElementById('projectselector');
	var projectfolder = projectselector.options[projectselector.selectedIndex].getAttribute('data-projectfolder');
	selectfilepopup.style.display = "block";
	var gauze = document.getElementById('gauze');
	gauze.style.display = "block";
	getFileBrowserData(projectfolder);
}

function folderSelect(element) {
	var selectfilepopup = document.getElementById('selectfilepopup');
	var filetarget = document.getElementById('filetarget');
	filetarget.value = element.getAttribute('data-target');
	selectfilepopup.style.display = "block";
	var projectselector = parent.parent.document.getElementById('projectselector');
	var projectfolder = projectselector.options[projectselector.selectedIndex].getAttribute('data-projectfolder');
	var gauze = document.getElementById('gauze');
	gauze.style.display = "block";
	getFolderBrowserData(projectfolder);
}

function showFileBrowserData(data){
	var JSONdata = JSON.parse(data);
	var directories = JSONdata.directories;
	var files = JSONdata.files;
	var selectfiledirectory = document.getElementById('selectfiledirectory');
	var selectfiletable = document.getElementById('selectfiletable');
	selectfiletable.innerHTML = "";
	
	var row = selectfiletable.insertRow(-1);
	var rootdir = JSONdata.rootdirectory.split('/');
	rootdir.pop()
	row.id = rootdir.join('/');
	var cell1 = row.insertCell(0);
	cell1.innerHTML = "<img src=../img/folder.png class=folderimage>";
	var cell2 = row.insertCell(1);
	cell2.innerHTML = "..";
	cell2.style.width = "100%";
	row.ondblclick = function(){getFileBrowserData(this.id)};
	if(!!directories){
		directories.sort();
		for (var i = 0; i < directories.length; i++) {
			if(directories[i][0] != "."){
				var row = selectfiletable.insertRow(-1);
				row.className = "filefolder";
				row.id = JSONdata.rootdirectory + "/" + directories[i];
				var cell1 = row.insertCell(0);
				cell1.innerHTML = "<img src=../img/folder.png class=folderimage>";
				var cell2 = row.insertCell(1);
				cell2.innerHTML = directories[i];
				cell2.style.width = "100%";
				row.ondblclick = function(){getFileBrowserData(this.id)};
			}
		}
	}
	if(!!files){
		files.sort();
		for (var i = 0; i < files.length; i++) {
			if(files[i][0] != "."){
				var row = selectfiletable.insertRow(-1);
				row.className = "filefolder";
				row.setAttribute("data-target",JSONdata.rootdirectory + "/" + files[i]);
				var cell1 = row.insertCell(0);
				var cell2 = row.insertCell(1);
				cell2.innerHTML = files[i];
				cell2.style.width = "100%";
				row.onclick = function(){
										var filefolders = document.getElementsByClassName('filefolder');
										for(var i = 0; i < filefolders.length; i++){
											filefolders[i].style.background = "";
										}
										this.style.background = "#6698ab"; 
										document.getElementById(document.getElementById('filetarget').value).value=this.getAttribute('data-target');
										};
			}
		}
	}
	selectfiledirectory.value = JSONdata.rootdirectory;
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
			row.className = "filefolder";
			row.id = JSONdata.rootdirectory + "/" + directories[i];
			row.setAttribute("data-target",JSONdata.rootdirectory + "/" + directories[i]);
			var cell1 = row.insertCell(0);
			cell1.innerHTML = "<img src=../img/folder.png class=folderimage>";
			var cell2 = row.insertCell(1);
			cell2.innerHTML = directories[i];
			cell2.style.width = "100%";
			row.ondblclick = function(){getFolderBrowserData(this.id)};
			row.onclick = function(){
									var filefolders = document.getElementsByClassName('filefolder');
									for(var i = 0; i < filefolders.length; i++){
										filefolders[i].style.background = "";
									}
									this.style.background = "#6698ab"; 
									document.getElementById(document.getElementById('filetarget').value).value=this.getAttribute('data-target');
									};
		}
	}
	selectfiledirectory.value = JSONdata.rootdirectory;
}

function getFileBrowserData(directory){
	var url = '../JSONhandler?function=listdir';
	if (directory){	
		url += "&directoryname=" + directory;
	}
	var filter = document.getElementById('filefilter').value;
	if (!!filter){	
		url += "&filefilter=" + filter;
	}
	getAjax(url, function(data){showFileBrowserData(data)});	
}

function getFolderBrowserData(directory){
	var url = '../JSONhandler?function=listdir';
	if (directory){	
		url += "&directoryname=" + directory;
	}
	getAjax(url, function(data){showFolderBrowserData(data)});	
}

function hideFileSelect () {
	var selectfilepopup = document.getElementById('selectfilepopup');
	var gauze = document.getElementById('gauze');
	selectfilepopup.style.display = "none";
	gauze.style.display = "none";
}
