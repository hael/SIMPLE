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

function postAjax (url, success) {
    var xhr = new XMLHttpRequest();
    xhr.open('POST', url, true);
    xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    xhr.onreadystatechange = function() {
        if (xhr.readyState>3 && xhr.status==200){
			success(xhr.responseText);
		}
    };
    xhr.send();
}

function showManualPicker() {
	var jobpane = document.getElementById('jobpane');
	var pickpane = document.getElementById('pickpane');
	jobpane.style.display = "none";
	pickpane.style.display = "inline-flex";
	
	var projecttable = document.getElementById('projecttable').value;
	var projectfolder = document.getElementById('projectfolder').value;
	var jobname = document.getElementById('jobname').value;
	var jobdescription = document.getElementById('jobdescription').value;
	var unbluroutput = document.getElementById('unbluroutput').value;
	var micrographsdirectory = document.getElementById('micrographsdirectory').value;
	var smpd = document.getElementById('smpd').value;
	
	var url = "../JSONhandler?function=manualpick";
	url += "&jobtype=manualpick";
	url += "&table=" + projecttable;
	url += "&jobfolder=" + projectfolder;
	url += "&jobname=" + jobname;
	url += "&jobdescription=" + jobdescription;
	url += "&unbluroutput=" + unbluroutput;
	url += "&micrographsdirectory=" + micrographsdirectory;
	url += "&smpd=" + smpd;
	
	getAjax(url, populateMicrographs);
}

function populateMicrographs(data){
	var JSONdata = JSON.parse(data);
	var rootdirectory = document.getElementById('rootdirectory');
	rootdirectory.value = JSONdata.rootdirectory;
	
	var pickmicrographs = document.getElementById('pickmicrographs');
	
	for (var i = 0; i < JSONdata.snapshots.length; i++) {
		var micrograph = document.createElement("div");
		micrograph.innerHTML = JSONdata.snapshots[i].micrograph;
		micrograph.onclick = function(ev){ev.stopPropagation(); pickMicrograph(this)};
		var keys = Object.keys(JSONdata.snapshots[i]);
		for(var j = 0; j < keys.length; j++) {
			micrograph.setAttribute('data-'+keys[j], JSONdata.snapshots[i][keys[j]]);
		}
		pickmicrographs.appendChild(micrograph);
	}
}

function pickMicrograph(element){
	var canvas = document.getElementById('pickcanvas');
	var context = canvas.getContext('2d');
	var backgroundimage = new Image();
	var rootdirectory = document.getElementById('rootdirectory').value;
	canvas.addEventListener('click', function(event) { refineBox(this, event)}, false);
	canvas.setAttribute('data-micrograph', rootdirectory + "/" + element.getAttribute('data-micrograph'));
	canvas.setAttribute('data-boxfile', element.getAttribute('data-boxfile'));
	backgroundimage.src = "../JPEGhandler?filename=" + rootdirectory + "/" + element.getAttribute('data-micrograph') + "&contrast=" + document.getElementById('pickcontrast').value + "&brightness=" + document.getElementById('pickbrightness').value + "&frameid=0";
	backgroundimage.onload = function(){
		context.clearRect(0, 0, canvas.width, canvas.height);
		canvas.width = backgroundimage.width;
		canvas.height = backgroundimage.height;
		context.drawImage(backgroundimage, 0, 0,backgroundimage.width, backgroundimage.height, 0, 0, canvas.width, canvas.height);
		var url = "../JSONhandler?function=boxfiledata"
		url += "&filename=" + canvas.getAttribute('data-boxfile');
		getAjax(url, drawBoxes);
	};
}

function reloadMicrograph(){
	var canvas = document.getElementById('pickcanvas');
	var context = canvas.getContext('2d');
	var backgroundimage = new Image();
	var rootdirectory = document.getElementById('rootdirectory').value;
	canvas.addEventListener('click', function(event) { refineBox(this, event)}, false);
	backgroundimage.src = "../JPEGhandler?filename=" + canvas.getAttribute('data-micrograph') + "&contrast=" + document.getElementById('pickcontrast').value + "&brightness=" + document.getElementById('pickbrightness').value + "&frameid=0";
	backgroundimage.onload = function(){
		context.clearRect(0, 0, canvas.width, canvas.height);
		canvas.width = backgroundimage.width;
		canvas.height = backgroundimage.height;
		context.drawImage(backgroundimage, 0, 0,backgroundimage.width, backgroundimage.height, 0, 0, canvas.width, canvas.height);
		var url = "../JSONhandler?function=boxfiledata"
		url += "&filename=" + canvas.getAttribute('data-boxfile');
		getAjax(url, drawBoxes);
	};
}

function refineBox(element, event){
    var totalOffsetX = 0;
    var totalOffsetY = 0;
    var canvasX = 0;
    var canvasY = 0;
	
  //  do{
  //      totalOffsetX += element.offsetLeft - element.scrollLeft;
   //     totalOffsetY += element.offsetTop - element.scrollTop;
  //  }
   // while(element = element.offsetParent)

	console.log("offsetleft " + totalOffsetX);
	console.log("scrollleft " + totalOffsetY);
	console.log("eventx " + event.pageX);
	
   // canvasX = event.pageX - totalOffsetX;
  //  canvasY = event.pageY - totalOffsetY;
	canvasX = event.pageX - element.offsetLeft + element.parentNode.scrollLeft;
	canvasY = event.pageY - element.offsetTop + element.parentNode.scrollTop;
	
	var canvas = document.getElementById('pickcanvas');

	var url = "../JSONhandler?function=refinebox"
	url += "&micrograph=" + canvas.getAttribute('data-micrograph');
	url += "&boxsize=" + document.getElementById('boxsize').value;
	url += "&xcoord=" + canvasX;
	url += "&ycoord=" + canvasY;
	url += "&particlediameter=" + document.getElementById('particlediameter').value;
	url += "&boxfile=" + canvas.getAttribute('data-boxfile');
	getAjax(url, drawBoxes);
}

function drawBoxes(data){
	var JSONdata = JSON.parse(data);
	var canvas = document.getElementById('pickcanvas');
	var context = canvas.getContext('2d');
	var boxsize = document.getElementById('boxsize').value;
	var particlediameter = document.getElementById('particlediameter').value;
	
	for(var i = 0; i < JSONdata.boxes.length; i++) {
		context.beginPath();
		context.arc(parseInt(JSONdata.boxes[i][0]) + parseInt(JSONdata.boxes[i][2]) / 2, parseInt(JSONdata.boxes[i][1]) + parseInt(JSONdata.boxes[i][2]) / 2, parseInt(particlediameter) / 2 ,0,2*Math.PI);
		context.lineWidth=8;
		context.strokeStyle=  "green";
		context.stroke();
	}
}
function fileSelect(element) {
	var selectfilepopup = document.getElementById('selectfilepopup');
	var filetarget = document.getElementById('filetarget');
	filetarget.value = element.getAttribute('data-target');
	
	selectfilepopup.style.display = "block";
	var gauze = document.getElementById('gauze');
	gauze.style.display = "block";
	getFileBrowserData();
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
				row.setAttribute("data-target",JSONdata.rootdirectory + "/" + files[i]);
				var cell1 = row.insertCell(0);
				var cell2 = row.insertCell(1);
				cell2.innerHTML = files[i];
				cell2.style.width = "100%";
				row.onclick = function(){this.style.background = "#6698ab"; document.getElementById(document.getElementById('filetarget').value).value=this.getAttribute('data-target')};
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

function showHideControlsPopup () {
	var controlspopup = document.getElementById('controlspopup');
	var gauze = document.getElementById('gauze');
	if (controlspopup.style.display == "block") {
		controlspopup.style.display = "none";
		gauze.style.display = "none";
	} else {
		controlspopup.style.display = "block";
		gauze.style.display = "block";
	}
}


