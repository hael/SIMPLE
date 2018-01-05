backgroundimage = new Image();

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
	var jobfolder = document.getElementById('jobfolder');
	jobfolder.value = JSONdata.jobfolder;
	
	var pickmicrographs = document.getElementById('pickmicrographs');
	
	for (var i = 0; i < JSONdata.snapshots.length; i++) {
		var row = pickmicrographs.insertRow();
		var cell1 = row.insertCell(0);
		cell1.className = "micrograph";
		var micrograph = JSONdata.snapshots[i].micrograph.split("/").pop();
		cell1.innerHTML = micrograph;
		cell1.onclick = function(ev){ev.stopPropagation(); pickMicrograph(this)};
		var keys = Object.keys(JSONdata.snapshots[i]);
		for(var j = 0; j < keys.length; j++) {
			cell1.setAttribute('data-'+keys[j], JSONdata.snapshots[i][keys[j]]);
		}
	}
	var canvas = document.getElementById('pickcanvas');
	canvas.addEventListener('click', function(event) { refineBox(this, event)}, false);
	canvas.oncontextmenu = function(event) { deleteBox(this, event); return false};
}

function pickMicrograph(element){
	var micrographs = document.getElementsByClassName('micrograph');
	for(var i = 0; i < micrographs.length; i++){
		micrographs[i].style.removeProperty("background-color");
	}
	element.style.backgroundColor = "#6698ab";
	var canvas = document.getElementById('pickcanvas');
	var context = canvas.getContext('2d');
//	var backgroundimage = new Image();
	var rootdirectory = document.getElementById('rootdirectory').value;
	//canvas.addEventListener('click', function(event) { refineBox(this, event)}, false);
	canvas.setAttribute('data-micrograph', rootdirectory + "/" + element.getAttribute('data-micrograph'));
	canvas.setAttribute('data-boxfile', element.getAttribute('data-boxfile'));
	backgroundimage.src = "../JPEGhandler?filename=" + rootdirectory + "/" + element.getAttribute('data-micrograph') + "&contrast=" + document.getElementById('pickcontrast').value + "&brightness=" + document.getElementById('pickbrightness').value + "&frameid=0";
	backgroundimage.onload = function(){
		context.clearRect(0, 0, canvas.width, canvas.height);
		canvas.width = backgroundimage.width;
		canvas.height = backgroundimage.height;
		var scale = document.getElementById('scale').value;
		canvas.style.width = canvas.width * scale;
		context.drawImage(backgroundimage, 0, 0,backgroundimage.width, backgroundimage.height, 0, 0, canvas.width, canvas.height);
		var url = "../JSONhandler?function=boxfiledata"
		url += "&filename=" + canvas.getAttribute('data-boxfile');
		getAjax(url, drawBoxes);
	};
}

function reloadMicrograph(){
	var canvas = document.getElementById('pickcanvas');
	var context = canvas.getContext('2d');
	//var backgroundimage = new Image();
	var rootdirectory = document.getElementById('rootdirectory').value;
	//canvas.addEventListener('click', function(event) { refineBox(this, event)}, false);
	backgroundimage.src = "../JPEGhandler?filename=" + canvas.getAttribute('data-micrograph') + "&contrast=" + document.getElementById('pickcontrast').value + "&brightness=" + document.getElementById('pickbrightness').value + "&frameid=0";
	backgroundimage.onload = function(){
		context.clearRect(0, 0, canvas.width, canvas.height);
		canvas.width = backgroundimage.width;
		canvas.height = backgroundimage.height;
		var scale = document.getElementById('scale').value;
		canvas.style.width = canvas.width * scale;
		context.drawImage(backgroundimage, 0, 0,backgroundimage.width, backgroundimage.height, 0, 0, canvas.width, canvas.height);
		document.getElementById('pickcanvasgauze').style.display = "none";
		var url = "../JSONhandler?function=boxfiledata"
		url += "&filename=" + canvas.getAttribute('data-boxfile');
		getAjax(url, drawBoxes);
	};
}

function refineBox(element, event){
	var scale = document.getElementById('scale').value;
    var totalOffsetX = 0;
    var totalOffsetY = 0;
    var canvasX = 0;
    var canvasY = 0;

	canvasX = (event.pageX - element.offsetLeft + element.parentNode.scrollLeft) / scale;
	canvasY = (event.pageY - element.offsetTop + element.parentNode.scrollTop) / scale;
	
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

function deleteBox(element, event){
	var scale = document.getElementById('scale').value;
    var totalOffsetX = 0;
    var totalOffsetY = 0;
    var canvasX = 0;
    var canvasY = 0;

	canvasX = (event.pageX - element.offsetLeft + element.parentNode.scrollLeft) / scale;
	canvasY = (event.pageY - element.offsetTop + element.parentNode.scrollTop) / scale;
	
	var canvas = document.getElementById('pickcanvas');

	var url = "../JSONhandler?function=deletebox"
	url += "&boxsize=" + document.getElementById('boxsize').value;
	url += "&xcoord=" + canvasX;
	url += "&ycoord=" + canvasY;
	url += "&boxfile=" + canvas.getAttribute('data-boxfile');
	getAjax(url, drawBoxes);
}

function drawBoxes(data){
	var gauze = document.getElementById('pickcanvasgauze');
	gauze.style.display = "none";
	var JSONdata = JSON.parse(data);
	var canvas = document.getElementById('pickcanvas');
	var context = canvas.getContext('2d');
	var boxsize = document.getElementById('boxsize').value;
	var particlediameter = document.getElementById('particlediameter').value;
	var particlecount =  document.getElementById('particlecount');
	particlecount.innerHTML = "Particle Count: " + JSONdata.boxes.length;
	canvas.width = backgroundimage.width;
	canvas.height = backgroundimage.height;
	var scale = document.getElementById('scale').value;
	canvas.style.width = canvas.width * scale;
	context.clearRect(0, 0, canvas.width, canvas.height);
	context.drawImage(backgroundimage, 0, 0,backgroundimage.width, backgroundimage.height, 0, 0, canvas.width, canvas.height);
	for(var i = 0; i < JSONdata.boxes.length; i++) {

		if(document.getElementById("pickbox").checked){
			context.beginPath();
			context.rect(parseInt(JSONdata.boxes[i][0]), parseInt(JSONdata.boxes[i][1]), parseInt(boxsize), parseInt(boxsize));
			context.lineWidth=8;
			context.strokeStyle=  "green";
			context.stroke();
		}
		if(document.getElementById("pickparticle").checked){
			context.beginPath();
			context.arc(parseInt(JSONdata.boxes[i][0]) + parseInt(JSONdata.boxes[i][2]) / 2, parseInt(JSONdata.boxes[i][1]) + parseInt(JSONdata.boxes[i][2]) / 2, parseInt(particlediameter) / 2 ,0,2*Math.PI);
			context.lineWidth=8;
			context.strokeStyle=  "green";
			context.stroke();
		}
	}
}

function clearBoxes(){
	var canvas = document.getElementById('pickcanvas');
	var context = canvas.getContext('2d');
	context.clearRect(0, 0, canvas.width, canvas.height);
	context.drawImage(backgroundimage, 0, 0,backgroundimage.width, backgroundimage.height, 0, 0, canvas.width, canvas.height);
	var url = "../JSONhandler?function=clearboxes"
	url += "&filename=" + canvas.getAttribute('data-boxfile');
	
	getAjax(url, drawBoxes);
}

function autoPick(){
	var canvas = document.getElementById('pickcanvas');
	var gauze = document.getElementById('pickcanvasgauze');
	gauze.style.top = canvas.offsetTop;
	gauze.style.left = canvas.offsetLeft;
	gauze.style.display = "block";
	var url = "../JSONhandler?function=autopick"
	url += "&micrograph=" + canvas.getAttribute('data-micrograph');
	url += "&pgrp=" + document.getElementById('pgrp').value;
	url += "&folder=" + document.getElementById('jobfolder').value;
	url += "&smpd=" + document.getElementById('smpd').value;
	if(document.getElementById('pickrefs').value != ""){
		url += "&pickrefs=" + document.getElementById('pickrefs').value;
	}else if(document.getElementById('pickvol').value != ""){
		url += "&pickvol=" + document.getElementById('pickvol').value;
	}
	if(document.getElementById('pcontrast').checked){
		url += "&pcontrast=white";
	}else{
		url += "&pcontrast=black";
	}
	url += "&thres=" + document.getElementById('thres').value;
	url += "&ndev=" + document.getElementById('ndev').value;
	getAjax(url, drawBoxes);
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


