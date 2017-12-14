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

function getParticleViewData () {
 getAjax('JSONhandler?function=extractview&' + window.location.toString().split('?')[1], function(data){viewParticleViewInit(data)});
}

function viewParticleViewInit (data) {
	var JSONdata = JSON.parse(data);
	
	var micrographs = JSONdata.iterations;
	var rootdirectory = document.getElementById('rootdirectory');
	rootdirectory.value = JSONdata.rootdirectory;
	var particlecount = document.getElementById('particlecount');
	particlecount.innerHTML = JSONdata.particlecount;
	var micrographselector = document.getElementById('micrographselector');
	
	for(var i = 0; i < micrographs.length; i++){
			var select = document.createElement("option");
			select.setAttribute("data-stackfile", micrographs[i].stackfile);
			select.setAttribute("data-particlecount", micrographs[i].particlecount);
			select.innerHTML = micrographs[i].stackfile;
			micrographselector.appendChild(select);
	}
}

function viewparticlesSelectMicrograph(caller) {
	var micrograph = caller.options[caller.selectedIndex];
	var rootdirectory = document.getElementById('rootdirectory').value;
	var snapshots = document.getElementById('snapshots');
	snapshots.innerHTML = "";
	for (var i = 0; i < micrograph.getAttribute("data-particlecount"); i++) {
		var div = document.createElement("div");
		div.className = 'snapshot';
		var img = document.createElement("img");
		img.setAttribute('data-rootdirectory', rootdirectory);
		img.setAttribute('data-stackfile', micrograph.getAttribute("data-stackfile"));
		img.setAttribute('data-frameid', i + 1);
		img.style.width = "200px";
		img.style.height = "200px";
		img.alt = "Loading ...";
		img.oncontextmenu = function(){
			var controlstarget = document.getElementById('controlstarget');
			controlstarget.value = this.className;
			showHideControlsPopup();
			return false;
		}
		div.appendChild(img);
		snapshots.appendChild(div);
	} 
	 loadImages2();
}


function loadImages(){
	var unloaded = document.querySelectorAll('[data-show="yes"]')[0];
	var brightness = document.getElementById('brightness').value;
	var contrast = document.getElementById('contrast').value;
	unloaded.removeAttribute('data-show');
	unloaded.addEventListener('load', loadImages);
	unloaded.src = "JPEGhandler?filename="+ unloaded.getAttribute('data-rootdirectory') + "/" + unloaded.getAttribute('data-stackfile') + "&contrast=" + contrast + "&brightness=" + brightness + "&frameid=" + unloaded.getAttribute('data-frameid');
}

function loadImages2(){

	var snapshots = document.getElementsByClassName('snapshot');
	var pane = document.getElementById('snapshots');
	paneBoundingClientRect = pane.getBoundingClientRect();
	var topcutoff = paneBoundingClientRect.top;
	var bottomcutoff = paneBoundingClientRect.bottom + (paneBoundingClientRect.bottom - paneBoundingClientRect.top);
	for (var i = 0; i < snapshots.length; i++){
		var snapshotBoundingClientRect = snapshots[i].getBoundingClientRect();
		var images = snapshots[i].getElementsByTagName('img');
		if(snapshotBoundingClientRect.bottom > topcutoff && snapshotBoundingClientRect.top < bottomcutoff) {
			for (var j = 0; j < images.length; j++){
				images[j].setAttribute('data-show',"yes");
			}
		} else 	{
			for (var j = 0; j < images.length; j++){
				images[j].removeAttribute('src');
				images[j].removeAttribute('data-show');
				images[j].removeEventListener("load", loadImages);
			}
		}
	}
	loadImages();
}

function updateControls () {
	var controlstarget = document.getElementById('controlstarget');
	var images = document.getElementsByClassName(controlstarget.value);
	var brightness = document.getElementById('brightness').value;
	var contrast = document.getElementById('contrast').value;
	
	for (var i = 0; i < images.length; i++){
		images[i].setAttribute('data-loaded',"no");
		
		images[i].setAttribute('data-src', "JPEGhandler?filename=" + images[i].id + "&contrast=" + contrast + "&brightness=" + brightness + "&frameid=" + images[i].getAttribute('data-frameid'));
	}
	showHideControlsPopup();
	loadImages2();
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

function showProjectHistory(){
	var mainpaneiframe = parent.document.getElementById('mainpaneiframe');
	mainpaneiframe.src = "projecthistory.html"
}

getParticleViewData();
