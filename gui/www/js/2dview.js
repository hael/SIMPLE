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

function get2DViewData () {
 getAjax('JSONhandler?function=2dview&' + window.location.toString().split('?')[1], function(data){view2DInit(data)});
}

function view2DInit (data) {
	var JSONdata = JSON.parse(data);
	
	var iterations = JSONdata.iterations;
	iterations.sort(function(a,b) {return (a.cavgs > b.cavgs) ? 1 : ((b.cavgs > a.cavgs) ? -1 : 0);} );
	var rootdirectory = document.getElementById('rootdirectory');
	rootdirectory.value = JSONdata.rootdirectory;
	var saveselectionfolder = document.getElementById('saveselectionfolder');
	saveselectionfolder.value = JSONdata.rootdirectory;
	var saveselectionparticlesfolder = document.getElementById('saveselectionparticlesfolder');
	saveselectionparticlesfolder.value = JSONdata.rootdirectory;
	
	var iterationselector = document.getElementById('iterationselector');
	
	for(var i = 0; i < iterations.length; i++){
		if(! iterations[i].cavgs.includes("ranked")){
			var select = document.createElement("option");
			select.setAttribute("data-classdoc", iterations[i].classdoc);
			select.setAttribute("data-cavgs", iterations[i].cavgs);
			select.setAttribute("data-prime2ddoc", iterations[i].prime2ddoc);
			select.innerHTML = iterations[i].classdoc.replace(".txt", "").replace("classdoc_","");
			iterationselector.appendChild(select);
		}
	}
}

function view2DSelectIteration (caller) {
	var iteration = caller.options[caller.selectedIndex];
	var rootdirectory = document.getElementById('rootdirectory').value;
	getAjax('JSONhandler?function=2dviewiteration&classdoc=' + rootdirectory + "/" + iteration.getAttribute("data-classdoc") + "&cavgs="  + rootdirectory + "/" + iteration.getAttribute("data-cavgs")  + "&prime2ddoc="  + rootdirectory + "/" + iteration.getAttribute("data-prime2ddoc"), function(data){showFileViewerData(data)});
}


function showFileViewerData (data){
	var JSONdata = JSON.parse(data);
	
	//var rootdirectory = document.getElementById('rootdirectory');
	//rootdirectory.value = JSONdata.rootdirectory;
	
	var sorton = document.getElementById('sorton');
	sorton.innerHTML = "";
	
	var sortoptions = ["class", "pop", "res","corr", "w"]; 
	
	for (var i = 0; i < sortoptions.length; i++) {
		var option = document.createElement("option");
		option.value = sortoptions[i];
		option.innerHTML = sortoptions[i];
		sorton.appendChild(option);
	}
	
	var selectionattribute = document.getElementById('selectionattribute');
	selectionattribute.innerHTML = "";
	for (var i = 0; i < sortoptions.length; i++) {
		var option = document.createElement("option");
		option.value = sortoptions[i];
		option.innerHTML = sortoptions[i];
		selectionattribute.appendChild(option);
	}
	
	var snapshots = document.getElementById('snapshots');
	snapshots.innerHTML = "";
	for (var i = 0; i < JSONdata.snapshots.length; i++) {
		var div = document.createElement("div");
		div.className = 'snapshot';
		div.setAttribute('data-selected', "yes");
		div.setAttribute('data-rootdirectory', JSONdata.rootdirectory);
		div.onclick = function() {
			selectUnselectSnapshot(this);
		}
		var attributesdiv = document.createElement("div");
		attributesdiv.className = 'attributes';	
			
		var keys = Object.keys(JSONdata.snapshots[0]);
		for(var j = 0; j < keys.length; j++) {
			div.setAttribute('data-'+keys[j], JSONdata.snapshots[i][keys[j]]);
			var attributediv = document.createElement("div");
			attributediv.innerHTML = keys[j];
			attributediv.innerHTML += " : ";
			attributediv.innerHTML += JSONdata.snapshots[i][keys[j]];
			attributesdiv.appendChild(attributediv);
		}
		div.setAttribute('data-prime2ddoc', JSONdata.prime2ddoc);
		div.setAttribute('data-id', JSONdata.snapshots[i].id);
		if(typeof(JSONdata.snapshots[i].frameid) != undefined){
			div.setAttribute('data-frameid', JSONdata.snapshots[i].frameid);
		}else{
			div.setAttribute('data-frameid', 0);
		}
		
			
		div.appendChild(attributesdiv);
		var gauze = document.createElement("div");
		gauze.className = "snapshotgauze";
		div.appendChild(gauze);
		snapshots.appendChild(div);
	} 

	updateFileViewer();
}

function updateFileViewer () {
	var snapshots = document.getElementsByClassName('snapshot');
	var selecteddisplayoptions = document.getElementsByClassName('selecteddisplayoptions');
	for (var i = 0; i < snapshots.length; i++){
		var images = snapshots[i].getElementsByTagName('img');
		for (var j = 0; j < images.length; j++){
			images[j].parentNode.removeChild(images[j]);
		}
		var fileimg = document.createElement("img");
		fileimg.src = "img/info.png";
		fileimg.className = "snapshotattributesbutton";
		fileimg.onclick = function(e){
			if (!e) var e = window.event;
			e.cancelBubble = true;
			if (e.stopPropagation) e.stopPropagation();
			var attributes = this.parentNode.getElementsByClassName('attributes')[0];
			if(attributes.style.visibility =="visible"){
				attributes.style.visibility = "hidden"
			}else{
				attributes.style.visibility = "visible";
			}
		};
		snapshots[i].appendChild(fileimg);
		var eyeimg = document.createElement("img");
		eyeimg.src = "img/ptcls.png";
		eyeimg.className = "snapshotviewparticlesbutton";
		eyeimg.onclick = function(e){
			if (!e) var e = window.event;
			e.cancelBubble = true;
			if (e.stopPropagation) e.stopPropagation();
			showHideParticleViewPopup();
			getAjax('JSONhandler?function=2dviewparticles&prime2ddoc=' + this.parentNode.getAttribute('data-prime2ddoc') + "&class="  + this.parentNode.getAttribute('data-class') , function(data){showParticleData(data)});
		};
		snapshots[i].appendChild(eyeimg);
		var img = document.createElement("img");
		img.setAttribute('data-src',"JPEGhandler?filename=" + snapshots[i].getAttribute('data-cavgs') + "&contrast=2&brightness=128&size=200&frameid=" + snapshots[i].getAttribute('data-frameid'));
		img.style.width = "200px";
		img.style.height = "200px";
		img.className = 'data-splreferenceimage';
		img.alt = "Loading ...";
		img.setAttribute('data-frameid', snapshots[i].getAttribute('data-frameid'));
		img.id = snapshots[i].getAttribute('data-cavgs');
		img.oncontextmenu = function(){
			var controlstarget = document.getElementById('controlstarget');
			controlstarget.value = this.className;
			showHideControlsPopup();
			return false;
		}
		snapshots[i].appendChild(img);
	}

	document.getElementById('snapshots').addEventListener('scroll', function(){loadImages2()});
	loadImages2();
}

function loadImages(){
	//var unloaded = document.querySelectorAll('[data-loaded="no"]')[0];
	var unloaded = document.querySelectorAll('[data-show="yes"]')[0];
	unloaded.removeAttribute('data-show');
	unloaded.addEventListener('load', loadImages)
	unloaded.src = unloaded.getAttribute('data-src');
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
			for (var j = 2; j < images.length; j++){
				images[j].setAttribute('data-show',"yes");
			}
		} else 	{
			for (var j = 2; j < images.length; j++){
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

function showHideMenu () {
	var headermenu = document.getElementById('headermenu');
	if (headermenu.style.display == "block") {
		headermenu.style.display = "none";
	} else {
		headermenu.style.display = "block";
	}
}


function showHideSortPopup () {
	var sortpopup = document.getElementById('sortpopup');
	var gauze = document.getElementById('gauze');
	if (sortpopup.style.display == "block") {
		sortpopup.style.display = "none";
		gauze.style.display = "none";
	} else {
		sortpopup.style.display = "block";
		gauze.style.display = "block";
	}
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

function showHideSelectFilePopup () {
	var selectfilepopup = document.getElementById('selectfilepopup');
	var gauze = document.getElementById('gauze');
	if (selectfilepopup.style.display == "none") {
		selectfilepopup.style.display = "block";
		gauze.style.display = "block";
	} else {
		selectfilepopup.style.display = "none";
		gauze.style.display = "none";
	}
}

function showHideSaveSelectionPopup () {
	var saveselectionpopup = document.getElementById('saveselectionpopup');
	var gauze = document.getElementById('gauze');
	document.getElementById('saveselectionbutton').innerHTML = "Save";
	if (saveselectionpopup.style.display == "block") {
		saveselectionpopup.style.display = "none";
		gauze.style.display = "none";
	} else {
		saveselectionpopup.style.display = "block";
		gauze.style.display = "block";
	}
}

function showHideSelectionPopup () {
	var selectionpopup = document.getElementById('selectionpopup');
	var gauze = document.getElementById('gauze');
	if (selectionpopup.style.display == "block") {
		selectionpopup.style.display = "none";
		gauze.style.display = "none";
	} else {
		selectionpopup.style.display = "block";
		gauze.style.display = "block";
	}
}

function showHideSaveSelectionParticlesPopup () {
	var saveselectionpopup = document.getElementById('saveselectionparticlespopup');
	var gauze = document.getElementById('gauze');
	document.getElementById('saveselectionparticlesbutton').innerHTML = "Save";
	if (saveselectionpopup.style.display == "block") {
		saveselectionpopup.style.display = "none";
		gauze.style.display = "none";
	} else {
		saveselectionpopup.style.display = "block";
		gauze.style.display = "block";
	}
}



function showHideParticleViewPopup () {
	var particleviewpopup = document.getElementById('particleviewpopup');
	var gauze = document.getElementById('gauze');
	if (particleviewpopup.style.display == "block") {
		particleviewpopup.style.display = "none";
		gauze.style.display = "none";
	} else {
		particleviewpopup.style.display = "block";
		gauze.style.display = "block";
	}
}

function showHideHeaderMenu (element) {
	var headermenu = element.parentNode.parentNode.getElementsByClassName('headermenu')[0];
	if (headermenu.style.display == "block") {
		headermenu.style.display = "none";
	} else {
		headermenu.style.display = "block";
	}
}


function selectUnselectSnapshot(element) {
	var snapshotgauze = element.getElementsByClassName('snapshotgauze')[0];
	if (snapshotgauze.style.display == "block"){
		snapshotgauze.style.display = "none";
		element.setAttribute('data-selected', "yes");
	} else {
		snapshotgauze.style.display = "block";
		element.setAttribute('data-selected', "no");
	}
}

function showFileBrowserData(data){
	var JSONdata = JSON.parse(data);
	var directories = JSONdata.directories;
	var files = JSONdata.files;
	var selectfiledirectory = document.getElementById('selectfiledirectory');
	var selectfiletable = document.getElementById('selectfiletable');
	selectfiletable.innerHTML = "";
	directories.sort();
	files.sort();
	for (var i = 0; i < directories.length; i++) {
		var row = selectfiletable.insertRow(-1);
		row.id = JSONdata.directory + "/" + directories[i];
		var cell1 = row.insertCell(0);
		cell1.innerHTML = "<img src=img/folder.png class=folderimage>";
		var cell2 = row.insertCell(1);
		cell2.innerHTML = directories[i];
		//row.ondblclick = function(){getFileBrowserData(this.id)};
		row.onclick = function(){getFileBrowserData(this.id)};
	}
	for (var i = 0; i < files.length; i++) {
		if(files[i].includes(".mrc") || files[i].includes(".star")) {
			var row = selectfiletable.insertRow(-1);
			row.id = JSONdata.directory + "/" + files[i];
			var cell1 = row.insertCell(0);
			var cell2 = row.insertCell(1);
			cell2.innerHTML = files[i];
			row.onclick = function(){
				document.getElementById('fileviewerfilename').value = this.id
				if(typeof document.getElementsByClassName('fileviewerfileselected')[0] !== "undefined") {
					document.getElementsByClassName('fileviewerfileselected')[0].style.background = "white";
					document.getElementsByClassName('fileviewerfileselected')[0].className = "";
				}
				this.className = "fileviewerfileselected";
				this.style.background = "#6698ab";
			};
		}
	}
	selectfiledirectory.value = JSONdata.directory;
}

function getFileBrowserData(directory, filter){
	var url = 'JSONhandler?function=listdir';
	if (directory){	
		url += "&directoryname=" + directory;
	}
	if (filter){	
		url += "&filefilter=" + filter;
	}
	getAjax(url, function(data){showFileBrowserData(data)});	
}

function showProjectHistory(){
	var mainpaneiframe = parent.document.getElementById('mainpaneiframe');
	mainpaneiframe.src = "projecthistory.html"
}

function sortFileViewer(){
	showHideSortPopup();
	var sorton = document.getElementById('sorton');
	var sortoption = sorton.options[sorton.selectedIndex].value;
	var snapshot = document.getElementsByClassName('snapshot');
	snapshot = Array.prototype.slice.call(snapshot, 0);
	var sortorder = document.getElementById('sortorder');
	if (sortorder.options[sortorder.selectedIndex].value == "ascending"){
		snapshot.sort(function(a,b){
			if(Number(a.getAttribute("data-" +sortoption.toLowerCase())) < Number(b.getAttribute("data-" +sortoption.toLowerCase()))) return -1;
			if(Number(a.getAttribute("data-" +sortoption.toLowerCase())) > Number(b.getAttribute("data-" +sortoption.toLowerCase()))) return 1;
			return 0;
		});
	} else {
		snapshot.sort(function(a,b){
			if(Number(a.getAttribute("data-" +sortoption.toLowerCase())) < Number(b.getAttribute("data-" +sortoption.toLowerCase()))) return 1;
			if(Number(a.getAttribute("data-" +sortoption.toLowerCase())) > Number(b.getAttribute("data-" +sortoption.toLowerCase()))) return -1;
			return 0;
		});
	}
	var snapshots = document.getElementById('snapshots');
	snapshots.innerHTML = "";
	for (var i = 0; i < snapshot.length; i++){
		snapshots.appendChild(snapshot[i]);
	}
}

function saveSelection() {
	var saveselectionfilename = document.getElementById('saveselectionfolder').value + "/" + document.getElementById('saveselectionfilename').value;
	var selected = document.querySelectorAll('[data-selected="yes"]');
	var rootdirectory = document.getElementById('rootdirectory').value;
	var inputfilename = document.getElementsByClassName('snapshot')[0].getAttribute('data-cavgs');
	var url = 'JSONhandler?function=2dsaveselection';
	if(inputfilename){
		url += "&inputfilename=" + inputfilename;
	}
	if(saveselectionfilename){
		url += "&outputfilename=" + saveselectionfilename;
	}
	url += "&selection=";
	for (var i = 0; i < selected.length; i++) {
		url += selected[i].getAttribute("data-frameid") + ",";
	} 
	document.getElementById('saveselectionbutton').innerHTML = "Saving ...";
	getAjax(url, function(data){showHideSaveSelectionPopup()});
}

function saveSelectionParticles() {
	var saveselectionfilename = document.getElementById('saveselectionparticlesfolder').value + "/" + document.getElementById('saveselectionparticlesfilename').value;
	var inverseselected = document.querySelectorAll('[data-selected="no"]');
	var inputfilename = document.getElementsByClassName('snapshot')[0].getAttribute('data-prime2ddoc');
	var url = 'JSONhandler?function=2dsaveselectionparticles';
	if(inputfilename){
		url += "&inputfilename=" + inputfilename;
	}
	if(saveselectionfilename){
		url += "&outputfilename=" + saveselectionfilename;
	}
	url += "&inverseselection=";
	for (var i = 0; i < inverseselected.length; i++) {
		url += inverseselected[i].getAttribute("data-class") + ",";
	} 
	document.getElementById('saveselectionparticlesbutton').innerHTML = "Saving ...";
	getAjax(url, function(data){showHideSaveSelectionParticlesPopup()});
}

function applySelection() {
	var selectionattribute = document.getElementById('selectionattribute');
	var selectioncriteria = document.getElementById('selectioncriteria');
	var selectionvalue = document.getElementById('selectionvalue').value;
	var selectionappend = document.getElementById('selectionappend');
	var selectionattributevalue = selectionattribute.options[selectionattribute.selectedIndex].value;
	var selectioncriteriavalue = selectioncriteria.options[selectioncriteria.selectedIndex].value;
	var snapshots = document.getElementsByClassName('snapshot');
	if(! selectionappend.checked){
		for (var i = 0; i < snapshots.length; i++){
			snapshots[i].setAttribute('data-selected', "yes");
			var snapshotgauze = snapshots[i].getElementsByClassName('snapshotgauze')[0];
			snapshotgauze.style.display = "none";
		}
	}
	for (var i = 0; i < snapshots.length; i++){
		if (selectioncriteriavalue == "greaterthan"){
			if(Number(snapshots[i].getAttribute("data-" +selectionattributevalue.toLowerCase())) < Number(selectionvalue)){
				var snapshotgauze = snapshots[i].getElementsByClassName('snapshotgauze')[0];
				snapshotgauze.style.display = "block";
				snapshots[i].setAttribute('data-selected', "no");
			}
		}else if (selectioncriteriavalue == "lessthan"){
			if(Number(snapshots[i].getAttribute("data-" +selectionattributevalue.toLowerCase())) > Number(selectionvalue)){
				var snapshotgauze = snapshots[i].getElementsByClassName('snapshotgauze')[0];
				snapshotgauze.style.display = "block";
				snapshots[i].setAttribute('data-selected', "no");
			}	
		}
	}
	showHideSelectionPopup();
}

function showParticleData (data){
	var JSONdata = JSON.parse(data);
	var particles = document.getElementById('particles');
	var rootdirectory = document.getElementById('rootdirectory').value;
	particles.innerHTML = "";
	for (var i = 0; i < JSONdata.snapshots.length; i++) {
		var div = document.createElement("div");
		div.setAttribute('data-frameid', JSONdata.snapshots[i].frameid);
		div.setAttribute('data-src', rootdirectory + "/" + JSONdata.snapshots[i].stackfile);
		div.className = 'particle';
		particles.appendChild(div);
	} 
	updateParticleViewer();
}

function updateParticleViewer () {
	var particles = document.getElementsByClassName('particle');
	var rootdirectory = document.getElementById('rootdirectory').value;
	var particlecontrast = document.getElementById('particlecontrast').value;
	var particlebrightness= document.getElementById('particlebrightness').value;
	var particlesize= document.getElementById('particlesize').value;
	
	for (var i = 0; i < particles.length; i++){
		var img = document.createElement("img");
		img.setAttribute('data-src',"JPEGhandler?filename=" + particles[i].getAttribute('data-src') + "&contrast=" +particlecontrast+"&brightness="+particlebrightness+"&size="+particlesize+"&frameid=" + particles[i].getAttribute('data-frameid'));
		img.style.width = particlesize+"px";
		img.style.height = particlesize+"px";
		img.alt = "Loading ...";
		img.setAttribute('data-frameid', particles[i].getAttribute('data-frameid'));
		particles[i].appendChild(img);
	}
	document.getElementById('particles').addEventListener('scroll', function(){loadParticles2()});
	loadParticles2();
}

function loadParticles(){
	//var unloaded = document.querySelectorAll('[data-loaded="no"]')[0];
	var unloaded = document.querySelectorAll('[data-particleshow="yes"]')[0];
	unloaded.removeAttribute('data-particleshow');
	unloaded.addEventListener('load', loadParticles)
	unloaded.src = unloaded.getAttribute('data-src');
}

function loadParticles2(){

	var particles= document.getElementsByClassName('particle');
	var pane = document.getElementById('particles');
	paneBoundingClientRect = pane.getBoundingClientRect();
	var topcutoff = paneBoundingClientRect.top;
	var bottomcutoff = paneBoundingClientRect.bottom + (paneBoundingClientRect.bottom - paneBoundingClientRect.top);
	for (var i = 0; i < particles.length; i++){
		var snapshotBoundingClientRect = particles[i].getBoundingClientRect();
		var images = particles[i].getElementsByTagName('img');
		if(snapshotBoundingClientRect.bottom > topcutoff && snapshotBoundingClientRect.top < bottomcutoff) {
			images[0].setAttribute('data-particleshow',"yes");
		} else 	{
			images[0].removeAttribute('src');
			images[0].removeAttribute('data-particleshow');
			images[0].removeEventListener("load", loadImages);
		}
	}
	loadParticles();
}

function updateParticlesControls () {
	var particles= document.getElementsByClassName('particle');
	var rootdirectory = document.getElementById('rootdirectory').value;
	var particlecontrast = document.getElementById('particlecontrast').value;
	var particlebrightness= document.getElementById('particlebrightness').value;
	var particlesize= document.getElementById('particlesize').value;
	
	for (var i = 0; i < particles.length; i++){
		var img = particles[i].getElementsByTagName('img')[0];
		img.setAttribute('data-loaded',"no");
		img.setAttribute('data-src',"JPEGhandler?filename=" + particles[i].getAttribute('data-src') + "&contrast=" +particlecontrast+"&brightness="+particlebrightness+"&size="+particlesize+"&frameid=" + particles[i].getAttribute('data-frameid'));
		img.style.width = particlesize+"px";
		img.style.height = particlesize+"px"
	}
	loadParticles2();
}


get2DViewData();
