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

function getFileviewerData () {
 var fileviewerfilename = document.getElementById('fileviewerfilename').value;
 getAjax('JSONhandler?function=viewfile&filename='+fileviewerfilename, function(data){showFileViewerData(data)});
}
 

function showFileViewerData (data){
	var JSONdata = JSON.parse(data);
	
	var rootdirectory = document.getElementById('rootdirectory');
	rootdirectory.value = JSONdata.rootdirectory;
	
	var displaytable = document.getElementById('displaytable');
	displaytable.innerHTML = "";
	for (var i = 0; i < JSONdata.displayoptions.length; i++) {
		var row = displaytable.insertRow(-1);
		var cell1 = row.insertCell(0);
		var cell2 = row.insertCell(1);
		var checkbox = document.createElement("input");
		checkbox.setAttribute('type', 'checkbox');
		checkbox.id = JSONdata.displayoptions[i];
		checkbox.onchange = function(){
			if(this.checked){
				this.className = "selecteddisplayoptions";
			} else {
				this.className = "";
			}
		}
		cell1.appendChild(checkbox);
		cell2.innerHTML = JSONdata.displayoptions[i];
	}
	
	var sorttable = document.getElementById('sorttable');
	sorttable.innerHTML = "";
	for (var i = 0; i < JSONdata.sortoptions.length; i++) {
		var row = sorttable.insertRow(-1);
		var cell1 = row.insertCell(0);
		var cell2 = row.insertCell(1);
		var radio = document.createElement("input");
		radio.setAttribute('type', 'radio');
		radio.name = "sortvalue";
		radio.value = JSONdata.sortoptions[i];
		cell1.appendChild(radio);
		cell2.innerHTML = JSONdata.sortoptions[i];
	}
	
	var selectionattribute = document.getElementById('selectionattribute');
	selectionattribute.innerHTML = "";
	for (var i = 0; i < JSONdata.sortoptions.length; i++) {
		var option = document.createElement("option");
		option.value = JSONdata.sortoptions[i];
		option.innerHTML = JSONdata.sortoptions[i];
		selectionattribute.appendChild(option);
	}
	
	var snapshots = document.getElementById('snapshots');
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
		for(var j = 0; j < JSONdata.sortoptions.length; j++) {
			div.setAttribute('data-'+JSONdata.sortoptions[j], JSONdata.snapshots[i][JSONdata.sortoptions[j]]);
			var attributediv = document.createElement("div");
			attributediv.innerHTML = JSONdata.sortoptions[j];
			attributediv.innerHTML += " : ";
			attributediv.innerHTML += JSONdata.snapshots[i][JSONdata.sortoptions[j]];
			attributesdiv.appendChild(attributediv);
		}
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
	showHideSelectFilePopup();
	showHideDisplayPopup ();
	updateFileViewer();
}

function updateFileViewer () {
	var snapshots = document.getElementsByClassName('snapshot');
	var selecteddisplayoptions = document.getElementsByClassName('selecteddisplayoptions');
	for (var i = 0; i < snapshots.length; i++){
		//snapshots[i].innerHTML = "";
		//var gauze = document.createElement("div");
		//gauze.className = "snapshotgauze";
		//snapshots[i].appendChild(gauze);
		var images = snapshots[i].getElementsByTagName('img');
		for (var j = 0; j < images.length; j++){
			images[j].parentNode.removeChild(images[j]);
		}
		var fileimg = document.createElement("img");
		fileimg.src = "img/file.png";
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
		for (var j = 0; j < selecteddisplayoptions.length; j++){
			var img = document.createElement("img");
			
			//if(snapshots[i].getAttribute('data-stackframe') != null){
				//img.setAttribute('data-src',"JPEGhandler?filename=" + snapshots[i].getAttribute('data-rootdirectory') + "/" + snapshots[i].getAttribute('data-frane') + "&contrast=2&brightness=128&size=200&nz=" + snapshots[i].getAttribute('data-id'));
			//} else {
			img.setAttribute('data-src',"JPEGhandler?filename=" + snapshots[i].getAttribute('data-rootdirectory') + "/" + snapshots[i].getAttribute('data-'+selecteddisplayoptions[j].id) + "&contrast=2&brightness=128&size=200&frameid=" + snapshots[i].getAttribute('data-frameid'));
			//}
			
		//	img.setAttribute('data-loaded',"no");
			img.style.width = "200px";
			img.style.height = "200px";
			img.className = 'data-'+selecteddisplayoptions[j].id;
			img.alt = "Loading ...";
			img.setAttribute('data-frameid', snapshots[i].getAttribute('data-frameid'));
			img.id = snapshots[i].getAttribute('data-rootdirectory') + "/" + snapshots[i].getAttribute('data-' + selecteddisplayoptions[j].id);
			img.oncontextmenu = function(){
				var controlstarget = document.getElementById('controlstarget');
				controlstarget.value = this.className;
				showHideControlsPopup();
				return false;
			}
			snapshots[i].appendChild(img);
		}
	}
	//loadImages();
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
			for (var j = 1; j < images.length; j++){
				images[j].setAttribute('data-show',"yes");
			}
		} else 	{
			for (var j = 1; j < images.length; j++){
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

function showHideDisplayPopup () {
	var displaypopup = document.getElementById('displaypopup');
	var gauze = document.getElementById('gauze');
	if (displaypopup.style.display == "block") {
		displaypopup.style.display = "none";
		gauze.style.display = "none";
	} else {
		displaypopup.style.display = "block";
		gauze.style.display = "block";
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

function sortFileViewer(){
	var sortvalues = document.getElementsByName('sortvalue');
	for(var i = 0; i < sortvalues.length; i++){
		if(sortvalues[i].checked){
			var sortoption = sortvalues[i].value;
			break;
		}
	}
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
	var saveselectionfilename = document.getElementById('saveselectionfilename').value;
	var inverseselected = document.querySelectorAll('[data-selected="no"]');
	var fileviewerfilename = document.getElementById('fileviewerfilename').value;
	var rootdirectory = document.getElementById('rootdirectory').value;
	var url = 'JSONhandler?function=savesel';
	if(fileviewerfilename){
		url += "&inputfilename=" + fileviewerfilename;
	}
	if(saveselectionfilename){
		url += "&outputfilename=" + rootdirectory + "/" + saveselectionfilename;
	}
	url += "&inverseselection=";
	for (var i = 0; i < inverseselected.length; i++) {
		url += inverseselected[i].getAttribute("data-id") + ",";
	} 
	getAjax(url, function(data){showHideSaveSelectionPopup()});
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

getFileBrowserData();
