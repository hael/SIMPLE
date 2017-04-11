

function pipelineView(jobdir){
	var pipelineviewtable = document.getElementById('pipelineviewtable');
	pipelineviewtable.innerHTML = "";
	var pipelinepspeccontrast = document.getElementById('pipelinepspeccontrast');
	var pipelinepspecbrightness = document.getElementById('pipelinepspecbrightness');
	var pipelinethumbcontrast = document.getElementById('pipelinethumbcontrast');
	var pipelinethumbbrightness = document.getElementById('pipelinethumbbrightness');
	var pipelinepspeczoom = document.getElementById('pipelinepspeczoom');
	var pipelinethumbzoom = document.getElementById('pipelinethumbzoom');
	pipelinepspeccontrast.onchange = function(){pipelineViewPage()};
	pipelinepspecbrightness.onchange = function(){pipelineViewPage()};
	pipelinethumbcontrast.onchange = function(){pipelineViewPage()};
	pipelinethumbbrightness.onchange = function(){pipelineViewPage()};
	pipelinepspeczoom.onchange = function(){updateZoom()};
	pipelinethumbzoom.onchange = function(){updateZoom()};
	showTask("pipelineviewtask");
	var command = 'cmd=pipelineview dir=' + jobdir + ' usr=local pspeccontrast='+pipelinepspeccontrast.value+' pspecbrightness='+ pipelinepspecbrightness.value+ ' thumbcontrast='+pipelinethumbcontrast.value+' thumbbrightness='+pipelinethumbbrightness.value;
	//var command = 'cmd=pipelineview dir=' + jobdir + ' usr=local';
	var pipelineviewfolder = document.getElementById('pipelineviewfolder');
	pipelineviewfolder.value = jobdir;
	ws.send(command);
}

function pipelineViewPage(){
	var pipelineviewtable = document.getElementById('pipelineviewtable');
	pipelineviewtable.innerHTML = "";
	var pipelinepspeccontrast = document.getElementById('pipelinepspeccontrast').value;
	var pipelinepspecbrightness = document.getElementById('pipelinepspecbrightness').value;
	var pipelinethumbcontrast = document.getElementById('pipelinethumbcontrast').value;
	var pipelinethumbbrightness = document.getElementById('pipelinethumbbrightness').value;
	var pipelinepagenumber = document.getElementById('pipelinepagenumber').value;
	var pipelineviewfolder = document.getElementById('pipelineviewfolder').value;
	var pipelineviewpagecounter = document.getElementsByName('pipelinepagecount')[0];
	var pipelineviewpagecount = pipelineviewpagecounter.options[pipelineviewpagecounter.selectedIndex].value;
	var command = 'cmd=pipelineviewpage dir=' + pipelineviewfolder + ' usr=local pagecount=' + pipelineviewpagecount + ' page='+pipelinepagenumber+' pspeccontrast='+pipelinepspeccontrast+' pspecbrightness='+pipelinepspecbrightness + ' thumbcontrast='+pipelinethumbcontrast+' thumbbrightness='+pipelinethumbbrightness;
	ws.send(command);
}

function addPipelineView(data){
	var pspeccontrast = document.getElementById('pipelinepspeccontrast').value;
	var pspecbrightness = document.getElementById('pipelinepspecbrightness').value;
	var thumbcontrast = document.getElementById('pipelinethumbcontrast').value;
	var thumbbrightness = document.getElementById('pipelinethumbbrightness').value;
	var jobdir = document.getElementById('pipelineviewfolder').value;
	var pipelineviewpagecounter = document.getElementsByName('pipelinepagecount')[0];
	var pipelineviewpagecount = pipelineviewpagecounter.options[pipelineviewpagecounter.selectedIndex].value;
	var micrographs = getArgument(data, 'micrographs');
	var micrographsarray = micrographs.split(";");
	micrographsselectedarray = [];
	console.log(pipelineviewpagecounter.options[0].value);
	for(var i = 0; i <  micrographsarray.length; i++) {
		micrographsselectedarray.push([micrographsarray[i],'false']);
	}
	var command = 'cmd=pipelineviewpage dir=' + jobdir + ' usr=local pagecount=' + pipelineviewpagecount +' pspeccontrast='+pspeccontrast+' pspecbrightness='+pspecbrightness + ' thumbcontrast='+thumbcontrast+' thumbbrightness='+thumbbrightness;
	ws.send(command);
}

function addPipelineViewData(data){
	var pipelineviewtable = document.getElementById('pipelineviewtable');
	var row = pipelineviewtable.insertRow();
	
	var cell1 = row.insertCell();
	var cell2 = row.insertCell();
	var cell3 = row.insertCell();
	var cell4 = row.insertCell();
	var cell5 = row.insertCell();
	var cell6 = row.insertCell();
	var cell7 = row.insertCell();
	var cell8 = row.insertCell();
	
	cell2.className= 'micrograph';
	cell3.className= 'dfx';
	cell4.className= 'dfy';
	cell5.className= 'angast';
	cell6.className= 'pspec pipelineimg';
	cell7.className = 'thumb pipelineimg';
	cell8.className = 'boxcount';
	
	if (getArgument(data, 'selected') == 'true'){
		for(var i = 0; i <  micrographsselectedarray.length; i++) {
			if(micrographsselectedarray[i][0] == getArgument(data, 'micrograph')){
				micrographsselectedarray[i][1] = 'true';
			}
		}
		cell1.innerHTML = "<input class=selected type=checkbox checked=checked>";
		row.className = "pipelineviewselected";
	} else {
		cell1.innerHTML = "<input class=selected type=checkbox>";
	}
	
	
	if (getArgument(data, 'micrograph') != ''){
		cell2.innerHTML = getArgument(data, 'micrograph');
	}
	
	if (getArgument(data, 'dfx') != ''){
		cell3.innerHTML = getArgument(data, 'dfx');
	}
	
	if (getArgument(data, 'dfy') != ''){
		cell4.innerHTML = getArgument(data, 'dfy');
	}
		
	if (getArgument(data, 'angast') != ''){
		cell5.innerHTML = getArgument(data, 'angast');
	}
	
	var pspecdiv = document.createElement("div");
	pspecdiv.className = "pspecdiv";
	var pspecimg = document.createElement("img");
	if (getArgument(data, 'pspec') != ''){
		var pipelinepspeczoom = document.getElementById('pipelinepspeczoom');
		var pipelinepspeczoomval = pipelinepspeczoom.value;
		pspecimg.style.transform = "scale("+(1 / Number(pipelinepspeczoomval)) + "," +(1 / Number(pipelinepspeczoomval)) + ")";
		pspecimg.className = "pspecimg";
		pspecimg.setAttribute("src","data:image/png;base64," + getArgument(data, 'pspec'));
	}
	pspecdiv.appendChild(pspecimg);
	cell6.appendChild(pspecdiv);
	
	var thumbdiv = document.createElement("div");
	thumbdiv.className = "thumbdiv";
	var thumbimg = document.createElement("img");
	if (getArgument(data, 'thumb') != ''){
		var pipelinethumbzoom = document.getElementById('pipelinethumbzoom');
		var pipelinethumbzoomval = pipelinethumbzoom.value;
		thumbimg.style.transform = "scale("+(1 / Number(pipelinethumbzoomval)) + "," +(1 / Number(pipelinethumbzoomval)) + ")";
		thumbimg.className = "thumbimg"
		thumbimg.setAttribute("src","data:image/png;base64," + getArgument(data, 'thumb'));
	}
	thumbdiv.appendChild(thumbimg);
	cell7.appendChild(thumbdiv);
	
	if (getArgument(data, 'boxcount') != ''){
		cell8.innerHTML = getArgument(data, 'boxcount');
		var thumbimg = document.createElement("img");
		thumbimg.setAttribute("src","img/eye.png");
		cell8.appendChild(thumbimg);
		cell8.onclick = function(){
			pipelineViewMicrograph();
			var pipelineviewboxselectmicrograph = document.getElementById('pipelineviewboxselectmicrograph');
			pipelineviewboxselectmicrograph.innerHTML = "";
			var jobdir = document.getElementById('pipelineviewfolder').value;
			ws.send("cmd=pipelineviewmicrograph dir="+jobdir + " micrograph=" +getArgument(data, 'micrograph'));
		}
	}
	
	//row.onclick = function(){
	cell1.onclick = function(){
		if(row.className == "pipelineviewselected"){
			row.className = "";
			row.getElementsByClassName('selected')[0].checked = false;
			for(var i = 0; i <  micrographsselectedarray.length; i++) {
				if(micrographsselectedarray[i][0] == getArgument(data, 'micrograph')){
					micrographsselectedarray[i][1] = 'false';
				}
			}
		} else {
			row.className = "pipelineviewselected";
			row.getElementsByClassName('selected')[0].checked = true;
			for(var i = 0; i <  micrographsselectedarray.length; i++) {
				if(micrographsselectedarray[i][0] == getArgument(data, 'micrograph')){
					micrographsselectedarray[i][1] = 'true';
				}
			}
		}
		
	};
	var pipelinepagenumber = document.getElementById('pipelinepagenumber');
	if(getArgument(data, 'page') != '' && getArgument(data, 'page') != pipelinepagenumber.value){
		pipelinepagenumber.value = getArgument(data, 'page');
	}
	var pipelinepagemax = document.getElementById('pipelinepagemax');
	if(getArgument(data, 'pagemax') != '' && getArgument(data, 'pagemax') != pipelinepagemax.innerHTML){
		pipelinepagemax.innerHTML = getArgument(data, 'pagemax');
	}
	var numberrows = pipelineviewtable.rows.length;
	var pipelinecount = document.getElementById('pipelinecount');
	pipelinecount.innerHTML = "Display count " + numberrows;
}

function updateZoom(){
		var pipelinepspeczoom = document.getElementById('pipelinepspeczoom');
		var pipelinepspeczoomval = pipelinepspeczoom.value;
		var pipelinethumbzoom = document.getElementById('pipelinethumbzoom');
		var pipelinethumbzoomval = pipelinethumbzoom.value;
		var pspecimgs = document.getElementsByClassName('pspecimg');
		for(var i = 0; i <  pspecimgs.length; i++) {
			pspecimgs[i].style.transform = "scale("+(1 / Number(pipelinepspeczoomval)) + "," +(1 / Number(pipelinepspeczoomval)) + ")";
		}
		var thumbimgs = document.getElementsByClassName('thumbimg');
		for(var i = 0; i <  thumbimgs.length; i++) {
			thumbimgs[i].style.transform = "scale("+(1 / Number(pipelinethumbzoomval)) + "," +(1 / Number(pipelinethumbzoomval)) + ")";
		}
}

function pipelineViewPageUp(){
	var pipelinepagenumber = document.getElementById('pipelinepagenumber');
	var currentpage = pipelinepagenumber.value;
	pipelinepagenumber.value = Number(currentpage) + 1;
	pipelineViewPage();
}

function pipelineViewPageDown(){
	var pipelinepagenumber = document.getElementById('pipelinepagenumber');
	var currentpage = pipelinepagenumber.value;
	pipelinepagenumber.value = Number(currentpage) - 1;
	pipelineViewPage();
}

function pipelineViewSelectAll(){
	var pipelineviewtable = document.getElementById('pipelineviewtable');
	var rows = pipelineviewtable.getElementsByTagName('tr');
	for (var i = 0; i < rows.length; i++) {
		rows[i].className = "pipelineviewselected";
		rows[i].getElementsByClassName('selected')[0].checked = true;
	}
}

function pipelineViewClearAll(){
	var pipelineviewtable = document.getElementById('pipelineviewtable');
	var rows = pipelineviewtable.getElementsByTagName('tr');
	for (var i = 0; i < rows.length; i++) {
		rows[i].className = "";
		rows[i].getElementsByClassName('selected')[0].checked = false;
	}
}

function pipelineViewSave(){
	document.getElementById('pipelineviewsave').style.display = "block";
}

function cancelPipelineViewSave(){
	document.getElementById('pipelineviewsave').style.display = "none";
}

function confirmPipelineViewSave(){
	var dir = document.getElementById('pipelineviewfolder').value;
	var filename = document.getElementById('pipelineviewsavefilename').value;
	var outdir = document.getElementById('pipelineviewsavefolder').value;
	var command = 'cmd=pipelinesaveselection dir=' + dir + ' outdir=' + outdir + ' filename='+ filename + ' selections=';
	for(var i = 0; i <  micrographsselectedarray.length; i++) {
		if(micrographsselectedarray[i][1] == 'true'){
			command += micrographsselectedarray[i][0];
			command += ";";
		}
	}
	console.log(command);
	ws.send(command);
	document.getElementById('pipelineviewsave').style.display = "none";
}

function pipelineViewMicrograph(){
	document.getElementById('pipelineviewboxselect').style.display = "block";
	var pipelineviewboxselectmicrograph = document.getElementById('pipelineviewboxselectmicrograph');
	pipelineviewboxselectmicrograph.innerHTML = "";
}

function addPipelineViewMicrograph(data){
	var pipelineviewboxselectmicrograph = document.getElementById('pipelineviewboxselectmicrograph');
	pipelineviewboxselectmicrograph.innerHTML = "";
	var micrographimg = document.createElement("img");
	micrographimg.setAttribute("src","data:image/png;base64," + getArgument(data, 'micrographdata'));
	pipelineviewboxselectmicrograph.appendChild(micrographimg);
	var micrographcanvas = document.createElement("canvas");
	micrographcanvas.id = "micrographcanvas";
	micrographcanvas.width = "800";
	micrographcanvas.height = "800";
	micrographcanvas.addEventListener('click', function(e) {
        pipelineViewdrawBoxes(e);
	});
	pipelineviewboxselectmicrograph.appendChild(micrographcanvas);
	var canvas = micrographcanvas.getContext("2d");
	xscale = micrographimg.width/getArgument(data, 'nx');
	yscale = micrographimg.height/getArgument(data, 'ny');
	var boxes = getArgument(data, 'boxes');
	boxesarray = boxes.split(";");
	pipelineViewdrawBoxes();
}

function pipelineViewdrawBoxes(e){
	var micrographcanvas = document.getElementById("micrographcanvas");
	var canvas = micrographcanvas.getContext("2d");
	canvas.clearRect(0,0,800,800);
	for(box in boxesarray){
		canvas.beginPath();
		var boxarray = boxesarray[box].split(",");
		var xcoord = boxarray[0] * xscale;
		var ycoord = boxarray[1] * yscale;
		var xboxsize = boxarray[2] * xscale
		var yboxsize = boxarray[2] * yscale
		var xmax = xcoord + xboxsize;
		var ymax = ycoord + yboxsize;
		if(typeof e !== 'undefined'){
			if(e.offsetX > xcoord && e.offsetX < xmax && e.offsetY > ycoord && e.offsetY < ymax){
				delete boxesarray[box];
			} else {
				canvas.rect(xcoord,ycoord,xboxsize,yboxsize);
			}
		} else {
			canvas.rect(xcoord,ycoord,xboxsize,yboxsize);
		}
		canvas.strokeStyle = "green";
		canvas.stroke();
	}
}

function pipelineViewClearBoxes(){
	boxesarray = [];
	pipelineViewdrawBoxes();
}

function pipelineViewSaveBoxes(){
	for(box in boxesarray){
		var boxarray = boxesarray[box].split(",");
		console.log(boxarray[0]);
	}
	pipelineViewClose();
}

function pipelineViewClose(){
	document.getElementById('pipelineviewboxselect').style.display = "none";
}
