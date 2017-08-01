class Websocket{
	
	constructor() {
		var self = this;
		this.socket = new WebSocket('ws://localhost:8080/');
		this.socket.onopen = function(){self.socket.close()};
	}
	
	setupSocket(command, handler){
		var self = this;
		var state = this.socket.readyState;
		if(state == 3){
			this.socket = new WebSocket('ws://localhost:8080/');
			this.socket.onopen = function(){self.socket.send(command)};
			this.socket.onclose = function(){console.log("close")};
			this.socket.onmessage = function(event){window[handler](event.data)};
			clearInterval(this.status);
		}else{
			this.socket.close();
		}
	}
	
	execute(command, handler){
		this.status = setInterval(this.setupSocket.bind(this), 100, command, handler);
	}
}

var projectSelectorWS = new Websocket();
var projectHistoryWS = new Websocket();
var newProjectWS = new Websocket();
var selectFolderWS = new Websocket();
var selectFileWS = new Websocket();
var fileBrowserWS = new Websocket();
var pipelineViewWS = new Websocket();
var runTaskWS = new Websocket();
var filesViewWS = new Websocket();

function getArgument(argumentstring, argumentname){
	var argumentarray = argumentstring.split(' ');
	var returnvalue = "";
	for (var i = 0; i < argumentarray.length; i++) {
		var keyvalue = argumentarray[i].split('=');
		if (keyvalue[0] == argumentname){
			returnvalue = keyvalue[1].replace(/\^/g, " ").replace(/\|/g, "=");
		}
	}
	return returnvalue;
}

function includeHTML() {
  var z, i, elmnt, file, xhttp;
  z = document.getElementsByTagName("*");
  for (i = 0; i < z.length; i++) {
    elmnt = z[i];
    file = elmnt.getAttribute("include-html");
    if (file) {
      xhttp = new XMLHttpRequest();
      xhttp.onreadystatechange = function() {
        if (this.readyState == 4 && this.status == 200) {
          elmnt.innerHTML = this.responseText;
          elmnt.removeAttribute("include-html");
          includeHTML();
        }
      }
      xhttp.open("GET", file, true);
      xhttp.send();
      return;
    }
  }
}

function projectSelector(){
	projectSelectorWS.execute("cmd=projectlist", "addProject");
}

function selectProject(){
	var projectselector = document.getElementById('projectselector');
	var project = projectselector.options[projectselector.selectedIndex].value;
	if (project == 'createnew'){
		var newprojectpopup = document.getElementById('newprojectpopup');
		showGauze();
		newprojectpopup.style.display = 'block';
	} else {
		document.getElementById('historytable').innerHTML = "";
		projectHistoryWS.execute("cmd=jobslist projectname="+project, "addJob");
		projectHistoryFull();
	}
}

function refreshProject(){
	var projectselector = document.getElementById('projectselector');
	var project = projectselector.options[projectselector.selectedIndex].value;
	if (project != 'createnew' && project != ''){
		document.getElementById('historytable').innerHTML = "";
		projectHistoryWS.execute("cmd=jobslist projectname="+project, "addJob");
	}	
}

function addProject(data){
	var projectselector = document.getElementById('projectselector');
	var projectname = getArgument(data, 'projectname');
	if(projectname !=""){
		var option = document.createElement("option");
		option.text = getArgument(data, 'projectname');
		option.value = getArgument(data, 'projectname');
		projectselector.add(option);
	}
}

function addJob(data){
	var historytable = document.getElementById('historytable');
	var row = historytable.insertRow();
	
	if(getArgument(data, 'jobid') != ""){
		var cell1 = row.insertCell();
		cell1.className = "historyjobid";
		var cell2 = row.insertCell();
		cell2.className = "historyjobstatus";
		var cell3 = row.insertCell();
		cell3.className = "historyjobtype";
		var cell4 = row.insertCell();
		cell4.className = "historyjobname";
		var cell5 = row.insertCell();
		cell5.className = "historyjobdescription";
		var cell6 = row.insertCell();
		cell6.className = "historyjobactions";
		
		cell1.innerHTML = getArgument(data, 'jobid') ;
		cell2.innerHTML = getArgument(data, 'jobstatus');
		cell3.innerHTML = getArgument(data, 'jobtype');
		cell4.innerHTML = getArgument(data, 'jobname');
		cell5.innerHTML = getArgument(data, 'jobdescription');
		var actionmenudiv = document.createElement('div');
		
		actionmenudiv.className = "bodybutton";
		actionmenudiv.innerHTML += "<img src=img/dropdown.png></img>";
	
		var hoverdiv = document.createElement('div');
		hoverdiv.className = "hoverdiv";
		
		if(getArgument(data, 'jobtype') == "pipeline" && getArgument(data, 'jobstatus') == "Complete"){
			hoverdiv.innerHTML += "<div onclick=pipelineView('"+getArgument(data, 'jobdir')+"')>View Output From Job</div>";
		} else if (getArgument(data, 'jobtype') == "unblur_ctffind" && getArgument(data, 'jobstatus') == "Complete"){
			hoverdiv.innerHTML += "<div onclick=pipelineView('"+getArgument(data, 'jobdir')+"')>View Output From Job</div>";
		}else if (getArgument(data, 'jobtype') == "prime2d" && getArgument(data, 'jobstatus') == "Complete"){
			hoverdiv.innerHTML += "<div onclick=viewFile('"+getArgument(data, 'jobdir')+"/prime2Dcavgs_final_ranked.mrc','" + getArgument(data, 'jobid') + "_" + getArgument(data, 'jobtype') + "')>View Output From Job</div>";
		}
		hoverdiv.innerHTML += "<div onclick=logfileView('"+getArgument(data, 'jobdir')+"/job.log'," + "'" + getArgument(data, 'jobid') + "_" + getArgument(data, 'jobtype') +"')>View Log File</div>";
		hoverdiv.innerHTML += "<div onclick=filesView('"+getArgument(data, 'jobdir')+"'," + "'" + getArgument(data, 'jobid') + "_" + getArgument(data, 'jobtype') +"')>View Files From Job</div>";
		if(getArgument(data, 'taskname') != ""){
			hoverdiv.innerHTML += "<div onclick=rerunJob('"+getArgument(data, 'jobid')+"')>Rerun Job</div>";
		}
		//hoverdiv.innerHTML += "<div onclick=cancelJob('"+getArgument(data, 'jobid')+"')>Cancel Job</div>";
		hoverdiv.innerHTML += "<div onclick=deleteJob('"+getArgument(data, 'jobid')+"')>Delete Job</div>";
		actionmenudiv.appendChild(hoverdiv);
		cell6.appendChild(actionmenudiv);
	}
}

function projectHistoryFull(){
	hideTasks();
	var taskdiv = document.getElementById('taskdiv');
	taskdiv.style.height = "0px";
	var historytable = document.getElementsByClassName('historytable')[0]; 
	historytable.style.height = "calc(100vh - 205px)";
	var historyicon = document.getElementById('historyicon');
	historyicon.setAttribute("src","img/up_chevron.png");
}

function projectHistorySmall(){
	var historytable = document.getElementsByClassName('historytable')[0]; 
	historytable.style.height = "0px";
	var taskdiv = document.getElementById('taskdiv');
	taskdiv.style.height = "calc(100vh - 210px)";
	var historyicon = document.getElementById('historyicon');
	historyicon.setAttribute("src","img/right_chevron.png");
}

function showTask(taskname){
	hideTasks();
	projectHistorySmall();
	var taskpane = document.getElementById(taskname);
	taskpane.style.display = 'block';
	
	var lines = taskpane.getElementsByClassName('taskheaderline');
	for (var i = 0; i < lines.length; i++) { 
		lines[i].style.display = 'table-row';
	}
	
	var lines = taskpane.getElementsByClassName('taskrequiredline');
	for (var i = 0; i < lines.length; i++) { 
		lines[i].style.display = 'table-row';
	}
	
	var lines = taskpane.getElementsByClassName('taskoptionalline');
	for (var i = 0; i < lines.length; i++) { 
		lines[i].style.display = 'none';
	}
}

function hideTasks(){
	var taskpanes = document.getElementsByClassName('taskwindow');
	for (var i = 0; i <  taskpanes.length; i++){
		taskpanes[i].style.display = 'none';
	}
}

function fileBrowser(callingelement){
	showGauze();
	var filesviewcaller = document.getElementById('filesviewcaller');
	filesviewcaller.value = callingelement;
	var filebrowserpopup = document.getElementById('filebrowserpopup');
	filebrowserpopup.style.display = "block";
	fileBrowserWS.execute('cmd=dirlist', "addFileBrowserData");
}

function addFileBrowserData(data){
	if(getArgument(data, 'parentdirectory') != ''){
		addFileBrowserHome(getArgument(data, 'parentdirectory'));
	}else if(getArgument(data, 'directory') != ''){
		addFileBrowserFolderLine(getArgument(data, 'directory'));
	}else if(getArgument(data, 'file') != ''){
		addFileBrowserFileLine(getArgument(data, 'file'));
	}
}

function addFileBrowserHome(foldername){
	var filebrowserdirectorytable = document.getElementById('filebrowserdirectorytable');
	var row = filebrowserdirectorytable.insertRow();
	var cell1 = row.insertCell();
	var cell2 = row.insertCell();
	cell2.innerHTML = "<input id=filebrowsercurrentfolder type=text value=" + foldername + "></input>";
	cell1.innerHTML = "<img src=img/uparrow.png></img>"
	cell1.onclick = function(){
		clearFileBrowserLines();
		var foldersplit = foldername.split('/');
		var parentfolder = '';
		for (var i = 1; i <  (foldersplit.length - 1); i++){
			parentfolder += '/';
			parentfolder += foldersplit[i];
		}
		fileBrowserWS.execute('cmd=dirlist dir='+parentfolder, "addFileBrowserData");
	}
	cell2.onkeydown = function(event){
		if (event.keyCode == 13){
			var filebrowsercurrentfolder= document.getElementById('filebrowsercurrentfolder');
			clearFileBrowserLines();
			fileBrowserWS.execute('cmd=dirlist dir='+filebrowsercurrentfolder.value, "addFileBrowserData");
		}
	}
}

function addFileBrowserFolderLine(foldername){
	var filebrowserlisttable = document.getElementById('filebrowserlisttable');
	var row = filebrowserlisttable.insertRow();
	var cell1 = row.insertCell();
	var cell2 = row.insertCell();
	var folderarray = foldername.split('/');
	cell1.innerHTML = "<img src=img/folder.png></img>";
	cell2.innerHTML = folderarray[folderarray.length - 1];
	row.onclick = function(){
		clearFileBrowserLines();
		fileBrowserWS.execute('cmd=dirlist dir='+foldername, "addFileBrowserData");
	}
}

function addFileBrowserFileLine(filename){
	var filebrowserlisttable = document.getElementById('filebrowserlisttable');
	var row = filebrowserlisttable.insertRow();
	var cell1 = row.insertCell();
	var cell2 = row.insertCell();
	var filearray = filename.split('/');
	cell1.innerHTML = "<img src=img/file.png></img>";
	cell2.innerHTML = filearray[filearray.length - 1];
	cell2.style = "width:100%";
	row.onclick = function(){
		var rows = this.parentElement.childNodes;
		for (var i = 0; i <  rows.length; i++){
			if(rows[i].style.backgroundColor == "rgb(106, 133, 168)"){
				rows[i].removeAttribute("style");
			}
		}
		this.style.backgroundColor = "rgb(106, 133, 168)";
		document.getElementById('filesviewfilename').value = filename;
		document.getElementById('filebrowsercurrentfolder').value = filename;
	}
}

function clearFileBrowserLines(){
	var filebrowserlisttable = document.getElementById('filebrowserlisttable');
	filebrowserlisttable.innerHTML = "";
	var filebrowserdirectorytable = document.getElementById('filebrowserdirectorytable');
	filebrowserdirectorytable.innerHTML = "";
}

function cancelFileBrowser(){
	clearFileBrowserLines();
	var folderselector = document.getElementById('filebrowserpopup');
	folderselector.style.display = 'none';
	hideGauze();
}

function confirmFileBrowser(){ 
	clearFileBrowserLines();
	var filesviewcaller = document.getElementById('filesviewcaller').value;
	if(filesviewcaller == "fileviewer"){
		showTask('fileviewer');
		filesViewFetch();
	} else {
		document.getElementById(filesviewcaller).value = document.getElementById('filesviewfilename').value;
	}
	
//	clearFileBrowserLines();
	var folderselector = document.getElementById('filebrowserpopup');		
	folderselector.style.display = 'none';
	hideGauze();
//	filesViewFetch();
}

function filesViewFetch(){
	var filesviewdiv = document.getElementById('filesviewdiv');
	filesviewdiv.innerHTML = "";
	var filename = document.getElementById('filesviewfilename').value;
	var contrast = document.getElementById('filesviewcontrast').value;
	var brightness = document.getElementById('filesviewbrightness').value;
	console.log(filename);
	filesViewWS.execute("cmd=showfile file=" + filename + " contrast=" + contrast + " brightness=" + brightness, "filesViewFile");
}

function filesViewFile(data){
	if(getArgument(data, 'image') != ""){
		var filesviewdiv = document.getElementById('filesviewdiv');
		var micrographimg = document.createElement("div");
		filesviewdiv.appendChild(micrographimg);
		var image = "url(data:image/png;base64," + getArgument(data, 'image') + ")";
		micrographimg.style.backgroundImage = image;
		micrographimg.style.height = getArgument(data, 'nx') + "px";
		var originalheight = document.getElementById('filesvieworiginalheight');
		originalheight.value = getArgument(data, 'nx');
		micrographimg.style.width = getArgument(data, 'ny') + "px";
		var originalwidth = document.getElementById('filesvieworiginalwidth')
		originalwidth.value = getArgument(data, 'ny');
		micrographimg.className = "micrographimgselected";
		micrographimg.id = getArgument(data, 'nz');
		micrographimg.onclick = function(){
			if(this.className == "micrographimgdeselected"){
				this.className = "micrographimgselected";
			} else {
				this.className = "micrographimgdeselected";
			}
		}
	}else if(getArgument(data, 'log') != ""){
		var filesviewdiv = document.getElementById('filesviewdiv');
		var logdiv = document.createElement("div");
		logdiv.innerHTML = getArgument(data, 'log');
		filesviewdiv.appendChild(logdiv);
	}else if(getArgument(data, 'txt') != ""){
		var filesviewdiv = document.getElementById('filesviewdiv');
		var txtdiv = document.createElement("div");
		txtdiv.innerHTML = getArgument(data, 'txt');
		filesviewdiv.appendChild(txtdiv);
	}
}

function showGauze(){
	var gauze = document.getElementById('gauze');
	gauze.style.display = 'block';
}

function hideGauze(){
	var gauze = document.getElementById('gauze');
	gauze.style.display = 'none';
}

function selectAll(element){
	var container = document.getElementById(element);
	var micrographs = container.children;
	for (var i = 0; i < micrographs.length; i++) {
		micrographs[i].className = "micrographimgselected";
	}
}

function deselectAll(element){
	var container = document.getElementById(element);
	var micrographs = container.children;
	for (var i = 0; i < micrographs.length; i++) {
		micrographs[i].className = "micrographimgdeselected";
	}
}

function selectMode(){
	var modeselector = document.getElementById('modeselector');
	var mode = modeselector.options[modeselector.selectedIndex].value;
	var tasklists = document.getElementsByClassName('tasklist');
	for (var i = 0; i <  tasklists.length; i++){
			tasklists[i].style.display = 'none';
	}
	var tasklist = document.getElementById(mode);
	tasklist.style.display = 'block';
}

function filesViewSave(){
	showGauze();
	document.getElementById('filesviewsave').style.display = "block";
}

function cancelFilesViewSave(){
	document.getElementById('filesviewsave').style.display = "none";
	hideGauze();
}

function confirmFilesViewSave(){
	var filename = document.getElementById('filesviewfilename').value;
	var savefolder = document.getElementById('filesviewsavefolder').value;
	var savefilename = document.getElementById('filesviewsavefilename').value;
	
	var command = 'cmd=savestackselection filename=' + filename + ' outfile=' + savefolder + "/" + savefilename + ' selections=';
	var selections = document.getElementsByClassName('micrographimgselected');
	for(var i = 0; i <  selections.length; i++) {
		command += selections[i].id;
		command += ";";
	}
	fileBrowserWS.execute(command, "cancelFilesViewSave" );
}

function selectFolder(callingelement){
	var folderselector = document.getElementById('folderselectorpopup');
	var folderselectorcaller = document.getElementById('folderselectorcaller');
	folderselectorcaller.value = callingelement;
	showGauze();
	folderselector.style.display = 'block';
	selectFolderWS.execute('cmd=dirlist', "addFolderViewData"); //ADD USER SELECTION
}

function hideSelectFolder(){
	var folderselector = document.getElementById('folderselectorpopup');
	folderselector.style.display = 'none';
	clearFolderLines();
}

function addFolderViewData(data){
	if(getArgument(data, 'parentdirectory') != ''){
		addFolderHome(getArgument(data, 'parentdirectory'));
	}else if(getArgument(data, 'directory') != ''){
		addFolderLine(getArgument(data, 'directory'));
	}
}

function addFolderHome(foldername){
	var folderdirectorytable = document.getElementById('folderdirectorytable');
	var row = folderdirectorytable.insertRow();
	var cell1 = row.insertCell();
	var cell2 = row.insertCell();
	cell2.innerHTML = "<input id=currentfolder type=text value=" + foldername + "></input>";
	cell1.innerHTML = "<img src=img/uparrow.png></img>"
	var folderselectorfolder = document.getElementById('folderselectorfolder');
	folderselectorfolder.value = foldername;
	cell1.onclick = function(){
		clearFolderLines();
		var foldersplit = foldername.split('/');
		var parentfolder = '';
		for (var i = 1; i <  (foldersplit.length - 1); i++){
			parentfolder += '/';
			parentfolder += foldersplit[i];
		}
		selectFolderWS.execute('cmd=dirlist dir='+parentfolder, "addFolderViewData"); //ADD USER SELECTION
	}
	cell2.onkeydown = function(event){
		if (event.keyCode == 13){
			var currentfolder = document.getElementById('currentfolder')
			clearFolderLines();
			selectFolderWS.execute('cmd=dirlist dir='+currentfolder.value, "addFolderViewData"); //ADD USER SELECTION
		}
	}
}

function addFolderLine(foldername){
	var folderlisttable = document.getElementById('folderlisttable');
	var row = folderlisttable.insertRow();
	var cell1 = row.insertCell();
	var cell2 = row.insertCell();
	var folderarray = foldername.split('/');
	cell1.innerHTML = "<img src=img/folder.png></img>";
	cell2.innerHTML = folderarray[folderarray.length - 1];
	cell2.style = "width:100%";
	cell1.onclick = function(){
		clearFolderLines();
		selectFolderWS.execute('cmd=dirlist dir='+foldername, "addFolderViewData"); //ADD USER SELECTION
	}
	cell2.onclick = function(){
		var rows = this.parentElement.parentElement.childNodes;
		for (var i = 0; i <  rows.length; i++){
			if(rows[i].style.backgroundColor == "rgb(106, 133, 168)"){
				rows[i].removeAttribute("style");
			}
		}
		this.parentElement.style.backgroundColor = "rgb(106, 133, 168)";
		var currentfolder = document.getElementById('currentfolder');
		currentfolder.value = foldername;
		var folderselectorfolder = document.getElementById('folderselectorfolder');
		folderselectorfolder.value = foldername;
	}
}

function clearFolderLines(){
	var folderlisttable = document.getElementById('folderlisttable');
	folderlisttable.innerHTML = "";
	var folderdirectorytable = document.getElementById('folderdirectorytable');
	folderdirectorytable.innerHTML = "";
}

function confirmFolderSelection(){
	hideSelectFolder();
	var folderselectorcaller = document.getElementById('folderselectorcaller').value;
	var callingelement = document.getElementById(folderselectorcaller);
	callingelement.value = document.getElementById('folderselectorfolder').value;
	hideGauze();
}

function cancelFolderSelection(){
	hideSelectFolder();
	hideGauze();
}

function filesViewZoom(){
		var filesviewscaleval = document.getElementById('filesviewscale').value;
		var originalheight = document.getElementById('filesvieworiginalheight').value;
		var originalwidth = document.getElementById('filesvieworiginalwidth').value
		var micrographimgselected = document.getElementsByClassName('micrographimgselected');
		var micrographimgdeselected = document.getElementsByClassName('micrographimgdeselected')
		for(var i = 0; i <  micrographimgselected.length; i++) {
			micrographimgselected[i].style.height =  originalheight * (filesviewscaleval/5) + "px";
			micrographimgselected[i].style.width = originalwidth * (filesviewscaleval/5) + "px";
		}
		for(var i = 0; i <  micrographimgdeselected.length; i++) {
			micrographimgdeselected[i].style.height =  originalheight* (filesviewscaleval/5) + "px";
			micrographimgdeselected[i].style.width =  originalwidth* (filesviewscaleval/5) + "px";
		}	
}

function pipelineView(directory){
	var pipelineviewdiv = document.getElementById('pipelineviewdiv');
	pipelineviewdiv.innerHTML = "";
	showTask('pipelineview');
	var pipelineviewdirectory = document.getElementById('pipelineviewdirectory');
	pipelineviewdirectory.value = directory;
	pipelineViewWS.execute("cmd=pipelineviewmaster dir=" + directory, "pipelineViewPopulateArray");
}

function pipelineViewPopulateArray(data){
	if(getArgument(data, 'micrographs') != ""){
		globalmicrographsarray = []; 
		var micrographsarray = getArgument(data, 'micrographs').split(";");
		for (var i = 0; i < micrographsarray.length; i++) {
			var micrographdata = micrographsarray[i].split(",");
			if(micrographdata[0] != ""){
				micrographdata.push("no");
				globalmicrographsarray.push(micrographdata);
			}
		}
		pipelineViewSortColumn(0);
		//globalmicrographsarray.sort();
		
		if(getArgument(data, 'selected') != ""){
			var selectedarray = getArgument(data, 'selected').split(";");
			for (var i = 0; i < selectedarray.length; i++) {
				for (var j = 0; j < globalmicrographsarray.length; j++) {
					if(globalmicrographsarray[j][0] == selectedarray[i]){
						globalmicrographsarray[j][6] = "yes";
					}
				}
			}
		}
		
		pipelineViewAddLines();
	}
}

function pipelineViewAddLines(){
	var pipelineviewdiv = document.getElementById('pipelineviewdiv');
	pipelineviewdiv.innerHTML = "";
		
	for (var i = 0; i < globalmicrographsarray.length; i++) {
		var micrograph = globalmicrographsarray[i];
		var micrographdiv = document.createElement("div");
		var paramsdiv = document.createElement("div");
		var pspecdiv = document.createElement("div");
		var thumbdiv = document.createElement("div");
		
		micrographdiv.appendChild(paramsdiv);
		micrographdiv.appendChild(pspecdiv);
		micrographdiv.appendChild(thumbdiv);
		
		micrographdiv.id = micrograph[0];
		micrographdiv.className = "pipelineviewmicrograph";
		
		if(micrograph[4] == "yes"){
			micrographdiv.setAttribute("name","pipelineviewselected");
		}
		
		paramsdiv.className = 'paramsdiv';
		pspecdiv.className = 'pspecdiv pipelineimg';
		thumbdiv.className = 'thumbdiv pipelineimg';
		
		var paramstable = document.createElement("table");
		paramstable.className = 'paramstable';
		
		var micrographrow = paramstable.insertRow();
		var cell = micrographrow.insertCell();
		cell.innerHTML = "<b>" + micrograph[0] + "</b>";
		cell.colSpan = 2;
		
		var dfxrow = paramstable.insertRow();
		var dfxcell = dfxrow.insertCell()
		dfxcell.innerHTML = "DFX";
		var dfxvaluecell = dfxrow.insertCell()
		dfxvaluecell.innerHTML = micrograph[3]; 
		
		var dfyrow = paramstable.insertRow();
		var dfycell = dfyrow.insertCell()
		dfycell.innerHTML = "DFY";
		var dfyvaluecell = dfyrow.insertCell()
		dfyvaluecell.innerHTML = micrograph[4]; 
		
		var angastrow = paramstable.insertRow();
		var angastcell = angastrow.insertCell()
		angastcell.innerHTML = "Angast";
		var angastvaluecell = angastrow.insertCell()
		angastvaluecell.innerHTML = micrograph[5];
		
		var boxrow = paramstable.insertRow();
		var boxcell = boxrow.insertCell()
		boxcell.innerHTML = "Boxcount";
		var boxvaluecell = boxrow.insertCell()
		boxvaluecell.className = "boxvaluecell";
		
		var viewboxrow = paramstable.insertRow();
		var viewboxcell = viewboxrow.insertCell()
		viewboxcell.innerHTML = "View Boxes";
		var viewboximgcell = viewboxrow.insertCell()
		var viewboximg = document.createElement("img");
		viewboximg.setAttribute("src","img/eye.png");
		viewboximgcell.className = micrograph[0]
		viewboximgcell.appendChild(viewboximg);
		viewboximgcell.onclick = function(){
			pipelineViewMicrograph(this.className);
		}
		
		paramsdiv.appendChild(paramstable);
		
		micrographdiv.onclick = function(){
			if(this.getAttribute("name") == "pipelineviewselected"){
				this.removeAttribute("name");
			} else {
				this.setAttribute("name","pipelineviewselected");
			}
		}
		
		var pspecimg = document.createElement("img");
		var pipelinepspeczoom = document.getElementById('pipelinepspeczoom').value;
		pspecimg.style.transform = "scale("+(1 / Number(pipelinepspeczoom)) + "," +(1 / Number(pipelinepspeczoom)) + ")";
		pspecimg.className = "pspecimg";
		pspecdiv.appendChild(pspecimg);
		
		var thumbimg = document.createElement("img");
		var pipelinethumbzoom = document.getElementById('pipelinethumbzoom').value;
		thumbimg.style.transform = "scale("+(1 / Number(pipelinethumbzoom)) + "," +(1 / Number(pipelinethumbzoom)) + ")";
		thumbimg.className = "thumbimg";
		thumbdiv.appendChild(thumbimg);
		
		pipelineviewdiv.appendChild(micrographdiv);
	}
	
	var pipelinepagenumber = document.getElementById('pipelinepagenumber');
	pipelinepagenumber.value = 1;
	pipelineViewShowPage();
	
}

function pipelineViewAddData(data){
	if(getArgument(data, 'thumb') != ""){
		var line = document.getElementById(getArgument(data, 'thumb'));
		line.getElementsByClassName('thumbimg')[0].setAttribute("src","data:image/png;base64," + getArgument(data, 'image'));
	}
	
	if(getArgument(data, 'pspec') != ""){
		var line = document.getElementById(getArgument(data, 'pspec'));
		line.getElementsByClassName('pspecimg')[0].setAttribute("src","data:image/png;base64," + getArgument(data, 'image'));
	}
	
	if(getArgument(data, 'box') != ""){
		var line = document.getElementById(getArgument(data, 'box'));
		var boxvaluecell = line.getElementsByClassName('boxvaluecell')[0];
		boxvaluecell.innerHTML = getArgument(data, 'count');
	}
	
}

function pipelineViewShowPage(){
	var pipelineview = document.getElementById('pipelineview');
	var pipelinepagenumber = document.getElementById('pipelinepagenumber').value;
	var pipelineviewpagecounter = document.getElementsByName('pipelinepagecount')[0];
	var pipelineviewpagecount = pipelineviewpagecounter.options[pipelineviewpagecounter.selectedIndex].value;
	var pipelineviewtable = document.getElementById('pipelineviewtable');
	var pipelineviewdirectory = document.getElementById('pipelineviewdirectory').value;
	var pipelinepspeccontrast = document.getElementById('pipelinepspeccontrast').value;
	var pipelinepspecbrightness = document.getElementById('pipelinepspecbrightness').value;
	var pipelinethumbcontrast= document.getElementById('pipelinethumbcontrast').value;
	var pipelinethumbbrightness = document.getElementById('pipelinethumbbrightness').value;
	var pipelinepagemax = document.getElementById('pipelinepagemax');
	var pipelineviewmicrographs = document.getElementsByClassName('pipelineviewmicrograph');
	pipelinepagemax.innerHTML = Math.ceil(pipelineviewmicrographs.length / pipelineviewpagecount);
	var command = "cmd=pipelineviewslave dir=" + pipelineviewdirectory + 
				  " pspecbrightness=" + pipelinepspecbrightness + 
				  " pspeccontrast=" + pipelinepspeccontrast +
				  " thumbbrightness=" + pipelinethumbbrightness + 
				  " thumbcontrast=" + pipelinethumbcontrast +
				  "  micrographs=";
	for (var i = 0; i < pipelineviewmicrographs.length; i++) {
		pipelineviewmicrographs[i].style.display = "none";
		pipelineviewmicrographs[i].getElementsByClassName('thumbimg')[0].removeAttribute("src");
		pipelineviewmicrographs[i].getElementsByClassName('pspecimg')[0].removeAttribute("src");
	}
	
	for (var i = 0; i < pipelineviewmicrographs.length; i++) {
		if( i < (pipelinepagenumber * pipelineviewpagecount) && i >= ((pipelinepagenumber - 1) * pipelineviewpagecount)){
			pipelineviewmicrographs[i].style.display = null;
			command += pipelineviewmicrographs[i].id;
			command += ";"
		}
	}
	
	var selections = document.getElementsByName('pipelineviewselected');
	command += " selections=";
	for(var i = 0; i <  selections.length; i++) {
		command += selections[i].id;
		command += ";";
	}
	
	pipelineview.scrollTop = 0;
	pipelineViewWS.execute(command, "pipelineViewAddData");
}

function pipelineViewPageUp(){
	var pipelinepagenumber = document.getElementById('pipelinepagenumber');
	var currentpage = pipelinepagenumber.value;
	pipelinepagenumber.value = Number(currentpage) + 1;
	pipelineViewShowPage();
}

function pipelineViewPageDown(){
	var pipelinepagenumber = document.getElementById('pipelinepagenumber');
	var currentpage = pipelinepagenumber.value;
	pipelinepagenumber.value = Number(currentpage) - 1;
	pipelineViewShowPage();
}

function pipelineViewUpdateZoom(){
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

function pipelineViewSave(){
	showGauze();
	document.getElementById('pipelineviewsave').style.display = "block";
}

function pipelineViewSort(sortelement){ 
	var pipelinepagenumber = document.getElementById('pipelinepagenumber');
	pipelinepagenumber.value = 1;

	if(sortelement == "name"){
		pipelineViewSortColumn(0);
	} else if (sortelement == "dfx"){
		pipelineViewSortColumn(3);
	} else if (sortelement == "dfy"){
		pipelineViewSortColumn(4);
	} else if (sortelement == "angast"){
		pipelineViewSortColumn(5);
	}
	pipelineViewAddLines();
}

function pipelineViewSortColumn(i){ 
	globalmicrographsarray.sort(function(a,b) {
        return a[i]-b[i]
    });
}

function pipelineViewSelectAll(){
	var pipelineviewtable = document.getElementById('pipelineview');
	var micrographs = pipelineviewtable.getElementsByClassName('pipelineviewmicrograph');
	for (var i = 0; i < micrographs.length; i++) {
		micrographs[i].setAttribute("name","pipelineviewselected");
	}
}

function pipelineViewClearAll(){
	var pipelineviewtable = document.getElementById('pipelineview');
	var micrographs = pipelineviewtable.getElementsByClassName('pipelineviewmicrograph');
	for (var i = 0; i < micrographs.length; i++) {
		micrographs[i].removeAttribute("name");
	}
}

function cancelPipelineViewSave(){
	hideGauze();
	document.getElementById('pipelineviewsave').style.display = "none";
}

function confirmPipelineViewSave(){
	var dir = document.getElementById('pipelineviewdirectory').value;
	var pipelineviewsavefilename = document.getElementById('pipelineviewsavefilename');
	var filename = pipelineviewsavefilename.value;
	var pipelineviewsavefolder = document.getElementById('pipelineviewsavefolder');
	var outdir = pipelineviewsavefolder.value;
	var command = 'cmd=pipelinesaveselection dir=' + dir + ' filename=' + outdir + '/' + filename + ' selections=';
	var selections = document.getElementsByName('pipelineviewselected');
	for(var i = 0; i <  selections.length; i++) {
		command += selections[i].id;
		command += ";";
	}
	pipelineViewWS.execute(command, "cancelPipelineViewSave" );
	pipelineviewsavefilename.value = "";
	pipelineviewsavefolder.value = "";
	
}

function pipelineViewMicrograph(micrograph){
	document.getElementById('pipelineviewboxselect').style.display = "block";
	var pipelineviewboxselectmicrograph = document.getElementById('pipelineviewboxselectmicrograph');
	pipelineviewboxselectmicrograph.innerHTML = "";
	var pipelineviewboxname = document.getElementById('pipelineviewboxname');
	pipelineviewboxname.innerHTML = micrograph;
	var pipelinethumbcontrast= document.getElementById('pipelinethumbcontrast').value;
	var pipelinethumbbrightness = document.getElementById('pipelinethumbbrightness').value;
	
	var dir = document.getElementById('pipelineviewdirectory').value;
	pipelineViewWS.execute("cmd=pipelineviewboxes dir=" + dir + " thumbbrightness=" + pipelinethumbbrightness + " thumbcontrast=" + pipelinethumbcontrast + " micrograph=" + micrograph, "pipelineViewAddMicrograph");
}

function pipelineHideMicrograph(){
	document.getElementById('pipelineviewboxselect').style.display = "none";
}

function pipelineViewAddMicrograph(data){
	if(getArgument(data, 'image') != ""){
		var pipelineviewboxselectmicrograph = document.getElementById('pipelineviewboxselectmicrograph');
		pipelineviewboxselectmicrograph.innerHTML = "";
		var micrographimg = document.createElement("img");
		pipelineviewboxselectmicrograph.appendChild(micrographimg);
		micrographimg.setAttribute("src","data:image/png;base64," + getArgument(data, 'image'));
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
		console.log(boxarray[0] + "," + boxarray[1]);
	}
	pipelineViewClose();
}

function hideNewProjectPopup(){
	var newprojectpopup = document.getElementById('newprojectpopup');
	newprojectpopup.style.display = 'none';
	clearNewProjectPopupLines();
}

function clearNewProjectPopupLines(){
	var newprojectname = document.getElementById('newprojectname');
	newprojectname.value = "";
	var newprojectfolder = document.getElementById('newprojectfolder');
	newprojectfolder.value = "";
}

function confirmNewProjectPopup(){
	var newprojectname = document.getElementById('newprojectname').value;
	var newprojectfolder = document.getElementById('newprojectfolder').value;
	newProjectWS.execute('cmd=addproject projectname=' + newprojectname + ' projectdir=' + newprojectfolder, "pageReload");
	hideGauze();
	hideNewProjectPopup();
	projectSelector();
}

function pageReload(){
	location.reload(); 
}

function cancelNewProjectPopup(){
	hideNewProjectPopup();
	hideGauze();
}

function rerunJob(jobid){
	var projectselector = document.getElementById('projectselector');
	var project = projectselector.options[projectselector.selectedIndex].value;	
	projectHistoryWS.execute("cmd=jobparameters projectname="+project+" jobid="+jobid, "rerunJobProcess");
}

function rerunJobProcess(data){
	if(getArgument(data, 'taskname') != ""){
		showTask(getArgument(data, 'taskname'));
		var argsarray = data.split(" ");
		var taskpane = document.getElementById(getArgument(data, 'taskname'));
		for (var i = 0; i < argsarray.length; i++) {
			var key = argsarray[i].split("=")[0];
			var value = argsarray[i].split("=")[1].replace(/\^/g, " ");
			var jobargs = taskpane.getElementsByClassName('jobarg');
			for (var j = 0; j < jobargs.length; j++) {
				if(jobargs[j].id ==  key){
					jobargs[j].value = value;
				}
			}
		}
	}
}

function cancelJob(jobid){
	var projectselector = document.getElementById('projectselector');
	var project = projectselector.options[projectselector.selectedIndex].value;	
	projectHistoryWS.execute("cmd=killjob projectname="+project+" jobid="+jobid, "selectProject");
}

function deleteJob(jobid){
	var projectselector = document.getElementById('projectselector');
	var project = projectselector.options[projectselector.selectedIndex].value;	
	projectHistoryWS.execute("cmd=deletejob projectname="+project+" jobid="+jobid, "selectProject");
}

function runTask(task){
	var commandstring = "";
	var cleanstring = "";
	var taskwindow = document.getElementById(task);
	var jobargs = taskwindow.getElementsByClassName('jobarg');
	for(var i = 0; i < jobargs.length; i++){
		if(jobargs[i].value != ""){
			cleanstring = jobargs[i].value.replace(/ /g, "^");
			commandstring += jobargs[i].id + "=" + cleanstring + " ";
		}
	}
	
	var projectselector = document.getElementById('projectselector');
	var project = projectselector.options[projectselector.selectedIndex].value;
	
	var jobrunparameterspopup = document.getElementById('jobrunparameterspopup');
	var jobargs = jobrunparameterspopup.getElementsByClassName('jobarg');
	
	for(var i = 0; i < jobargs.length; i++){
		if(jobargs[i].value != ""){
			cleanstring = jobargs[i].value.replace(/ /g, "^");
			commandstring += jobargs[i].id + "=" + cleanstring + " ";
		}
	}
	
	commandstring += "projectname=" + project;
	runTaskWS.execute(commandstring, "Null");
	var xmlHttp = new XMLHttpRequest();
	var url = "http://simplecryoem.com/tracking/index.php?jobtype=" + getArgument(commandstring, "jobtype");
	xmlHttp.open( "GET", url, true );
    xmlHttp.send( null );
	hideRunParameters();
	
	var jobargs = document.getElementsByClassName('jobarg');
	for(var i = 0; i < jobargs.length; i++){
		jobargs[i].value = ""; 
	}
	
	selectProject();
	projectHistoryFull();
}

function taskHeaderShowHide(element){
	var parent = element.parentElement;
	var lines = parent.getElementsByClassName('taskheaderline');
	for (var i = 0; i < lines.length; i++) { 
		if(lines[i].style.display == 'table-row'){
			lines[i].style.display = 'none';
		}else{
			lines[i].style.display = 'table-row';
		}
	}
}

function taskRequiredShowHide(element){
	var parent = element.parentElement;
	var lines = parent.getElementsByClassName('taskrequiredline');
	for (var i = 0; i < lines.length; i++) { 
		if(lines[i].style.display == 'table-row'){
			lines[i].style.display = 'none';
		}else{
			lines[i].style.display = 'table-row';
		}
	}
}

function taskOptionalShowHide(element){
	var parent = element.parentElement;
	var lines = parent.getElementsByClassName('taskoptionalline');
	for (var i = 0; i < lines.length; i++) { 
		if(lines[i].style.display == 'table-row'){
			lines[i].style.display = 'none';
		}else{
			lines[i].style.display = 'table-row';
		}
	}
}

function showRunParameters(task){
	showGauze();
	var jobrunparameterspopup = document.getElementById('jobrunparameterspopup');
	jobrunparameterspopup.style.display = 'block';
	var jobrunparametersrunbutton = document.getElementById('jobrunparametersrunbutton');
	jobrunparametersrunbutton.onclick = function(){
		runTask(task);
	};
}

function hideRunParameters(){
	var jobrunparameterspopup = document.getElementById('jobrunparameterspopup');
	jobrunparameterspopup.style.display = 'none';
	hideGauze();
}

function filesView(){
	showTask('filesviewtask');
	var filesviewtable = document.getElementById('viewjobfilestable');
	filesviewtable.innerHTML = "";
	var filesviewdescription = document.getElementById('filesviewdescription');
	filesviewdescription.innerHTML = "View Files From Job: " + arguments[1];
	filesViewWS.execute("cmd=dirlist dir=" + arguments[0], "filesViewAddData");
}

function filesViewAddData(data){
	if(getArgument(data, 'directory') != ''){
		filesViewAddFolderLine(getArgument(data, 'directory'));
	}else if(getArgument(data, 'file') != ''){
		filesViewAddFileLine(getArgument(data, 'file'));
	}
}

function filesViewAddFolderLine(foldername){
	var filesviewtable = document.getElementById('viewjobfilestable');
	var row = filesviewtable.insertRow();
	var cell1 = row.insertCell();
	var cell2 = row.insertCell();
	var folderarray = foldername.split('/');
	cell1.innerHTML = "<img src=img/folder.png></img>";
	cell2.innerHTML = folderarray[folderarray.length - 1];
	cell1.onclick = function(){
		clearFilesViewLines();
		filesViewWS.execute('cmd=dirlist dir='+foldername, "filesViewAddData");
	}
}

function filesViewAddFileLine(filename){
	var filesviewtable = document.getElementById('viewjobfilestable');
	var row = filesviewtable.insertRow();
	var cell1 = row.insertCell();
	var cell2 = row.insertCell();
	var filearray = filename.split('/');
	cell1.innerHTML = "<img src=img/file.png></img>";
	cell2.innerHTML = filearray[filearray.length - 1];
	cell2.onclick = function(){
		viewFile(filename, filename);
		//var filesviewfilename = document.getElementById('filesviewfilename');
		//filesviewfilename.value = filename;
		//showTask('fileviewer');
		//filesViewFetch();
	}
}

function logfileView(){
	var filesviewfilename = document.getElementById('filesviewfilename');
	filesviewfilename.value = arguments[0];
	var fileviewerdescription = document.getElementById('fileviewerdescription');
	fileviewerdescription.innerHTML = "<b>View Log File: " + arguments[1] + "</b>";
	showTask('fileviewer');
	filesViewFetch();
}

function viewFile(){
	var filesviewfilename = document.getElementById('filesviewfilename');fileviewerdescription
	filesviewfilename.value = arguments[0];
	var fileviewerdescription = document.getElementById('fileviewerdescription');
	fileviewerdescription.innerHTML = "<b>View Files From Job: " + arguments[1] + "</b>";
	showTask('fileviewer');
	filesViewFetch();
}

includeHTML();
projectSelector();
setInterval(function(){refreshProject()}, 30000);
