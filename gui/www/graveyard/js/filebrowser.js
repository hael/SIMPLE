function fileBrowser(){
	var fileBrowser = document.getElementById('filebrowserpopup');
	showGauze();
	fileBrowser.style.display = 'block';
	// Add twirl here
	ws_filebrowser = new WebSocket('ws://localhost:8080/');
	ws_filebrowser.onopen = function(){
		ws.send('cmd=listdir usr=' + USER);
	}
	ws_filebrowser.onerror = function(event){
		console.log(event.data);
	}
	ws_filebrowser.onmessage = function(event) {
		if (getArgument(event.data, 'parentdirectory') != ''){
			addFileBrowserHome(getArgument(event.data, 'parentdirectory'));
		} else if (getArgument(event.data, 'directory') != ''){
			addFileBrowserFolderLine(getArgument(event.data, 'directory'));
		} else if (getArgument(event.data, 'file') != ''){
			addFileBrowserFileLine(getArgument(event.data, 'file'));
		}
	};
}

function hideFileBrowser(){
	var folderselector = document.getElementById('filebrowserpopup');
	folderselector.style.display = 'none';
	clearFolderLines();
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
		ws_filebrowser.send('cmd=listdir loc='+parentfolder);
	}
	cell2.onkeydown = function(event){
		if (event.keyCode == 13){
			var filebrowsercurrentfolder= document.getElementById('filebrowsercurrentfolder');
			clearFileBrowserLines();
			ws_filebrowser.send('cmd=listdir loc='+filebrowsercurrentfolder.value);
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
	cell1.onclick = function(){
		clearFileBrowserLines();
		ws_filebrowser.send('cmd=listdir loc='+foldername);
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
	cell2.onclick = function(){
		filebrowser_value = filename;
		// Add select colour change here
	}
}

function clearFileBrowserLines(){
	var filebrowserlisttable = document.getElementById('filebrowserlisttable');
	filebrowserlisttable.innerHTML = "";
	var filebrowserdirectorytable = document.getElementById('filebrowserdirectorytable');
	filebrowserdirectorytable.innerHTML = "";
}

function confirmFileBrowser(){
	clearFileBrowserLines();
	ws_filebrowser.close();
	showTask('fileviewertask');
	hideFileBrowser();
	hideGauze();
	ws_fileviewer = new WebSocket('ws://localhost:8080/');
	ws_fileviewer.onopen = function(){
		ws_fileviewer.send('cmd=viewfile file=/crescent/lea/Ribosome/pickrefs_ribo.mrc'); //change to filebrowser_value
	}
	ws_fileviewer.onerror = function(event){
		console.log(event.data);
	}
	ws_fileviewer.onmessage = function(event) {
		if (getArgument(event.data, 'count') != ''){
			var fileviewertask = document.getElementById('fileviewertask');
			var image = document.createElement("img");
			image.setAttribute("src","data:image/png;base64," + getArgument(event.data, 'image'));
			fileviewertask.appendChild(image);
		}
	};
}

function cancelFileBrowser(){
	clearFileBrowserLines();
	ws_filebrowser.close();
	hideFileBrowser();
	hideGauze();
}
	
