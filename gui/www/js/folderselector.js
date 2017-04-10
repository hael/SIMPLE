function selectFolder(callingelement){
	selectfolder_caller = callingelement;
	var folderselector = document.getElementById('folderselectorpopup');
	showGauze();
	folderselector.style.display = 'block';
	// Add twirl here
	folderselection = '';
	ws.send('cmd=dir usr=local dir=home'); //ADD USER SELECTION
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
	cell1.onclick = function(){
		clearFolderLines();
		var foldersplit = foldername.split('/');
		var parentfolder = '';
		for (var i = 1; i <  (foldersplit.length - 1); i++){
			parentfolder += '/';
			parentfolder += foldersplit[i];
		}
		ws.send('cmd=dir user=local dir='+parentfolder);//ADD USER SELECTION
	}
	cell2.onkeydown = function(event){
		if (event.keyCode == 13){
			var currentfolder = document.getElementById('currentfolder')
			clearFolderLines();
			ws.send('cmd=dir user=local dir='+currentfolder.value);//ADD USER SELECTION
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
	cell1.onclick = function(){
		clearFolderLines();
		ws.send('cmd=dir user=local dir='+foldername);//ADD USER SELECTION
	}
	cell2.onclick = function(){
		selectfolder_value = foldername;
		// Add select colour change here
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
	var callingelement = document.getElementById(selectfolder_caller);
	callingelement.value = selectfolder_value;
	hideGauze();
}

function cancelFolderSelection(){
	hideSelectFolder();
	hideGauze();
}
	
