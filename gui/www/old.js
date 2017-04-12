
function statusPopup(data){
	var status = getArgument(data, 'status');
	if(status == "complete"){
		console.log("complete")
	}else{
		console.log("Error");
		console.log(data);
	}
}

function filesView(directory){
	showTask('filesviewtask');
	var filesviewdirectory = document.getElementById('filesviewdirectory');
	filesviewdirectory.value = directory;
	var filesviewtable = document.getElementById('filesviewtable');
	filesviewtable.innerHTML = "";
	filesViewWS.execute("cmd=dirlist dir=" + directory, "filesViewAddData");
}

function filesViewAddData(data){
	if(getArgument(data, 'directory') != ''){
		filesViewAddFolderLine(getArgument(data, 'directory'));
	}else if(getArgument(data, 'file') != ''){
		filesViewAddFileLine(getArgument(data, 'file'));
	}
}

function filesViewAddFolderLine(foldername){
	var filesviewtable = document.getElementById('filesviewtable');
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
	var filesviewtable = document.getElementById('filesviewtable');
	var row = filesviewtable.insertRow();
	var cell1 = row.insertCell();
	var cell2 = row.insertCell();
	var filearray = filename.split('/');
	cell1.innerHTML = "<img src=img/file.png></img>";
	cell2.innerHTML = filearray[filearray.length - 1];
	cell2.onclick = function(){
		var filesviewfilename = document.getElementById('filesviewfilename');
		filesviewfilename.value = filename;
		filesViewFetch();
		// Add select colour change here
	}
}





function filesViewBack(){
	var filesviewtable = document.getElementById('filesviewtable');
	filesviewtable.style.display = "block";
	var filesviewdiv = document.getElementById('filesviewdiv');
	filesviewdiv.innerHTML = "";
}





function Null(data){
	console.log(data);
}
