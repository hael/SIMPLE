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

function getLogFileviewerData () {
 
 getAjax('JSONhandler?function=viewlogfile&'+ window.location.toString().split('?')[1], function(data){showLogFileViewerData(data)});
}
 

function showLogFileViewerData (data){
	var JSONdata = JSON.parse(data);
	var snapshots = document.getElementById('snapshots');
	var logfile = JSONdata.logfile;
	var logfileformatted;
	logfileformatted = logfile.replace(/\^/g, "\\");
	logfileformatted = logfile.replace(/@/g, "<br>");
	snapshots.innerHTML = logfileformatted;
}

function showProjectHistory(){
	var mainpaneiframe = parent.document.getElementById('mainpaneiframe');
	mainpaneiframe.src = "projecthistory.html"
}

getLogFileviewerData();
