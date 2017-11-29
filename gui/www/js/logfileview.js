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
 
 getAjax('JSONhandler?function=viewlogfile&='+ window.location.toString().split('?')[1], function(data){showLogFileViewerData(data)});
}
 

function showLogFileViewerData (data){

}

function showProjectHistory(){
	var mainpaneiframe = parent.document.getElementById('mainpaneiframe');
	mainpaneiframe.src = "projecthistory.html"
}

getLogFileviewerData();
