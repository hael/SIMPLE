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

function setTitle() {
	var headerproject = document.getElementById('headerproject');
	var projectselector = parent.document.getElementById('projectselector');
	var project = projectselector.options[projectselector.selectedIndex].getAttribute('data-projectname');
	if (project){
		headerproject.innerHTML = projectselector.options[projectselector.selectedIndex].getAttribute('data-projectname');
	} else {
		document.getElementById('header').innerHTML = "";
	}
	var projectdescription = projectselector.options[projectselector.selectedIndex].getAttribute('data-projectdescription');
	if (projectdescription){
		headerproject.title = decodeURIComponent(projectselector.options[projectselector.selectedIndex].getAttribute('data-projectdescription'));
	}
	
}

function getJobs () {
	var projectselector = parent.document.getElementById('projectselector');
	getAjax('JSONhandler?function=getjobs&table=' + projectselector.options[projectselector.selectedIndex].getAttribute('data-projecttable'), function(data){addJobs(data)});
}

function addJobs(data){
	var JSONdata = JSON.parse(data);
	var projectselector = parent.document.getElementById('projectselector');
	var historytable = document.getElementById('historytable');
	for (var i = 0; i < JSONdata.jobs.length; i++) {
		var row = historytable.insertRow();
		var cell1 = row.insertCell();
		var cell2 = row.insertCell();
		var cell3 = row.insertCell();
		var cell4 = row.insertCell();
		var cell5 = row.insertCell();
		cell1.innerHTML = JSONdata.jobs[i].id;
		cell2.innerHTML = JSONdata.jobs[i].jobtype;
		cell3.innerHTML = JSONdata.jobs[i].jobstatus;
		cell3.innerHTML = JSONdata.jobs[i].jobstatus;
		cell4.innerHTML = JSONdata.jobs[i].jobname;
		cell4.title = decodeURIComponent(JSONdata.jobs[i].jobdescription);
		cell4.style.width = "100%";
		var gearimage = document.createElement("img");
		gearimage.src = "img/gear.png";
		gearimage.onclick = function(){
			var jobmenu = this.parentNode.getElementsByClassName('jobmenu')[0];
			if(jobmenu.style.display == "block"){
				jobmenu.style.display = "none";
			} else {
				jobmenu.style.display = "block";
			}
		}
		cell5.appendChild(gearimage);
		var jobmenu = document.createElement("div");
		jobmenu.className = "jobmenu";
		
		var viewoutput = document.createElement("div");
		viewoutput.innerHTML = "View Output";
		viewoutput.setAttribute('data-jobfolder', JSONdata.jobs[i].jobfolder);
		viewoutput.setAttribute('data-jobtype', JSONdata.jobs[i].jobtype);
		viewoutput.onclick = function(){viewOutput(this)};
		jobmenu.appendChild(viewoutput);
		
		var viewlogfile = document.createElement("div");
		viewlogfile.innerHTML = "View Log File";
		viewlogfile.setAttribute('data-jobfolder', JSONdata.jobs[i].jobfolder);
		viewlogfile.onclick = function(){viewLogfile(this)};
		jobmenu.appendChild(viewlogfile);
		
		var rerunjob = document.createElement("div");
		
		rerunjob.innerHTML = "Rerun Job";
		rerunjob.setAttribute('data-jobrerun', JSONdata.jobs[i].jobrerun);
		rerunjob.setAttribute('data-jobtype', JSONdata.jobs[i].jobtype);
		rerunjob.onclick = function(){rerunJob(this)};
		jobmenu.appendChild(rerunjob);
		
		var deletejob = document.createElement("div");
		
		deletejob.innerHTML = "Delete Job";
		deletejob.setAttribute('data-jobid', JSONdata.jobs[i].id);
		deletejob.setAttribute('data-projecttable', projectselector.options[projectselector.selectedIndex].getAttribute('data-projecttable'));
		deletejob.onclick = function(){deleteJob(this)};
		jobmenu.appendChild(deletejob);
		
		cell5.appendChild(jobmenu);
	}
}

function rerunJob(element){		
	var mainpaneiframe = parent.document.getElementById('mainpaneiframe');
	mainpaneiframe.src = "jobselector.html?" + element.getAttribute('data-jobrerun'); 
}

function viewOutput(element){		
	var mainpaneiframe = parent.document.getElementById('mainpaneiframe');
	if (element.getAttribute('data-jobtype') == "preproc"){
		mainpaneiframe.src = "pipelineview.html?folder=" + element.getAttribute('data-jobfolder');
	} else if (element.getAttribute('data-jobtype') == "stream2d"){
		mainpaneiframe.src = "2dview.html?folder=" + element.getAttribute('data-jobfolder');
	} else if (element.getAttribute('data-jobtype') == "prime2D_stream"){
		mainpaneiframe.src = "2dview.html?folder=" + element.getAttribute('data-jobfolder')
	} else if (element.getAttribute('data-jobtype') == "prime2D"){
		mainpaneiframe.src = "2dview.html?folder=" + element.getAttribute('data-jobfolder');
	} else if (element.getAttribute('data-jobtype') == "unblur"){
		mainpaneiframe.src = "unblurview.html?folder=" + element.getAttribute('data-jobfolder');
	} else if (element.getAttribute('data-jobtype') == "ctffind"){
		mainpaneiframe.src = "ctffindview.html?folder=" + element.getAttribute('data-jobfolder');
	} else if (element.getAttribute('data-jobtype') == "ctffit"){
		mainpaneiframe.src = "ctffindview.html?folder=" + element.getAttribute('data-jobfolder');	
	}else if (element.getAttribute('data-jobtype') == "extract"){
		mainpaneiframe.src = "particleview.html?folder=" + element.getAttribute('data-jobfolder');
	}else if (element.getAttribute('data-jobtype') == "prime3D"){
		mainpaneiframe.src = "3dview.html?folder=" + element.getAttribute('data-jobfolder');
	}else if (element.getAttribute('data-jobtype') == "ini3D"){
		console.log("JILJ");
		mainpaneiframe.src = "ini3dview.html?folder=" + element.getAttribute('data-jobfolder');
	}
}

function viewLogfile(element){		
	var mainpaneiframe = parent.document.getElementById('mainpaneiframe');
	mainpaneiframe.src = "logfileview.html?folder=" + element.getAttribute('data-jobfolder');

}

function deleteJob(element){
	
	var url = 'JSONhandler?function=deletejob&table=' + element.getAttribute('data-projecttable') + '&jobid=' + element.getAttribute('data-jobid');
	
	if (confirm("Do you wish to delete this job? Please note, you will have to remove the job directory yourself") == true) {
		getAjax(url, function(data){window.location.reload(1)});
	}
}

setTitle();
getJobs();

setTimeout(function(){
   window.location.reload(1);
}, 20000);
