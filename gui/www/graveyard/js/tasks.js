
function scrubSpaces(string){
	var returnstring = string.replace(/\s+/g, "^");
	return returnstring;
}

function submitTask(parentid, type){
	var project = document.getElementById('projectselector').value;
	var element = document.getElementById(parentid);
	var jobname = element.getElementsByClassName('jobname')[0].value;
	var jobdescription = element.getElementsByClassName('jobdescription')[0].value;
	var jobtype = element.getElementsByClassName('jobtype')[0].value;
	var externaljobfolder;
	
	if(jobtype == "external"){
		jobtype = element.getElementsByClassName('externaljobtype')[0].value;
		externaljobfolder = element.getElementsByClassName('externaljobfolder')[0].value;
	}
	
	var commandstring = "cmd=addjob " + 
						 "usr=local " +
						 "project=" + scrubSpaces(project) + " " +
						 "jobtype=" + scrubSpaces(jobtype) + " " +
						 "jobdir=" + scrubSpaces(externaljobfolder) + " " +
						 "jobdescription=" + scrubSpaces(jobdescription) + " " + 
						 "jobname=" + scrubSpaces(jobname);
	console.log(commandstring);
	ws.send(commandstring);	
}
