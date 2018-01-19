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

function get2DViewData () {
 getAjax('JSONhandler?function=2dview&' + window.location.toString().split('?')[1], function(data){view2DInit(data)});
}

function showParticleData (data){
	var JSONdata = JSON.parse(data);
	var particles = document.getElementById('particles');
	particles.innerHTML = "";
	for (var i = 0; i < JSONdata.snapshots.length; i++) {
		var div = document.createElement("div");
		div.setAttribute('data-frameid', JSONdata.snapshots[i].frameid);
		div.setAttribute('data-src', JSONdata.snapshots[i].stackfile);
		div.className = 'particle';
		particles.appendChild(div);
	} 
	updateParticleViewer();
}

function updateParticleViewer () {
	var particles = document.getElementsByClassName('particle');
	var rootdirectory = document.getElementById('rootdirectory').value;
	var particlecontrast = document.getElementById('particlecontrast').value;
	var particlebrightness= document.getElementById('particlebrightness').value;
	var particlesize= document.getElementById('particlesize').value;
	
	for (var i = 0; i < particles.length; i++){
		var img = document.createElement("img");
		img.setAttribute('data-src',"JPEGhandler?filename=" + particles[i].getAttribute('data-src') + "&contrast=" +particlecontrast+"&brightness="+particlebrightness+"&size="+particlesize+"&frameid=" + particles[i].getAttribute('data-frameid'));
		img.style.width = particlesize+"px";
		img.style.height = particlesize+"px";
		img.alt = "Loading ...";
		img.setAttribute('data-frameid', particles[i].getAttribute('data-frameid'));
		particles[i].appendChild(img);
	}
	document.getElementById('particles').addEventListener('scroll', function(){loadParticles2()});
	loadParticles2();
}

function loadParticles(){
	//var unloaded = document.querySelectorAll('[data-loaded="no"]')[0];
	var unloaded = document.querySelectorAll('[data-particleshow="yes"]')[0];
	unloaded.removeAttribute('data-particleshow');
	unloaded.addEventListener('load', loadParticles)
	unloaded.src = unloaded.getAttribute('data-src');
}

function loadParticles2(){

	var particles= document.getElementsByClassName('particle');
	var pane = document.getElementById('particles');
	paneBoundingClientRect = pane.getBoundingClientRect();
	var topcutoff = paneBoundingClientRect.top;
	var bottomcutoff = paneBoundingClientRect.bottom + (paneBoundingClientRect.bottom - paneBoundingClientRect.top);
	for (var i = 0; i < particles.length; i++){
		var snapshotBoundingClientRect = particles[i].getBoundingClientRect();
		var images = particles[i].getElementsByTagName('img');
		if(snapshotBoundingClientRect.bottom > topcutoff && snapshotBoundingClientRect.top < bottomcutoff) {
			images[0].setAttribute('data-particleshow',"yes");
		} else 	{
			images[0].removeAttribute('src');
			images[0].removeAttribute('data-particleshow');
			images[0].removeEventListener("load", loadImages);
		}
	}
	loadParticles();
}

function updateParticlesControls () {
	var particles= document.getElementsByClassName('particle');
	var rootdirectory = document.getElementById('rootdirectory').value;
	var particlecontrast = document.getElementById('particlecontrast').value;
	var particlebrightness= document.getElementById('particlebrightness').value;
	var particlesize= document.getElementById('particlesize').value;
	
	for (var i = 0; i < particles.length; i++){
		var img = particles[i].getElementsByTagName('img')[0];
		img.setAttribute('data-loaded',"no");
		img.setAttribute('data-src',"JPEGhandler?filename=" + particles[i].getAttribute('data-src') + "&contrast=" +particlecontrast+"&brightness="+particlebrightness+"&size="+particlesize+"&frameid=" + particles[i].getAttribute('data-frameid'));
		img.style.width = particlesize+"px";
		img.style.height = particlesize+"px"
	}
	loadParticles2();
}


get2DViewData();
