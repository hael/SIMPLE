var sourcefiles;
var targetflags;
var targetscale;
var rootdir;
var loadqueue;
var currentproject;
	  
function getAjax (url, success) {
	var xhr = new XMLHttpRequest();
	xhr.open('GET', url, true);
	xhr.onreadystatechange = function() {
		if (xhr.readyState>3 && xhr.status==200){
			var JSONdata = JSON.parse(xhr.responseText);
			// TEST FOR ERROR STING IN JSON
			success(JSONdata);
		}
	};
	xhr.send();
}
	  
function addSideMenuTab(elementName, elementFunction){
	var div = document.createElement("div");
	div.className = "sidemenutab";
	div.innerHTML = elementName;
	div.onclick = function(){
		var sidemenudisplay = document.getElementById('sidemenudisplay');
		if(this.className == "sidemenutabselected"){
			this.className = "sidemenutab";
			sidemenudisplay.innerHTML = "";
			sidemenudisplay.style.display = "none";
		} else {
			this.className = "sidemenutabselected";
			sidemenudisplay.style.display = "block";
			elementFunction();			
		}
	}
	document.getElementById('sidemenutabs').appendChild(div);
}
  
  function createSideMenuDisplay(title) {
	  var sidemenudisplay = document.getElementById('sidemenudisplay');
	  sidemenudisplay.innerHTML = "";
	  var header = document.createElement("div");
	  var simpleicon = document.createElement("img");
	  header.className = "sidemenuheader";
	  simpleicon.className = "simpleicon";
	  simpleicon.src = 'img/simple_favicon.ico';
	  header.appendChild(simpleicon);
	  header.innerHTML += title;
	  sidemenudisplay.appendChild(header);
	  return sidemenudisplay;
  }
  
function createInput(type, label, id, minimum, maximum, step){
	  var div = document.createElement("div");
	var input = document.createElement("input");
	  var helpicon = document.createElement("img");
	  input.id = id;
	  input.type = type;
	  input.min = minimum;
	  input.max = maximum;
	  input.step = step;
	  div.innerHTML = label;
	  div.appendChild(input);
	  div.appendChild(helpicon);
	 return div;
}
  
  function createSelect(label, id){
	  var div = document.createElement("div");
	  var select = document.createElement("select");
	  var helpicon = document.createElement("img");
	  select.id = id;
	  div.innerHTML = label;
	  div.appendChild(select);
	  div.appendChild(helpicon);
	 return div;
  }
  
  function createOption(label, value){
	  var option = document.createElement("option");
	  option.innerHTML = label;
	  option.value = value
	 return option;
  }
  
  function createButton(label, elementFunction){
	  var button = document.createElement("button");
	  button.innerHTML = label;
	  button.onclick = elementFunction;
	 return button;
  }
  
  function createLoader(){
		var div = document.createElement("div");
		div.className = "loader";
	 return div;
  }
  
	function createLargePopup(title){
		document.getElementById('gauze').style.display = "block";
		document.getElementById('largepopuptitle').innerHTML = title;
		document.getElementById('sidemenutabs').innerHTML = "";
		document.getElementById('sidemenudisplay').innerHTML = "";
		document.getElementById('largepopupcontent').innerHTML = "";
		document.getElementById('largepopup').style.display = "block";
		var div = document.createElement("div");
		div.className = "loader";
	 return div;
  }
  
  function createThumbnails(thumbnails){
	  var div = document.createElement("div");
	  div.className = "thumbnails";
	  for(var i = 0; i < thumbnails.length ; i++){
		var thumbnail = document.createElement("div");
		thumbnail.className = "thumbnail";
		for(var j = 0; j < thumbnails[i].length; j++){
			var imagediv = document.createElement("div");
			imagediv.className = "thumbnailimagecontainer";
			var image = document.createElement("img");
			image.className = "thumbnailimage";
			image.alt = "Loading ...";
			var keys = Object.keys(thumbnails[i][j]);
			for(var k = 0; k < keys.length; k++){
				image.setAttribute("data-" + keys[k], thumbnails[i][j][keys[k]]);
			}
			imagediv.appendChild(image);
			thumbnail.appendChild(imagediv);
		}
		div.appendChild(thumbnail);
	  }
	 return div;
  }
  
  function showViewControls(){
	 var display = createSideMenuDisplay("View Controls");
	 var optiondiv = document.createElement("div");
	 optiondiv.className = "optiondiv";
	 if(targetflags.length > 0){
		var select = createSelect("Target", "viewcontrolstarget");
		for(var i = 0; i < targetflags.length; i++){
			select.getElementsByTagName('select')[0].appendChild(createOption(targetflags[i].id, i)); 
		}
		optiondiv.appendChild(select); 
		select.onchange = function(){setViewControls()};
	 }

	
	 optiondiv.appendChild(createInput("range", "Contrast", "viewcontrolscontrast", 0.5, 10, 0.5));
	 optiondiv.appendChild(createInput("range", "Brightness", "viewcontrolsbrightness", 0.5, 256, 0.5));
	 optiondiv.appendChild(createInput("range", "Zoom", "viewcontrolszoom", 1, 10, 0.1));
	 display.appendChild(optiondiv);
	 display.appendChild(createInput("range", "Scale", "viewcontrolsscale",0.2, 5, 0.1));
	 display.appendChild(createButton("Apply", updateViewControls));
	 document.body.appendChild(display); 
	 setViewControls();
  }
  
  function setViewControls(){
	 var viewcontrolstarget = document.getElementById('viewcontrolstarget');
	 var targetflagid = viewcontrolstarget.options[viewcontrolstarget.selectedIndex].value;
	 document.getElementById('viewcontrolscontrast').value = targetflags[targetflagid]['contrast'];
	 document.getElementById('viewcontrolsbrightness').value = targetflags[targetflagid]['brightness'];
	 document.getElementById('viewcontrolszoom').value = targetflags[targetflagid]['zoom'];
	 document.getElementById('viewcontrolsscale').value = targetscale;
  }
   
  function updateViewControls(){
	 var viewcontrolstarget = document.getElementById('viewcontrolstarget');
	 var targetflagid = viewcontrolstarget.options[viewcontrolstarget.selectedIndex].value;
	 targetflags[targetflagid]['contrast'] = document.getElementById('viewcontrolscontrast').value;
	 targetflags[targetflagid]['brightness'] = document.getElementById('viewcontrolsbrightness').value;
	 targetflags[targetflagid]['zoom'] = document.getElementById('viewcontrolszoom').value;
	 targetscale = document.getElementById('viewcontrolsscale').value;
	 setThumbnailSources();
  }
  
  function createProjectManagerProjects(data){ 
	  var projectsarray = data.projects;
	  var largepopupcontent = document.getElementById('largepopupcontent');
	  largepopupcontent.innerHTML = "";
	  for (var i = 0; i < projectsarray.length; i++){
		var div = document.createElement("div");
		var div1 = document.createElement("div");
		var div2 = document.createElement("div");
		var div3 = document.createElement("div");
		div1.innerHTML = projectsarray[i][0];
		div2.innerHTML = projectsarray[i][2];
		div.appendChild(div1);
		div.appendChild(div2);
		div.appendChild(div3);
		largepopupcontent.appendChild(div);
	}
  }
  
  function showProjectManager(){
	  hideMenu();
	  createLargePopup("Project Manager");
	  document.getElementById('largepopupcontent').appendChild(createLoader());
	  addSideMenuTab("add", showViewControls);
	  getAjax("handler?function=listprojects", createProjectManagerProjects);
  }
  
  function hidePopup(){
	  document.getElementById('gauze').style.display = "none";
	  document.getElementById('largepopup').style.display = "none";
  }
  
  function hideMenu(){
	  document.getElementById('headermenuitems').style.display = "none";
  }
  
  function createLoadQueue() {
	 loadqueue = [];
	 var images = document.getElementsByClassName('thumbnailimage');
	 var mainpane = document.getElementById('mainpane');
	 var paneBoundingClientRect = mainpane.getBoundingClientRect(); 
	 var topcutoff = paneBoundingClientRect.top;
	 var bottomcutoff = paneBoundingClientRect.bottom + (paneBoundingClientRect.bottom - paneBoundingClientRect.top);
	 for (var i = 0; i < images.length; i++){
		 var snapshotBoundingClientRect = images[i].getBoundingClientRect();
		 if(snapshotBoundingClientRect.bottom > topcutoff && snapshotBoundingClientRect.top < bottomcutoff) {
			 loadqueue.push(images[i]);
		 } else {
			 images[i].removeAttribute('src');
		 }
	 }
	 loadImages();
  }
  
  function loadImages(){
	if(loadqueue.length > 0){
		loadqueue[0].src = loadqueue[0].getAttribute('data-src');
		loadqueue[0].onload = loadqueue[0].onerror = function() { loadImages(); }
		loadqueue.shift();
	}
  }
  
  function setThumbnailSources() {
	  var images = document.getElementsByClassName('thumbnailimage');
	  for(var i = 0; i < images.length; i++){
		  var sourcefile = sourcefiles[images[i].getAttribute('data-sourcefile')];
		  var image = images[i].getAttribute('data-image');
		  var contrast = targetflags[images[i].getAttribute('data-targetflag')]["contrast"];
		  var brightness = targetflags[images[i].getAttribute('data-targetflag')]["brightness"];
		  var zoom = targetflags[images[i].getAttribute('data-targetflag')]["zoom"];
		  images[i].parentElement.style.width = 100 * targetscale + "px";
		  images[i].parentElement.style.height = 100 * targetscale + "px";
		  images[i].style.width = 100 * targetscale + "px";
		  images[i].style.height = 100 * targetscale + "px";
		  images[i].style.transform = "scale(" + zoom + ")";
		  images[i].setAttribute('data-src', "JPEGhandler?filename="+ rootdir + "/" + sourcefile + "&contrast=" + contrast + "&brightness=" + brightness + "&frameid=" + image);
	  }
	  createLoadQueue();
  }
  
  function createProjects(data){
	  var projectsarray = data.projects;
	  var projects = document.getElementById('projects');
	  projects.innerHTML = "";
	  for (var i = 0; i < projectsarray.length; i++){
		var div = document.createElement("div");
		div.innerHTML = projectsarray[i][0];
		div.id = projectsarray[i][0];
		div.setAttribute("data-projectdescription", projectsarray[i][1]);
		div.setAttribute("data-projectfolder", projectsarray[i][2]);
		div.className = "project";
		div.onclick = function() { selectProject(this) };
		projects.appendChild(div);
	}
  }
  
  function selectProject(projectdiv){
	  document.getElementById('projects').style.display = "none";
	  currentproject = projectdiv.id;
	  document.getElementById('projectselectortitle').innerHTML = currentproject;  
  }
  
  function getProjects(){
	  var projects = document.getElementById('projects');
	  projects.innerHTML = "";
	  projects.appendChild(createLoader());
	  projects.style.display = "block";
	  getAjax("handler?function=listprojects", createProjects);
	  
  }
  
  function getMenu() {
	  document.getElementById('headermenuitems').style.display = "block";
  }
  
  function ini3dViewInit(data){
	var JSONdata = JSON.parse(data);
	sourcefiles = JSONdata.sourcefiles;
	targetflags = JSONdata.targetflags;
	targetscale = 1;
	rootdir = JSONdata.rootdir;
	document.getElementById('mainpane').innerHTML = "";
	document.getElementById('bodytitle').innerHTML = "Ini3D View";
	addSideMenuTab("View", showViewControls);
	document.getElementById('mainpane').appendChild(createThumbnails(JSONdata.thumbnails));
	setThumbnailSources();
	document.getElementById('mainpane').addEventListener('scroll', function(){createLoadQueue()});
  }
  
//  getAjax("JSONhandler?function=ini3dview&" + window.location.toString().split('?')[1], ini3dViewInit);
