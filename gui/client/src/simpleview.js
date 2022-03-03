class SimpleView {

  constructor(){}
  
  showLoader(){
	  var loader = document.createElement('div')
	  loader.className = "loader"
	  var mainpane = document.getElementById('mainpane')
	  mainpane.innerHTML = ""
	  mainpane.appendChild(loader) 
  }
  
  toggleViewtype(element){
	  for(var viewtype of document.getElementsByClassName('viewtype')){
		  viewtype.dataset.selected = false
	  }
	  element.dataset.selected = true
  }
  
  updateTitle(viewtype){
	  document.getElementById('simpleviewtitle').innerHTML = "View - " + viewtype
  }
  
  resetSideBar(){
	  document.getElementById('sidebar').innerHTML = ""
  }
  
  getRaw(projfile){
	  this.projfile = projfile
	  this.updateTitle('Raw')
	  this.resetSideBar()
	  this.showLoader()
	  var request = {
		  cls : "view",
		  fnc : "getViewSimpleRaw",
		  arg : {
			  projfile : projfile
			  }
		}
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			document.getElementById('mainpane').innerHTML = json.html
		})  
  }
  
  getManualPick(projfile){

	  this.projfile = projfile
	  this.updateTitle('Manual Particle Picking')
	  this.resetSideBar()
	  this.showLoader()
	  var request = {
		  cls : "view",
		  fnc : "getViewSimpleManualPick",
		  arg : {
			  projfile : projfile
			  }
		}
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			document.getElementById('mainpane').innerHTML = json.html
			this.addImageControls(['picking'], true, false)
			this.addBoxControls()
			this.addSaveControls('manualpick')
			if(json.boxes != undefined){
				this.boxes = json.boxes['boxes']
				document.getElementById('boxsize').value = json.boxes['boxsize']
				document.getElementById('ptclradius').value = json.boxes['ptclradius']
				document.getElementById('boxstep').value = json.boxes['boxstep']
			}else{
				this.boxes = {}
			}
			if(json.boxes != undefined && json.boxes['tracks'] != undefined){
				this.tracks = json.boxes['tracks']
			}else{
				this.tracks = {}
			}
			document.getElementById("pickingcanvas").addEventListener("click", (evt) => {
					var helical = document.getElementsByName('boxmode')
					if(!helical[0].checked){
						if(!this.boxes[this.currentmic]){
							  this.boxes[this.currentmic] = []
						}
						var canvas = document.getElementById("pickingcanvas")
						var rect = canvas.getBoundingClientRect();
						var x = evt.clientX - rect.left
						var y = evt.clientY - rect.top
						var xscale = canvas.width/rect.width
						var yscale = canvas.height/rect.height
						var xcoord = Math.round(x * xscale)
						var ycoord = Math.round(y * yscale)
						var newbox = true
						for(var box = this.boxes[this.currentmic].length - 1; box >= 0; box--){
							var distance = Math.sqrt(Math.pow(xcoord - this.boxes[this.currentmic][box][0] ,2) + Math.pow(ycoord - this.boxes[this.currentmic][box][1] ,2))
							if(distance < document.getElementById('ptclradius').value / 2){
								this.boxes[this.currentmic].splice(box,1)
								newbox = false
							}
						}
						if(newbox){
							this.boxes[this.currentmic].push([xcoord, ycoord])
						}
						this.drawBoxes()
					}
			})
			this.getManualPickMicrograph(document.getElementById('iterationselector'))
		})
  }
  

  manualPickHelix(){
	this.boxes[this.currentmic] = []
	this.findBoxes()
	document.getElementById("pickingcanvas").onmousedown = function() {
		if(!simpleview.boxes[simpleview.currentmic]){
				simpleview.boxes[simpleview.currentmic] = []
		}
		var ptclradius = document.getElementById('ptclradius').value
		var canvas = document.getElementById("pickingcanvas")
		var rect = canvas.getBoundingClientRect();
		var x = event.clientX - rect.left
		var y = event.clientY - rect.top
		var xscale = canvas.width/rect.width
		var yscale = canvas.height/rect.height
		var xcoord = Math.round(x * xscale)
		var ycoord = Math.round(y * yscale)	
		
		if(!simpleview.tracks[simpleview.currentmic]){
			simpleview.tracks[simpleview.currentmic] =  []
		}
		
		document.getElementById("pickingcanvas").onmouseup = function(){
			simpleview.mouseUp({"x":xcoord,"y":ycoord}, false)
		}
		
		var newpoint = true
		
		for ( var i = 0 ; i < simpleview.tracks[simpleview.currentmic].length ; i++ ){
			for ( var j = 0; j < simpleview.tracks[simpleview.currentmic][i].length ; j++ ){
				if ( (xcoord-simpleview.tracks[simpleview.currentmic][i][j].xcoord)**2 + (ycoord-simpleview.tracks[simpleview.currentmic][i][j].ycoord)**2 < ptclradius**2 ){
					var point = {"index":i,"pointno":j}
					document.getElementById("pickingcanvas").onmouseup = function(){
						simpleview.mouseUp(point, false)
					}
					document.getElementById("pickingcanvas").onmousemove = function(event){
						simpleview.checkMousePos(point)
					}
					
					break
				}
			}
		}
	};
	
  }

 checkMousePos(point){

	var ptclradius = document.getElementById('ptclradius').value
	var canvas = document.getElementById("pickingcanvas")
	var rect = canvas.getBoundingClientRect();
	var x = event.clientX - rect.left
	var y = event.clientY - rect.top
	var xscale = canvas.width/rect.width
	var yscale = canvas.height/rect.height
	var xcoord = Math.round(x * xscale)
	var ycoord = Math.round(y * yscale)	
	if ( (xcoord - this.tracks[this.currentmic][point.index][point.pointno].xcoord)**2 + (ycoord - this.tracks[this.currentmic][point.index][point.pointno].ycoord)**2 < ptclradius**2 ){
		document.getElementById("pickingcanvas").onmouseup = function(){
			simpleview.mouseUp(point, true)
		}
		document.getElementById("pickingcanvas").onmousemove = function(){
			simpleview.changeCoords(point)
		}
	}
  }
  
  changeCoords(point){
//	var x = event.offsetX;
//	var y = event.offsetY;
	var canvas = document.getElementById("pickingcanvas")
	var rect = canvas.getBoundingClientRect();
	var x = event.clientX - rect.left
	var y = event.clientY - rect.top
	var xscale = canvas.width/rect.width
	var yscale = canvas.height/rect.height
	var xcoord = Math.round(x * xscale)
	var ycoord = Math.round(y * yscale)	
	
	this.tracks[this.currentmic][point.index][point.pointno].xcoord = xcoord
	this.tracks[this.currentmic][point.index][point.pointno].ycoord = ycoord
	this.findBoxes()
	this.drawBoxes()
}

  mouseUp(point, dragging){
	document.getElementById("pickingcanvas").onmousemove = function(){}

	// clicked on a point, so either deleting or moving it
	if ( point.hasOwnProperty("index") ){
		// clicked on a point and didn't drag
		if ( !dragging ) {
			this.tracks[this.currentmic].splice(point.index,1)
			this.findBoxes()
			this.drawBoxes()
		}
	} else {
		// new point
		if ( this.tracks[this.currentmic].length !== 0 && this.tracks[this.currentmic][this.tracks[this.currentmic].length - 1].length == 1 ) {
			// if previous coord is unpaired, pair it, otherwise add on its own
			this.tracks[this.currentmic][this.tracks[this.currentmic].length - 1].push({"xcoord":point.x, "ycoord":point.y})
			this.findBoxes()
		} else {
			this.tracks[this.currentmic].push([{"xcoord":point.x, "ycoord":point.y}])
		}
		
		this.drawBoxes()
	}
}

  findBoxes(){
	var helical = document.getElementsByName('boxmode')
	if(helical[0].checked){
		this.boxes[this.currentmic] = []
		if(!this.tracks[this.currentmic]){
			this.tracks[this.currentmic] = []
		}
		
		var step = Number(document.getElementById("boxstep").value);
		
		for ( var coords of this.tracks[this.currentmic] ){
			if ( coords.length == 2 ){
				var linelength = Math.sqrt( (coords[0].xcoord-coords[1].xcoord)**2 + (coords[0].ycoord-coords[1].ycoord)**2 )
				
				var num = Math.round(linelength/step);
				
				var xseg = (coords[1].xcoord-coords[0].xcoord)/num
				var yseg = (coords[1].ycoord-coords[0].ycoord)/num
				
				for (var box = 0; box <= num; box++){
					var xc = coords[0].xcoord + (box*xseg)
					var yc = coords[0].ycoord + (box*yseg)
					this.boxes[this.currentmic].push([xc, yc])
				}
			}
		}
	}
  }
  
  findAllBoxes(){
	var helical = document.getElementsByName('boxmode')
	if(helical[0].checked){
		for(var key of Object.keys(this.boxes)){
			this.boxes[key] = []
			if(!this.tracks[key]){
				this.tracks[key] = []
			}
		
			var step = Number(document.getElementById("boxstep").value);
		
				for ( var coords of this.tracks[key] ){
				if ( coords.length == 2 ){
					var linelength = Math.sqrt( (coords[0].xcoord-coords[1].xcoord)**2 + (coords[0].ycoord-coords[1].ycoord)**2 )
					
					var num = Math.round(linelength/step);
					
					var xseg = (coords[1].xcoord-coords[0].xcoord)/num
					var yseg = (coords[1].ycoord-coords[0].ycoord)/num
					
					for (var box = 0; box <= num; box++){
						var xc = coords[0].xcoord + (box*xseg)
						var yc = coords[0].ycoord + (box*yseg)
						this.boxes[key].push([xc, yc])
					}
				}
			}
		}
	}
  } 
  

  getAutoPick(projfile){
	  this.projfile = projfile
	  this.updateTitle('Automatic Particle Picking')
	  this.resetSideBar()
	  this.showLoader()
	  var request = {
		  cls : "view",
		  fnc : "getViewSimplePick",
		  arg : {
			  projfile : projfile
			  }
		}
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			document.getElementById('mainpane').innerHTML = json.html
			this.addImageControls(['picking'], true, false)
			this.addBoxControls()
			this.addSaveControls('manualpick')
			if(json.boxes != undefined){
				this.boxes = json.boxes['boxes']
				document.getElementById('boxsize').value = json.boxes['boxsize']
				document.getElementById('ptclradius').value = json.boxes['ptclradius']
			}else{
				this.boxes = {}
			}
			this.getManualPickMicrograph(document.getElementById('iterationselector'))
		})
  }
  
  getManualPickMicrograph(selector, next, previous){
	  var selectedindex = selector.selectedIndex
	  if(next != undefined && next == true && selectedindex < selector.length - 1){
		  selector.selectedIndex = selectedindex + 1
	  }else if(previous != undefined && previous == true && selectedindex > 0){
		  selector.selectedIndex = selectedindex - 1
	  }
	  var intg = selector.options[selector.selectedIndex].dataset.intg
	  var xdim = selector.options[selector.selectedIndex].dataset.xdim
	  var ydim = selector.options[selector.selectedIndex].dataset.ydim
	  var canvas = document.getElementById("pickingcanvas")
	  var ctx = canvas.getContext("2d")
	  canvas.width = xdim
      canvas.height = ydim
      var rect = document.getElementById('thumbnails').getBoundingClientRect()
      var mindim = (rect.width < rect.height) ? rect.width : rect.height
      //canvas.style.width = mindim + 'px'
      canvas.style.height = mindim + 'px'
      this.currentmic = selector.options[selector.selectedIndex].value
	  this.background = new Image()
	  var url = window.location.href + "/image?stackfile=" + intg + "&frame=0&width=" + xdim + "&pick=true"
	  this.background.src = url.replace('//image','/image')
	  if(!this.boxes[this.currentmic]){
		this.boxes[this.currentmic] = []
	  }
	  var boxtotal = 0
	  for(var key in this.boxes){
			if(key != 'boxsize' && key != 'ptclradius'){
				boxtotal += this.boxes[key].length
			}
	  }  
	  
	  document.getElementById("particlecounter").innerHTML = 'Boxes ' + this.boxes[this.currentmic].length + '/' + boxtotal
	  this.background.onload = function(){
		  simpleview.drawBoxes()
		//  simpleview.drawTracks()
		  
	  }
  }

  drawBoxes(){
		var boxtotal = 0
		for(var key in this.boxes){
		if(key != 'boxsize' && key != 'ptclradius')
			boxtotal += this.boxes[key].length
		}   
				
	  document.getElementById("particlecounter").innerHTML = 'Boxes ' + this.boxes[this.currentmic].length + '/' + boxtotal
	  
	  var canvas = document.getElementById("pickingcanvas")
	  var ctx = canvas.getContext("2d")
	  ctx.clearRect(0, 0, canvas.width, canvas.height)
	  ctx.drawImage(this.background,0,0)
	  var boxsize = document.getElementById('boxsize').value
	  var ptclradius = document.getElementById('ptclradius').value
	  if(this.boxes[this.currentmic]){
		  for(var box of this.boxes[this.currentmic]){
				  ctx.beginPath()
				  ctx.lineWidth = document.getElementById('linewidth').value
				  ctx.strokeStyle = 'green'
				  if(document.getElementById('showradii').checked){
					ctx.arc(box[0], box[1], ptclradius, 0, 2 * Math.PI)
				  }
				  if(document.getElementById('showboxes').checked){
					ctx.rect(box[0] - boxsize/2, box[1] - boxsize/2, boxsize, boxsize)
				  }
				  ctx.stroke()
			 }
		 }
	  var helical = document.getElementsByName('boxmode')
	  if(helical[0].checked){
		  if(this.tracks[this.currentmic]){
			for ( var coords of this.tracks[this.currentmic] ){
				ctx.beginPath();
				ctx.arc(coords[0].xcoord,coords[0].ycoord, ptclradius,0,2*Math.PI)
				ctx.lineWidth = document.getElementById('linewidth').value
				ctx.strokeStyle="DodgerBlue"
				ctx.stroke()

				if ( coords.length == 2 ){
					ctx.beginPath()
					ctx.arc(coords[1].xcoord,coords[1].ycoord, ptclradius,0,2*Math.PI)
					ctx.stroke()
					
					ctx.beginPath()
					ctx.moveTo(coords[0].xcoord, coords[0].ycoord)
					ctx.lineTo(coords[1].xcoord, coords[1].ycoord)
					ctx.strokeStyle="Red"
					ctx.stroke()
				}
			}
		}
	}
  }
  
  saveBoxes(){
	  document.getElementById('savebutton').innerHTML = 'Saving'
	  var request = {
		  cls : "view",
		  fnc : "saveManualPickBoxes",
		  arg : {
			  projfile : this.projfile,
			  boxes: this.boxes,
			  tracks: this.tracks,
			  boxsize:document.getElementById('boxsize').value,
			  ptclradius:document.getElementById('ptclradius').value,
			  boxstep:document.getElementById('boxstep').value
			  }
	  }
	fetcher.fetchJSON(request)
	.then(response => response.json())
	.then (json => {
		 alert("Boxes Saved")
		 document.getElementById('savebutton').innerHTML = 'Save'
	})
	.catch(() => {
		 alert("Error Saving Boxes")
		 document.getElementById('savebutton').innerHTML = 'Save'
	})
  }
  
  clearBoxes(){
	 this.boxes = {} 
	 this.drawBoxes()
  }
  
  getMicrographs(projfile){
	  fetcher.killWorker()
	  this.projfile = projfile
	  this.updateTitle('Micrographs')
	  this.resetSideBar()
	  this.showLoader()
	  var request = {
		  cls : "view",
		  fnc : "getViewSimpleMicrographs",
		  arg : {
			  projfile : projfile
			  }
	  }
	 fetcher.fetchJSON(request)
	 .then(response => response.json())
		.then (json => {
			var viewmaster = document.createElement('div')
			var micrographcount = document.createElement('h')
			micrographcount.id = "micrographcount"
			micrographcount.className = "count"
			viewmaster.appendChild(micrographcount)
			var particlecount = document.createElement('h')
			particlecount.id = "particlecount"
			particlecount.className = "count"
			viewmaster.appendChild(particlecount)
			var thumbnails = document.createElement('div')
			thumbnails.id = "thumbnails"
			var trackviewer = document.createElement('div')
			var trackviewergauze = document.createElement('div')
			var trackviewerimg= document.createElement('img')
			trackviewergauze.className = "gauze"
			trackviewerimg.id = "trackviewerimg"
			trackviewer.appendChild(trackviewergauze)
			trackviewer.appendChild(trackviewerimg)
			trackviewer.id = "trackviewer"
			
			this.resetSideBar()
			var targets = []
			for(var micrograph of json.stats){
				var thumbnail = document.createElement('div')
				thumbnail.className = "thumbnail"
				thumbnail.dataset.number = micrograph['number']
				
				var container = document.createElement('div')
				var img = document.createElement('img')
				img.src = "img/log.png"
				container.className = "thumbnailbutton"
				container.dataset.path = window.location.href + "/image?trackfile=" + micrograph.intg + "&stackfile=.jpg"
				container.setAttribute("onclick", "event.stopPropagation(); simpleview.viewTracks(this)")
				container.title = "View motion tracks"
				img.className = "motiontracksicon"
				container.appendChild(img)
				
				thumbnail.appendChild(container)
				
				if(!micrograph.state){
					thumbnail.dataset.state = 1
				}
				thumbnail.setAttribute("onclick", "simpleview.toggleThumbnail(this)")
				thumbnails.appendChild(thumbnail)
				
				
				if(micrograph.thumb){
					targets.push("micrograph")
					targets.push("pspec")
					var container = document.createElement('div')
					container.className = "thumbnailcontainer"
					container.style.width = "200px"
					container.style.height = "200px"
					var containerloader = document.createElement('div')
					containerloader.className = "miniloader"
					container.appendChild(containerloader)
					var containerimg = document.createElement('img')
					containerimg.className = "thumbnailcontainerimg pspec"
					containerimg.dataset.path = window.location.href + "/image?stackfile=" + micrograph.thumb + "&spriteframe=0"
					container.appendChild(containerimg)
					thumbnail.appendChild(container)
					var container = document.createElement('div')
					container.className = "thumbnailcontainer"
					container.style.width = "200px"
					container.style.height = "200px"
					var containerloader = document.createElement('div')
					containerloader.className = "miniloader"
					container.appendChild(containerloader)
					var containerimg = document.createElement('img')
					containerimg.className = "thumbnailcontainerimg micrograph"
					containerimg.dataset.path = window.location.href + "/image?stackfile=" + micrograph.thumb + "&spriteframe=1"
					container.appendChild(containerimg)
					thumbnail.appendChild(container)
				}else{
					targets.push("micrograph")
					var container = document.createElement('div')
					container.className = "thumbnailcontainer"
					container.style.width = "200px"
					container.style.height = "200px"
					var containerloader = document.createElement('div')
					containerloader.className = "miniloader"
					container.appendChild(containerloader)
					var containerimg = document.createElement('img')
					containerimg.className = "thumbnailcontainerimg micrograph"
					containerimg.dataset.path = window.location.href + "/image?stackfile=" + micrograph.intg + "&frame=0"
					container.appendChild(containerimg)
					thumbnail.appendChild(container)
				}
				if(micrograph.ctfjpg){
					targets.push("ctffit")
					var container = document.createElement('div')
					container.className = "thumbnailcontainer"
					container.style.width = "200px"
					container.style.height = "200px"
					var containerloader = document.createElement('div')
					containerloader.className = "miniloader"
					container.appendChild(containerloader)
					var containerimg = document.createElement('img')
					containerimg.className = "thumbnailcontainerimg ctffit"
					containerimg.dataset.path = window.location.href + "/image?stackfile=" + micrograph.ctfjpg + "&frame=0"
					container.appendChild(containerimg)
					thumbnail.appendChild(container)
				}
				for(var key of Object.keys(micrograph)){
					thumbnail.title += key + " : " + micrograph[key] + "\n"
					thumbnail.dataset[key] = Number(micrograph[key])
				}
				if(micrograph.boxfile){
					targets.push("boxes")
					var container = document.createElement('div')
					container.className = "thumbnailcontainer"
					container.style.width = "200px"
					container.style.height = "200px"
					var containerloader = document.createElement('div')
					containerloader.className = "miniloader"
					container.appendChild(containerloader)
					var containerimg = document.createElement('img')
					containerimg.className = "thumbnailcontainerimg boxes"
					containerimg.dataset.path = window.location.href + "/image?stackfile=" + micrograph.thumb + "&frame=0&boxfile=" + micrograph.boxfile + "&intg=" + micrograph.intg
					container.appendChild(containerimg)
					thumbnail.appendChild(container)
				}
				
			}
			var mainpane = document.getElementById('mainpane')
			mainpane.innerHTML = ""
			mainpane.appendChild(viewmaster)
			mainpane.appendChild(thumbnails)
			trackviewer.style.display = "none"
			mainpane.appendChild(trackviewer)
			this.addImageControls(Array.from(new Set(targets)), true, false)
			this.addSortControls(json.stats)
			this.addSelectControls(json.stats)
			this.addPlotControls(json.stats)
			this.addSaveControls("mic")  
			this.updateMicrographCount()
			this.updateParticleCount()
			fetcher.imagewidth = 200
			fetcher.fetchImages()
		})

  }
  
  viewTracks(element){
	  var trackviewer = document.getElementById("trackviewer")
	  var trackviewerimg = document.getElementById("trackviewerimg")
	  trackviewerimg.src = ""
	  
	  trackviewerimg.src = element.dataset.path.replace("//image","/image")
	  trackviewer.style.display =  "block"
	  trackviewer.setAttribute("onclick", "document.getElementById('trackviewer').style.display = 'none'")
  }
  
  
  
  getIni3D(projfile){
	  fetcher.killWorker()
	  this.projfile = projfile
	  this.updateTitle('Initial 3D Model')
	  this.resetSideBar()
	  this.showLoader()
	  var request = {
		  cls : "view",
		  fnc : "getViewSimpleIni3D",
		  arg : {
			  projfile : projfile
			  }
	  }
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			document.getElementById('mainpane').innerHTML = json.html
			var iterationselector = document.getElementById('iterationselector')
			iterationselector.selectedIndex = iterationselector.options.length -1
			simpleview.getIni3DIteration(iterationselector)
		})  
  }
  
  getIni3DIteration(element){
	  fetcher.killWorker()
	  var selected = element.options[element.selectedIndex]
	  var request = {
		  cls : "view",
		  fnc : "prepareMRCVolume",
		  arg : {
			  volume : selected.dataset['file']
		   }
	  }
	  fetcher.fetchJSON(request)
	  .then(response => response.json())
	  .then (json => {
		  this.resetSideBar()
		  if(this.plugin){
			  this.plugin.destroy()
		  }
		  if(selected.dataset.fsc){
			var fsc = JSON.parse(selected.dataset.fsc)
			var plotimg = document.createElement("img")
			plotimg.src = "img/plot.png"
			plotimg.title = "Plot FSC Curve"
			plotimg.dataset.fsc = selected.dataset.fsc
			plotimg.setAttribute("onclick", "simpleview.plotFSC(this)")
			document.getElementById('resolution').innerHTML = "Estimated Resolution (FSC=0.500/FSC=0.143) : " + fsc['fsc500'] + " / " + fsc['fsc143']
			document.getElementById('resolution').appendChild(plotimg)
		  }else{
			 document.getElementById('resolution').innerHTML = "" 
		  }
		  if(selected.dataset.stats){
			  document.getElementById('reprojections').style.display = "block"
			  var thumbnails = document.getElementById('thumbnails')
			  thumbnails.style.display = "block"
			  thumbnails.innerHTML = ""
			  var header = JSON.parse(selected.dataset.header)
			  var stats = JSON.parse(selected.dataset.stats)
			  var xdim = (header.nx < 100) ? header.nx : 100;
			  var ydim = (header.ny < 100) ? header.ny : 100;
			  for(var thumbcount = 0; thumbcount < stats.length; thumbcount++){
				  var thumbnail = document.createElement('div')
				  thumbnail.className = "thumbnail"
				  var container = document.createElement('div')
				  container.className = "thumbnailcontainer"
				  container.style.width = xdim + "px"
				  container.style.height = ydim + "px"
				  var containerloader = document.createElement('div')
				  containerloader.className = "miniloader"
				  container.appendChild(containerloader)
				  var containerimg = document.createElement('img')
				  containerimg.className = "thumbnailcontainerimg cls2D"
				  containerimg.dataset.path = window.location.href + "/image?stackfile=" + selected.dataset.projections + "&frame=" + Number(thumbcount * 2)
				  container.appendChild(containerimg)
				  thumbnail.appendChild(container)
				  var container = document.createElement('div')
				  container.className = "thumbnailcontainer"
				  container.style.width = xdim + "px"
				  container.style.height = ydim + "px"
				  var containerloader = document.createElement('div')
				  containerloader.className = "miniloader"
				  container.appendChild(containerloader)
				  var containerimg = document.createElement('img')
				  containerimg.className = "thumbnailcontainerimg cls2D"
				  containerimg.dataset.path = window.location.href + "/image?stackfile=" + selected.dataset.projections + "&frame=" + Number((thumbcount * 2) + 1)
				  container.appendChild(containerimg)
				  thumbnail.appendChild(container)
				  if(stats){
					  for(var key of Object.keys(stats[thumbcount])){
						  thumbnail.title += key + " : " + stats[thumbcount][key] + "\n"
						  thumbnail.dataset[key] = Number(stats[thumbcount][key])
					  }
					 thumbnail.setAttribute("onclick", "simpleview.toggleThumbnail(this)")
				  }
				  thumbnails.appendChild(thumbnail)
			  }
			  this.addImageControls(['cls2D'], true, true)
			  this.addPlotControls(stats)
			  this.addSaveControls('cls2D')
			  var controls = this.createControlsDiv("fsccontrolpane")
			  controls.appendChild(this.createControlsHeader("FSC"))
			  var body = document.createElement('div')
			  var div = document.createElement('div')
			  div.id = "fscchart"
			  body.appendChild(div)
			  controls.appendChild(body)
			  document.getElementById('mainpane').appendChild(controls)
			  
			  fetcher.imagewidth = (header.nx < 100) ? header.nx : 100
			  fetcher.fetchImages()
		  }else{
			  this.addImageControls([], false, true)
			  document.getElementById('reprojections').style.display = "none"
			  document.getElementById('thumbnails').style.display = "none"
			  document.getElementById('resolution').innerHTML =""
		  }
		  this.plugin = LiteMol.Plugin.create({ target: '#litemol' })
		  var detail = document.getElementById('quality').value
		  this.mdb = json['mdb']
		  var url = window.location.href + "/DensityServer/local/"
		   console.log('HREF', url)
		  var surface = this.plugin.createTransform()
		  .add(this.plugin.root, LiteMol.Bootstrap.Entity.Transformer.Data.Download, { url: url.replace('//DensityServer', '/DensityServer') + json['mdb'] + "/cell?detail=" + detail, type: 'Binary', description: 'local Density', title: "joe2", id:"density3ddata"})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Data.ParseBinaryCif, {id:"density3dcif"}, {})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Density.CreateFromCif, {id:"density3dd", blockIndex: 1})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Density.CreateVisual, {
			style: LiteMol.Bootstrap.Visualization.Density.Style.create({
			  isoValue: 0,
			  isoValueType: LiteMol.Bootstrap.Visualization.Density.IsoValueType.Sigma,
			  color: LiteMol.Visualization.Color.fromHex(0xBB3333),
			  isWireframe: false,
			  transparency: { alpha: 1.0 }
			  })
			}, {ref : "density3d"})
		  this.plugin.applyTransform(surface);
		  LiteMol.Bootstrap.Command.Layout.SetViewportOptions.dispatch(this.plugin.context, { clearColor: LiteMol.Visualization.Color.fromRgb(255, 255, 255) });
		  LiteMol.Bootstrap.Command.Visual.ResetScene.dispatch(this.plugin.context, void 0)
	  })
  }
  
  getRefine3D(projfile){
	  fetcher.killWorker()
	  this.projfile = projfile
	  this.updateTitle('3D Volume')
	  this.resetSideBar()
	  this.showLoader()
	  var request = {
		  cls : "view",
		  fnc : "getViewSimpleRefine3D",
		  arg : {
			  projfile : projfile
			  }
	  }
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			document.getElementById('mainpane').innerHTML = json.html
			var iterationselector = document.getElementById('iterationselector')
			iterationselector.selectedIndex = iterationselector.options.length -1
			simpleview.getRefine3DIteration(iterationselector)
		})  
  }
  
  getRefine3DIteration(element){
	  fetcher.killWorker()
	  var selected = element.options[element.selectedIndex]
	  var request = {
		  cls : "view",
		  fnc : "prepareMRCVolume",
		  arg : {
			  volume : selected.dataset['file']
		   }
	  }
	  fetcher.fetchJSON(request)
	  .then(response => response.json())
	  .then (json => {
		  this.resetSideBar()
		  if(this.plugin){
			  this.plugin.destroy()
		  }
		  if(selected.dataset.fsc){
			var fsc = JSON.parse(selected.dataset.fsc)
			var plotimg = document.createElement("img")
			plotimg.src = "img/plot.png"
			plotimg.title = "Plot FSC Curve"
			plotimg.dataset.fsc = selected.dataset.fsc
			plotimg.setAttribute("onclick", "simpleview.plotFSC(this)")
			document.getElementById('resolution').innerHTML = "Estimated Resolution (FSC=0.500/FSC=0.143) : " + fsc['fsc500'] + " / " + fsc['fsc143']
			document.getElementById('resolution').appendChild(plotimg)
			this.addImageControls([], false, true)
			var controls = this.createControlsDiv("fsccontrolpane")
			controls.appendChild(this.createControlsHeader("FSC"))
			var body = document.createElement('div')
			var div = document.createElement('div')
			div.id = "fscchart"
			body.appendChild(div)
			controls.appendChild(body)
			document.getElementById('mainpane').appendChild(controls)
		  }else{
			 document.getElementById('resolution').innerHTML = "" 
		  }
		  this.plugin = LiteMol.Plugin.create({ target: '#litemol' })
		  var detail = document.getElementById('quality').value
		  this.mdb = json['mdb']
		  var url = window.location.href + "/DensityServer/local/"
		  var surface = this.plugin.createTransform()
		  .add(this.plugin.root, LiteMol.Bootstrap.Entity.Transformer.Data.Download, { url: url.replace('//DensityServer', '/DensityServer') + json['mdb'] + "/cell?detail=" + detail, type: 'Binary', description: 'local Density', title: "joe2", id:"density3ddata"})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Data.ParseBinaryCif, {id:"density3dcif"}, {})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Density.CreateFromCif, {id:"density3dd", blockIndex: 1})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Density.CreateVisual, {
			style: LiteMol.Bootstrap.Visualization.Density.Style.create({
			  isoValue: 0,
			  isoValueType: LiteMol.Bootstrap.Visualization.Density.IsoValueType.Sigma,
			  color: LiteMol.Visualization.Color.fromHex(0xBB3333),
			  isWireframe: false,
			  transparency: { alpha: 1.0 }
			  })
			}, {ref : "density3d"})
		  this.plugin.applyTransform(surface);
		  LiteMol.Bootstrap.Command.Layout.SetViewportOptions.dispatch(this.plugin.context, { clearColor: LiteMol.Visualization.Color.fromRgb(255, 255, 255) });
		  LiteMol.Bootstrap.Command.Visual.ResetScene.dispatch(this.plugin.context, void 0)
	  })
  }
  
  getPostprocess(projfile){
	  fetcher.killWorker()
	  this.projfile = projfile
	  this.updateTitle('3D Volume')
	  this.resetSideBar()
	  this.showLoader()
	  var request = {
		  cls : "view",
		  fnc : "getViewSimplePostprocess",
		  arg : {
			  projfile : projfile
			  }
	  }
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			document.getElementById('mainpane').innerHTML = json.html
			var iterationselector = document.getElementById('iterationselector')
			iterationselector.selectedIndex = iterationselector.options.length -1
			simpleview.getPostprocessVolume(iterationselector)
		})  
  }
  
  getPostprocessVolume(element){
	  fetcher.killWorker()
	  var selected = element.options[element.selectedIndex]
	  var request = {
		  cls : "view",
		  fnc : "prepareMRCVolume",
		  arg : {
			  volume : selected.dataset['file']
		   }
	  }
	  fetcher.fetchJSON(request)
	  .then(response => response.json())
	  .then (json => {
		  this.resetSideBar()
		  if(this.plugin){
			  this.plugin.destroy()
		  }
		  this.addImageControls([], false, true)
		  this.plugin = LiteMol.Plugin.create({ target: '#litemol' })
		  var detail = document.getElementById('quality').value
		  this.mdb = json['mdb']
		  var url = window.location.href + "/DensityServer/local/"
		  var surface = this.plugin.createTransform()
		  .add(this.plugin.root, LiteMol.Bootstrap.Entity.Transformer.Data.Download, { url: url.replace('//DensityServer', '/DensityServer') + json['mdb'] + "/cell?detail=" + detail, type: 'Binary', description: 'local Density', title: "joe2", id:"density3ddata"})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Data.ParseBinaryCif, {id:"density3dcif"}, {})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Density.CreateFromCif, {id:"density3dd", blockIndex: 1})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Density.CreateVisual, {
			style: LiteMol.Bootstrap.Visualization.Density.Style.create({
			  isoValue: 0,
			  isoValueType: LiteMol.Bootstrap.Visualization.Density.IsoValueType.Sigma,
			  color: LiteMol.Visualization.Color.fromHex(0xBB3333),
			  isWireframe: false,
			  transparency: { alpha: 1.0 }
			  })
			}, {ref : "density3d"})
		  this.plugin.applyTransform(surface);
		  LiteMol.Bootstrap.Command.Layout.SetViewportOptions.dispatch(this.plugin.context, { clearColor: LiteMol.Visualization.Color.fromRgb(255, 255, 255) });
		  LiteMol.Bootstrap.Command.Visual.ResetScene.dispatch(this.plugin.context, void 0)
	  })
  }
  
  updateIsoSurface(element) {
    var density = this.plugin.context.select('density3d')[0]
    var style = LiteMol.Bootstrap.Visualization.Density.Style.create({
      isoValue: Number(element.value),
      isoValueType: LiteMol.Bootstrap.Visualization.Density.IsoValueType.Sigma,
      color: LiteMol.Visualization.Color.fromHex(0xBB3333),
      isWireframe: false,
      transparency: { alpha: 1.0 }
    })
    var absolute = density['parent']['props']['data']['valuesInfo'].mean + density['parent']['props']['data']['valuesInfo'].sigma * Number(element.value)
    document.getElementById('absolute').value = absolute
    LiteMol.Bootstrap.Entity.Transformer.Density.CreateVisual.create({ style: style }, { ref: density.ref, isHidden: false }).update(this.plugin.context, density).run()
  }
  
  reloadIsoSurface() {
		  var detail = document.getElementById('quality').value
		  this.plugin.clear()
		  var url = window.location.href + "/DensityServer/local/"
		  var surface = this.plugin.createTransform()
		  .add(this.plugin.root, LiteMol.Bootstrap.Entity.Transformer.Data.Download, { url: url.replace('//DensityServer', '/DensityServer') + this.mdb + "/cell?detail=" + detail, type: 'Binary', description: 'local Density', title: "joe2", id:"density3ddata"})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Data.ParseBinaryCif, {id:"density3dcif"}, {})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Density.CreateFromCif, {id:"density3dd", blockIndex: 1})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Density.CreateVisual, {
			style: LiteMol.Bootstrap.Visualization.Density.Style.create({
			  isoValue: 0,
			  isoValueType: LiteMol.Bootstrap.Visualization.Density.IsoValueType.Sigma,
			  color: LiteMol.Visualization.Color.fromHex(0xBB3333),
			  isWireframe: false,
			  transparency: { alpha: 1.0 }
			  })
			}, {ref : "density3d"})
		  this.plugin.applyTransform(surface);
		  LiteMol.Bootstrap.Command.Layout.SetViewportOptions.dispatch(this.plugin.context, { clearColor: LiteMol.Visualization.Color.fromRgb(255, 255, 255) });
		  LiteMol.Bootstrap.Command.Visual.ResetScene.dispatch(this.plugin.context, void 0)
  }
  
  plotFSC(element){
	  var fsc = JSON.parse(element.dataset.fsc)
	  var data = []
	  for(var point of fsc['fsc']){
		  data.push({x:Number(1 / Number(point['resolution'])), y:point['correlation']})
	  }
	  this.toggleMenu(element, 'fsccontrolpane') 
	  var chart = new Chartist.Line('#fscchart', {
          series: [ data ]
        }, {
          axisX: {
            type: Chartist.AutoScaleAxis,
            onlyInteger: true
          },
          axisY: {
            type: Chartist.AutoScaleAxis,
            onlyInteger: false
          },
          fullWidth: true,
          showLine: true,
          chartPadding: {
            right: 40
          },
          showPoint: false,
          plugins: [
		Chartist.plugins.verticalLine({
			position: '0.143'
		}),
		 Chartist.plugins.ctAxisTitle({
      			axisX: {
        			axisTitle: 'RECIPROCAL RESOLUTION',
        			axisClass: 'ct-axis-title',
        			offset: {
          				x: 0,
          				y: 20
        			},
        			textAnchor: 'middle'
     	 		},
      			axisY: {
        			axisTitle: 'FSC',
        			axisClass: 'ct-axis-title',
        			offset: {
          				x: 0,
          				y: 0
        			},
        			textAnchor: 'middle',
        			flipTitle: false
      			}
    		})
	]
	})
  }
  
  getCluster2D(projfile){
	  fetcher.killWorker()
	  this.projfile = projfile
	  this.updateTitle('2D Clusters')
	  this.resetSideBar()
	  this.showLoader()
	  var request = {
		  cls : "view",
		  fnc : "getViewSimpleCluster2D",
		  arg : {
			  projfile : projfile
			  }
	  }
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			document.getElementById('mainpane').innerHTML = json.html
			var iterationselector = document.getElementById('iterationselector')
			iterationselector.selectedIndex = iterationselector.options.length -1
			simpleview.getCluster2DIteration(iterationselector, projfile)
		})  
  }
  
  getCluster2DIteration(caller, projfile){
	  fetcher.killWorker()
	  this.resetSideBar()
	  var thumbnails = document.getElementById('thumbnails')
	  thumbnails.innerHTML = ""
	  var element = document.getElementById('iterationselector')
	  var selected = element.options[element.selectedIndex]
	  var header = JSON.parse(selected.dataset.header)
	  var stats = selected.dataset.stats ? JSON.parse(selected.dataset.stats) : false
	  var xdim = (header.nx < 100) ? header.nx : 100;
	  var ydim = (header.ny < 100) ? header.ny : 100;
	  var wfilt = false
	  if(document.getElementById('wfilt') && document.getElementById('wfilt').checked){
		  var wfilt = true
	  }
	  for(var thumbcount = 0; thumbcount < header['nz']; thumbcount++){
		  var thumbnail = document.createElement('div')
		  thumbnail.className = "thumbnail"
		  var container = document.createElement('div')
		  container.className = "thumbnailcontainer"
		  container.style.width = xdim + "px"
		  container.style.height = ydim + "px"
		  var containerloader = document.createElement('div')
		  containerloader.className = "miniloader"
		  container.appendChild(containerloader)
		  var containerimg = document.createElement('img')
		  containerimg.className = "thumbnailcontainerimg cls2D"
		  if(wfilt){
			containerimg.dataset.path = window.location.href + "/image?stackfile=" + selected.dataset.file.replace(".mrc", "_wfilt.mrc") + "&frame=" + thumbcount
		  }else{
		  	containerimg.dataset.path = window.location.href + "/image?stackfile=" + selected.dataset.file + "&frame=" + thumbcount
		  }
	          container.appendChild(containerimg)
		  
		  if(stats){
			  for(var key of Object.keys(stats[thumbcount])){
				  thumbnail.title += key + " : " + stats[thumbcount][key] + "\n"
				  thumbnail.dataset[key] = Number(stats[thumbcount][key])
			  }
			 thumbnail.setAttribute("onclick", "simpleview.toggleThumbnail(this)")
			 var img = document.createElement('img')
			 img.src = "img/view-out.png"
			 img.className = "viewptcls"
			 img.dataset.classid = stats[thumbcount]['class']
			 img.dataset.projfile = projfile
			 img.setAttribute("onclick", "simpleview.getClassParticles(this)")
			 container.appendChild(img)
		  }
		  
		  var controls = this.createControlsDiv("clusterparticles")
		  controls.appendChild(this.createControlsHeader("Particles"))
		  var body = document.createElement('div')
		  var div = document.createElement('div')
		  div.id = "particles"
		  body.appendChild(div)
		  controls.appendChild(body)
		  document.getElementById('mainpane').appendChild(controls)
		  thumbnail.appendChild(container)
		  thumbnails.appendChild(thumbnail)
	  }
	  this.addImageControls(['cls2D'], true, false)
	  if(stats){
		  this.addSortControls(stats)
		  this.addSelectControls(stats)
		  this.addPlotControls(stats)
		  this.addSaveControls("cls2D")  
		  this.updateParticleCount()
		  document.getElementById('sortattribute').value = "res"
		  fetcher.imagewidth = xdim
		  this.sortThumbnails()
	  }
	  fetcher.imagewidth = xdim
	  fetcher.fetchImages()
  }
  
  getClassParticles(element){
	  this.getParticles(element.dataset.projfile, {ptclselector:'class', ptcloperator:'eq', ptclvalue:element.dataset.classid})
  }

  getParticles(projfile, selection){
	  fetcher.killWorker()
	  this.projfile = projfile
	  this.updateTitle('Particles')
	  this.resetSideBar()
	  this.showLoader()
	  var request = {
		  cls : "view",
		  fnc : "getViewSimpleParticles",
		  arg : {
			  projfile : projfile
			  }
	  }
		fetcher.fetchJSON(request)
		.then(response => response.json())
		.then (json => {
			document.getElementById('mainpane').innerHTML = json.html
		//	var selector = document.getElementById('stkselector')
		//	selector.selectedIndex = selector.options.length -1
		//	simpleview.getStkParticles(selector)
			this.particles = json.ptcls
			this.stks = json.stks
			document.getElementById('particlecount').innerHTML = "Showing 0 / " + this.particles.length
			if(selection){
				document.getElementById('ptclselector').value = selection['ptclselector']
				document.getElementById('ptcloperator').value = selection['ptcloperator']
				document.getElementById('ptclvalue').value = selection['ptclvalue']
				this.getStkParticles()
			}
		})  
  }
 
  getStkParticles(){
	  fetcher.killWorker()
	  this.resetSideBar()
	  var thumbnails = document.getElementById('thumbnails')
	  thumbnails.innerHTML = ""
	  var ptclselector = document.getElementById('ptclselector').value
	  var ptcloperator = document.getElementById('ptcloperator').value
	  var ptclvalue = document.getElementById('ptclvalue').value
	  
	  var xdim = (this.stks[0]['header']['nx'] < 100) ? this.stks[0]['header']['nx'] : 100;
	  var ydim = (this.stks[0]['header']['ny'] < 100) ? this.stks[0]['header']['ny'] : 100;
	  var particlecount = 0
	  
	  for(var particle in this.particles){
		  var view = false
		  if(ptcloperator == "eq" && Number(this.particles[particle][ptclselector]) == ptclvalue){
			  view = true
		  }else if(ptcloperator == "lt" && Number(this.particles[particle][ptclselector]) < ptclvalue){
			  view = true
		  }else if(ptcloperator == "gt" && Number(this.particles[particle][ptclselector]) > ptclvalue){
			  view = true
		  }
		  
		  if(view){
			  particlecount++
			  var thumbnail = document.createElement('div')
			  thumbnail.className = "thumbnail"
			  var container = document.createElement('div')
			  container.className = "thumbnailcontainer"
			  container.style.width = xdim + "px"
			  container.style.height = ydim + "px"
			  var containerloader = document.createElement('div')
			  containerloader.className = "miniloader"
			  container.appendChild(containerloader)
			  var containerimg = document.createElement('img')
			  containerimg.className = "thumbnailcontainerimg cls2D"
			  containerimg.dataset.path = window.location.href + "/image?stackfile=" + this.stks[Number(this.particles[particle]['stkind']) - 1]['file']  + "&frame=" + Number(particle - Number(this.stks[Number(this.particles[particle]['stkind']) - 1]['fromp']) + 1)
			  for(var key of Object.keys(this.particles[particle])){
				 thumbnail.title += key + " : " + this.particles[particle][key] + "\n"
				 thumbnail.dataset[key] = Number(this.particles[particle][key])
			  }
			  container.appendChild(containerimg)
			  thumbnail.appendChild(container)
			  thumbnails.appendChild(thumbnail)
			  if(particlecount > 2000){
				  alert("Too many particles (>2000) to view efficiently. Truncating view to first 2000 particles")
				  break
			  }
		  }
	  }
	//  if(particlecount > 2000){
	//	  alert("Too many particles (>2000) to view efficiently. Please refine your view criteria")
//	  }else{
		  document.getElementById('particlecount').innerHTML = "Showing " + particlecount + " / " + this.particles.length
		  this.addImageControls(['cls2D'], true, false)
		  this.addSortControls(this.particles)
		  this.addPlotControls(this.particles)
		  fetcher.imagewidth = xdim
		  fetcher.fetchImages()
	//  }
  }
  
  toggleThumbnail(element){
	if(Number(element.dataset.state) > 0){
		element.dataset.state = 0
	//	element.style.backgroundImage = "url('img/cancel.png')"
	}else{
		element.style.backgroundImage = ""
		element.dataset.state = 1
	}
	if(document.getElementById('particlecount')){
		this.updateParticleCount()
	}
	if(document.getElementById('micrographcount')){
		this.updateMicrographCount()
	}
  }
  
  updateParticleCount(){
	var particlecount = document.getElementById('particlecount')
	var total = 0
	var count = 0
	for(var thumbnail of document.getElementsByClassName('thumbnail')){
		if(thumbnail.dataset.pop){
			total += Number(thumbnail.dataset.pop)
			if(Number(thumbnail.dataset.state) > 0){
				count += Number(thumbnail.dataset.pop)
			}
		} 
		else if(thumbnail.dataset.nptcls){
			total += Number(thumbnail.dataset.nptcls)
			if(Number(thumbnail.dataset.state) > 0){
				count += Number(thumbnail.dataset.nptcls)
			}
		}
	}
	if(total != 0){
		particlecount.innerHTML = "Particle Count : " + count + " / " + total
	}
  }
    
  updateMicrographCount(){
	var micrographcount = document.getElementById('micrographcount')
	var total = 0
	var count = 0
	for(var thumbnail of document.getElementsByClassName('thumbnail')){
		if(thumbnail.dataset.state != 0){
			count++
		}
		total++
	}
	micrographcount.innerHTML = "Micrograph Count : " + count + " / " + total
  }

  hideMenus(){
	for(var menu of document.getElementsByClassName('controlpane')){
		menu.style.display = "none"
	}
	for(var menu of document.getElementById('sidebar').getElementsByTagName('div')){
		menu.dataset.selected = false
	}
  }
  
  toggleMenu(element, menu){
	  var thismenu = document.getElementById(menu)
	  if(thismenu.style.display == "block"){
		  this.hideMenus()
	  }else{
		  this.hideMenus()
		  thismenu.style.display = "block"
		  element.dataset.selected = true
	  }
  }
  
  createControlsIcon(image, title, target){
	  var div = document.createElement('div')
	  div.title = title
	  div.setAttribute("onclick" ,"simpleview.toggleMenu(this, '" + target + "')")
	  var img = document.createElement('img')
	  img.src = image
	  img.className = "menuicon"
	  div.appendChild(img)
	  return div
  }
  
  createControlsDiv(id){
	  var div = document.createElement('div')
	  div.className = "controlpane"
	  div.id = id
	  return div
  }
  
  createControlsHeader(text){
	  var div = document.createElement('div')
	  div.className = "viewheader"
	  var logo = document.createElement('img')
	  logo.src = "img/square_logo_small.png"
	  logo.id = "simplepagelogo"
	  var headertext = document.createElement('h')
	  headertext.innerHTML = text
	  div.appendChild(logo)
	  div.appendChild(headertext)
	  return div
  }
  
  createControlsButtons(){
	 var div = document.createElement('div')
	 div.className = "buttons"
	 return div
  }
  
  createControlsButton(text, func){
	  var div = document.createElement('div')
	  div.innerHTML = text
	  div.setAttribute("onclick", func)
	  return div
  }
  
  addImageControls(targets, d2controls, d3controls){
	  document.getElementById('sidebar').appendChild(this.createControlsIcon("img/view-out.png", "View", "imagecontrolpane"))
	  var controls = this.createControlsDiv("imagecontrolpane")
	  controls.appendChild(this.createControlsHeader("View"))
	  if(d2controls){
		  var body = document.createElement('div')
		  body.className = "controlbody"
		  var div = document.createElement('div')
		  div.className = "controlbodylabel"
		  div.innerHTML = "2D Controls"
		  body.appendChild(div)
		  var controlbodyelements = document.createElement('div')
		  controlbodyelements.className = "controlbodyelements"
		  
		  var div = document.createElement('div')
		  var label = document.createElement('label')
		  label.innerHTML = "Target"
		  var select = document.createElement('select')
		  select.id = "target"
		  select.setAttribute("onclick", "simpleview.updateTarget(this)")
		  for(var target of targets){
			  var option = document.createElement('option')
			  option.value = target
			  option.innerHTML = target
			  select.appendChild(option)
		  }
		  div.appendChild(label)
		  div.appendChild(select)
		  controlbodyelements.appendChild(div)
		  
		  var div = document.createElement('div')
		  var label = document.createElement('label')
		  label.innerHTML = "Contrast"
		  var input = document.createElement('input')
		  input.id = "contrast"
		  input.type="range"
		  input.min="0"
		  input.max="500"
		  input.step="1"
		  input.value="100"
		  input.setAttribute("oninput", "simpleview.updateImages()")
		  div.appendChild(label)
		  div.appendChild(input)
		  controlbodyelements.appendChild(div)
		  
		  var div = document.createElement('div')
		  var label = document.createElement('label')
		  label.innerHTML = "Brightness"
		  var input = document.createElement('input')
		  input.id = "brightness"
		  input.type="range"
		  input.min="0"
		  input.max="200"
		  input.step="1"
		  input.value="100"
		  input.setAttribute("oninput", "simpleview.updateImages()")
		  div.appendChild(label)
		  div.appendChild(input)
		  controlbodyelements.appendChild(div)
		  
		  var div = document.createElement('div')
		  var label = document.createElement('label')
		  label.innerHTML = "Zoom"
		  var input = document.createElement('input')
		  input.id = "zoom"
		  input.type="range"
		  input.min="0.2"
		  input.max="2"
		  input.step="0.1"
		  input.value="1"
		  input.setAttribute("oninput", "simpleview.updateImages()")
		  div.appendChild(label)
		  div.appendChild(input)
		  controlbodyelements.appendChild(div)
		  
		  var div = document.createElement('div')
		  var label = document.createElement('label')
		  label.innerHTML = "Scale"
		  var input = document.createElement('input')
		  input.id = "scale"
		  input.type="range"
		  input.min="50"
		  input.max="500"
		  input.step="10"
		  input.value="200"
		  input.setAttribute("oninput", "simpleview.scaleImages()")
		  div.appendChild(label)
		  div.appendChild(input)
		  controlbodyelements.appendChild(div)
		  
		  var div = document.createElement('div')
		  var label = document.createElement('label')
		  label.innerHTML = "View Deselected"
		  var input = document.createElement('input')
		  input.id = "showselected"
		  input.type="checkbox"
		  input.checked = true
		  input.setAttribute("oninput", "simpleview.updateImages()")
		  div.appendChild(label)
		  div.appendChild(input)
		  controlbodyelements.appendChild(div)
		  
		  for(var target of targets){
			  var div = document.createElement('div')
			  var label = document.createElement('label')
			  label.innerHTML = "Show " + target
			  var input = document.createElement('input')
			  input.className = "showtarget"
			  input.dataset.target = target
			  input.type="checkbox"
			  input.checked = true
			  input.setAttribute("oninput", "simpleview.updateImages()")
			  div.appendChild(label)
			  div.appendChild(input)
			  controlbodyelements.appendChild(div)
		  }
		  
		  body.appendChild(controlbodyelements)
		  controls.appendChild(body)
		  
		  
	  }
	  if(d3controls){
		  var body2 = document.createElement('div')
		  body2.className = "controlbody"		  
		  var div = document.createElement('div')
		  div.className = "controlbodylabel"
		  div.innerHTML = "3D Controls"
		  body2.appendChild(div)
		  var controlbodyelements = document.createElement('div')
		  controlbodyelements.className = "controlbodyelements"
		  
		  var div = document.createElement('div')
		  var label = document.createElement('label')
		  label.innerHTML = "Sigma"
		  var input = document.createElement('input')
		  input.id = "sigma"
		  input.type="range"
		  input.min="0.1"
		  input.max="20"
		  input.step="0.1"
		  input.value="1"
		  input.setAttribute("oninput", "simpleview.updateIsoSurface(this)")
		  div.appendChild(label)
		  div.appendChild(input)
		  controlbodyelements.appendChild(div)
		  
		  var div = document.createElement('div')
		  var label = document.createElement('label')
		  label.innerHTML = "Absolute"
		  var input = document.createElement('input')
		  input.id = "absolute"
		  input.type = "text"
		  input.disabled = true
		  div.appendChild(label)
		  div.appendChild(input)
		  controlbodyelements.appendChild(div)
		  
		  var div = document.createElement('div')
		  var label = document.createElement('label')
		  label.innerHTML = "Quality"
		  var input = document.createElement('input')
		  input.id = "quality"
		  input.type="range"
		  input.min="1"
		  input.max="6"
		  input.step="1"
		  input.value="4"
		  input.setAttribute("oninput", "simpleview.reloadIsoSurface()")
		  div.appendChild(label)
		  div.appendChild(input)
		  controlbodyelements.appendChild(div)
		 
		  body2.appendChild(controlbodyelements)
		  controls.appendChild(body2)
	  }
	
	
	  document.getElementById('mainpane').appendChild(controls)
	  
  }
  
  addSortControls(stats){
	  document.getElementById('sidebar').appendChild(this.createControlsIcon("img/sort.png", "Sort", "sortcontrolpane"))
	  var controls = this.createControlsDiv("sortcontrolpane")
	  controls.appendChild(this.createControlsHeader("Sort"))
	  
	  var controlbody = document.createElement('div')
	  controlbody.className = "controlbody"
	  
	  var controlbodyelements = document.createElement('div')
	  controlbodyelements.className = "controlbodyelements"
	  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Sort on"
	  var select = document.createElement('select')
	  select.id = "sortattribute"
	  for(var key of Object.keys(stats[0])){
		  var option = document.createElement('option')
		  option.value = key
		   option.innerHTML = key
		  select.appendChild(option)
	  }
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)
	  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Order"
	  var select = document.createElement('select')
	  select.id = "sortorder"
	  for(var key of ["ascending", "descending"]){
		  var option = document.createElement('option')
		  option.value = key
		  option.innerHTML = key
		  select.appendChild(option)
	  }
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)
	  
	  controlbody.appendChild(controlbodyelements)
	  
	  var buttons = this.createControlsButtons()
	  buttons.appendChild(this.createControlsButton("Sort", "simpleview.sortThumbnails()"))
	  controls.appendChild(controlbody)
	  controls.appendChild(buttons)
	  document.getElementById('mainpane').appendChild(controls)
  }
  
  addSelectControls(stats){
	  document.getElementById('sidebar').appendChild(this.createControlsIcon("img/select.png", "Select", "selectcontrolpane"))
	  var controls = this.createControlsDiv("selectcontrolpane")
	  controls.appendChild(this.createControlsHeader("Select"))
	  
	  var controlbody = document.createElement('div')
	  controlbody.className = "controlbody"
	  
	  var controlbodyelements = document.createElement('div')
	  controlbodyelements.className = "controlbodyelements"
	  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Attribute"
	  var select = document.createElement('select')
	  select.id = "selectattribute"
	  for(var key of Object.keys(stats[0])){
		  var option = document.createElement('option')
		  option.value = key
		  option.innerHTML = key
		  select.appendChild(option)
	  }
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)

	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Operator"
	  var select = document.createElement('select')
	  select.id = "selectoperator"
	  for(var key of ["gt", "lt", "eq"]){
		  var option = document.createElement('option')
		  option.value = key
		  option.innerHTML = key
		  select.appendChild(option)
	  }
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)
	  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Value"
	  var select = document.createElement('input')
	  select.id = "selectvalue"
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)
	    
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Append"
	  var select = document.createElement('input')
	  select.type = "checkbox"
	  select.id = "selectappend"
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)
	  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Invert Selection"
	  var select = document.createElement('input')
	  select.type = "button"
	  select.value = 'Invert'
	  select.onclick = (() => {
		  var thumbnails = document.getElementsByClassName('thumbnail')
		  for(var thumbnail of thumbnails){
			if(thumbnail.dataset.state == 0){
				thumbnail.dataset.state = 1
			}else {
				thumbnail.dataset.state = 0
			}
		  }
	  })
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)
	  
	  controlbody.appendChild(controlbodyelements)
	  
	  var buttons = this.createControlsButtons()
	  buttons.appendChild(this.createControlsButton("Select", "simpleview.selectThumbnails()")) 
	  controls.appendChild(controlbody)
	  controls.appendChild(buttons)
	  document.getElementById('mainpane').appendChild(controls)
  }
  
  addPlotControls(stats){
	  document.getElementById('sidebar').appendChild(this.createControlsIcon("img/plot.png", "Plot", "plotcontrolpane"))
	  var controls = this.createControlsDiv("plotcontrolpane")
	  controls.appendChild(this.createControlsHeader("Plot"))
	  
	  var controlbody = document.createElement('div')
	  controlbody.className = "controlbody"
	  
	  var controlbodyelements = document.createElement('div')
	  controlbodyelements.className = "controlbodyelements"
	  
	  var div = document.createElement('div')
	  div.id = 'xaxis'
	  div.style.display = 'none'
	  var label = document.createElement('label')
	  label.innerHTML = "X Axis"
	  var select = document.createElement('select')
	  select.id = "plotx"
	  for(var key of Object.keys(stats[0])){
		  var value = Number(stats[0][key])
		  if(!Number.isNaN(value)){
			  var option = document.createElement('option')
			  option.value = key
			  option.innerHTML = key
			  select.appendChild(option)
		  }
	  }
	  select.setAttribute("oninput", "simpleview.plot()")
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)
	  
	  var div = document.createElement('div')
	  div.id = 'yaxis'
	  var label = document.createElement('label')
	  label.innerHTML = "Y Axis"
	  var select = document.createElement('select')
	  select.id = "ploty"
	  for(var key of Object.keys(stats[0])){
		var value = Number(stats[0][key])
		if(!Number.isNaN(value)){ 
		  var option = document.createElement('option')
		  option.value = key
		  option.innerHTML = key
		  select.appendChild(option)
	    }
	  }
	  select.setAttribute("oninput", "simpleview.plot()")
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)
	  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Type"
	  var select = document.createElement('select')
	  select.id = "plottype"
	  for(var key of ["bar", "scatter"]){
		  var option = document.createElement('option')
		  option.value = key
		  option.innerHTML = key
		  select.appendChild(option)
	  }
	  select.setAttribute("oninput", "simpleview.togglePlot(this)")
	//  select.setAttribute("oninput", "simpleview.plot()")
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)
	   
	  var div = document.createElement('div')
	  div.id = 'resetzoom'
	  div.style.display = 'none'
	  var label = document.createElement('label')
	  var input = document.createElement('input')
	  label.innerHTML = "Reset Zoom"
	  div.appendChild(label)
	  input.type = 'button'
	  input.value = 'Reset'
	  input.setAttribute('onclick',"simpleview.plot()")
	  div.appendChild(input)
	  controlbodyelements.appendChild(div)
	  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Log Scale"
	  var input = document.createElement('input')
	  input.type = 'checkbox'
	  input.id = 'logscale'
	  input.setAttribute('onclick',"simpleview.plot()")
	  div.appendChild(label)
	  div.appendChild(input)
	  controlbodyelements.appendChild(div)
	  
	  var div = document.createElement('div')
	  div.id = "chart"

	  controlbody.appendChild(controlbodyelements)
	  
	/*  var buttons = this.createControlsButtons()
	  buttons.appendChild(this.createControlsButton("Plot", "simpleview.plot()"))
	  controls.appendChild(buttons)
	  */
	  controls.appendChild(controlbody)
	  controlbody.appendChild(div)
	  document.getElementById('mainpane').appendChild(controls)
  }
  
  togglePlot(element){
	  if(element.value == 'scatter'){
		  document.getElementById('xaxis').style.display = 'flex'
		  document.getElementById('resetzoom').style.display = 'flex'
	  }else{
		  document.getElementById('xaxis').style.display = 'none'
		  document.getElementById('resetzoom').style.display = 'none'
	  }
	  simpleview.plot()
  }
  
  addSaveControls(oritype){
	  document.getElementById('sidebar').appendChild(this.createControlsIcon("img/save.png", "Save", "savecontrolpane"))
	  var controls = this.createControlsDiv("savecontrolpane")
	  controls.appendChild(this.createControlsHeader("Save")) 
	  
	  var controlbody = document.createElement('div')
	  controlbody.className = "controlbody"
	  
	  var controlbodyelements = document.createElement('div')
	  controlbodyelements.className = "controlbodyelements"
	  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Description"
	  var select = document.createElement('textarea')
	  select.id = "description"
	  div.appendChild(label)
	  div.appendChild(select)
	  controlbodyelements.appendChild(div)
	  
	  if(oritype != 'manualpick'){
		  var div = document.createElement('div')
		  var label = document.createElement('label')
		  label.innerHTML = "Export Class Averages"
		  label.title = 'Saves selected averages as a .mrc stack file'
		  var input = document.createElement('input')
		  input.type = 'checkbox'
		  input.id = "saveclusters"
		  input.title = 'Saves selected averages as a .mrc stack file'
		  div.appendChild(label)
		  div.appendChild(input)
		  controlbodyelements.appendChild(div)
		  
		//  var div = document.createElement('div')
		//  var label = document.createElement('label')
		  //label.innerHTML = "Save STAR"
		//  var input = document.createElement('input')
		//  input.type = 'checkbox'
		//  input.id = "savestar"
		//  div.appendChild(label)
		//  div.appendChild(input)
		//  controlbodyelements.appendChild(div)
	  }
	  
	  controlbody.appendChild(controlbodyelements)
	  
	  var buttons = this.createControlsButtons()
	  if(oritype == 'manualpick'){
		 var button = this.createControlsButton("Save", "simpleview.saveBoxes()")
		 button.id = 'savebutton'
		 buttons.appendChild(button) 
	  }else{
		var button = this.createControlsButton("Save", "simpleview.saveSelection('" + oritype + "')")
		button.id = 'savebutton'
		buttons.appendChild(button)
	  }
	  controls.appendChild(controlbody)
	  controls.appendChild(buttons)
	  document.getElementById('mainpane').appendChild(controls)
  }
  
  addBoxControls(){
	  document.getElementById('sidebar').appendChild(this.createControlsIcon("img/select.png", "Box Controls", "boxcontrolpane"))
	  var controls = this.createControlsDiv("boxcontrolpane")
	  controls.appendChild(this.createControlsHeader("Box"))
	  
	  var body = document.createElement('div')
	  body.className = "controlbody"
		 
	  var controlbodyelements = document.createElement('div')
	  controlbodyelements.className = "controlbodyelements"

	  var div = document.createElement('div')
      var label = document.createElement('label')
      label.innerHTML = "Helical Mode"
      var input = document.createElement('input')
	  input.type = 'checkbox'
	  input.name = 'boxmode'
	  input.value = 'helical'
	  input.onclick = function(element){ 
		  if(element.target.checked){
				simpleview.manualPickHelix()
				document.getElementById("helixboxstep").style.display = 'flex'
				document.getElementById("showboxes").checked = false
				document.getElementById("showradii").checked = false
		  }else{
			  document.getElementById("pickingcanvas").onmouseup = function(){}
			  document.getElementById("pickingcanvas").onmousedown = function(){}
			  document.getElementById("pickingcanvas").onmousemove = function(){}
			  document.getElementById("helixboxstep").style.display = 'none'
			  document.getElementById("showboxes").checked = true
			  document.getElementById("showradii").checked = true
		  }
		  simpleview.drawBoxes()
	  }
	  div.appendChild(label)
          div.appendChild(input)
          controlbodyelements.appendChild(div)
		    
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Box Size (pixels)"
	  var input = document.createElement('input')
	  input.id = "boxsize"
	  input.value="100"
	  input.onchange = function () {
		simpleview.findBoxes()
		simpleview.drawBoxes()
	  }
	  div.appendChild(label)
	  div.appendChild(input)
	  controlbodyelements.appendChild(div)
	  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Box Step (pixels)"
	  var input = document.createElement('input')
	  input.id = "boxstep"
	  input.value="10"
	  input.onchange = function () {
		simpleview.findAllBoxes()
		simpleview.drawBoxes()
	  }
	  div.appendChild(label)
	  div.appendChild(input)
	  div.style.display = 'none'
	  div.id = 'helixboxstep'
	  controlbodyelements.appendChild(div)
	  
		  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Particle Radius (pixels)"
	  var input = document.createElement('input')
	  input.id = "ptclradius"
	  input.value="50" 
	  input.setAttribute("oninput", "simpleview.drawBoxes()")
	  div.appendChild(label)
	  div.appendChild(input)
	  controlbodyelements.appendChild(div)
		  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Line Width"
	  var input = document.createElement('input')
	  input.id = "linewidth"
	  input.type = 'range'
	  input.min="1"
	  input.max="10"
	  input.step="1"
	  input.value="5"
	  input.setAttribute("oninput", "simpleview.drawBoxes()")
	  div.appendChild(label)
	  div.appendChild(input)
	  controlbodyelements.appendChild(div)
	   
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Show Boxes"
	  var input = document.createElement('input')
	  input.id = "showboxes"
	  input.type="checkbox"
	  input.checked = true
	  input.setAttribute("oninput", "simpleview.drawBoxes()")
	  div.appendChild(label)
	  div.appendChild(input)
	  controlbodyelements.appendChild(div)
	  
	  var div = document.createElement('div')
	  var label = document.createElement('label')
	  label.innerHTML = "Show Particle Radius"
	  var input = document.createElement('input')
	  input.id = "showradii"
	  input.type="checkbox"
	  input.checked = true
	  input.setAttribute("oninput", "simpleview.drawBoxes()")
	  div.appendChild(label)
	  div.appendChild(input)
	  controlbodyelements.appendChild(div)
	  
	  var div = document.createElement('div')
	  div.className = 'buttons'
	  div.appendChild(this.createControlsButton("Clear All Boxes", "simpleview.clearBoxes()"))
	  controlbodyelements.appendChild(div)
	  
	  body.appendChild(controlbodyelements)
	  controls.appendChild(body)
	  document.getElementById('mainpane').appendChild(controls)
  
  }

  sortThumbnails() {
    if (fetcher.loaderworker != undefined){
      fetcher.loaderworker.terminate()
    }
    var attribute = document.getElementById('sortattribute').value
    var order = document.getElementById('sortorder').value
    var thumbnails = document.getElementsByClassName('thumbnail')
    var sorted = []
    for (var thumbnail of thumbnails){
      sorted.push([thumbnail.getAttribute('data-' + attribute), thumbnail])
    }
    for (var i = thumbnails.length - 1 ; i >= 0; i--){
      thumbnails[i].remove()
    }
    if(order == "ascending"){
      sorted.sort((x,y) => { return x[0] - y[0] })
    }else if(order == "descending"){
      sorted.sort((x,y) => { return y[0] - x[0] })
    }
    for(var element of sorted){
      document.getElementById('thumbnails').appendChild(element[1])
    }
    fetcher.fetchImages()
  }
  
  selectThumbnails() {
    var attribute = document.getElementById('selectattribute').value
    var operator = document.getElementById('selectoperator').value
    var value = Number(document.getElementById('selectvalue').value)
    var thumbnails = document.getElementsByClassName('thumbnail')

    for (var thumbnail of thumbnails){
      var thumbvalue = Number(thumbnail.getAttribute('data-' + attribute))
      if(document.getElementById('selectappend').checked === false){
        thumbnail.dataset.state = 0
      }
      if(operator == "gt" && thumbvalue > value){
        thumbnail.dataset.state = 1
      }else if (operator == "lt" && thumbvalue < value){
        thumbnail.dataset.state = 1
      }else if (operator == "eq" && thumbvalue == value){
       thumbnail.dataset.state = 1
      }
    }
  }
  
  saveSelection(oritype){
	 document.getElementById('savebutton').innerHTML = 'Saving'
	// var savestar = document.getElementById('savestar').checked ? true : false
	var savestar = false
	 var saveclusters = document.getElementById('saveclusters').checked ? true : false
	 var selection = []
	 var thumbnails = document.getElementsByClassName('thumbnail')
	 for(var thumbnail of thumbnails){
		 if(thumbnail.dataset.class){
			selection.push({class:thumbnail.dataset.class, state:thumbnail.dataset.state})
		}else if(thumbnail.dataset.number){
			selection.push({class:thumbnail.dataset.number, state:thumbnail.dataset.state})
		}
	 }
	 selection.sort(function(a, b){return a['class'] - b['class']});
	 var parent = simpleview.projfile.split('/')
	 parent.pop()
	 parent = parent.pop()
	 var request = {
		  cls : "view",
		  fnc : "saveSimpleSelection",
		  arg : {
				type : "selection",
				executable : "simple_exec",
				projectfolder : project.selectedfolder,
				projecttable : project.selectedtable,
				projfile : simpleview.projfile,
				parent : parent,
				name : "Selection from " + parent,
				description : document.getElementById('description').value,
				keys : {
					  keyinfile : simpleview.projfile.replace(".simple", "_selected.txt"),
					  keyoritype : oritype
				},
				selection : selection,
				savestar : savestar,
				saveclusters : saveclusters
		  }
	 }
	fetcher.fetchJSON(request)
	.then(response => response.json())
	.then (json => {
		 alert("Selection job started")
		 document.getElementById('savebutton').innerHTML = 'Save'
		 project.refreshHistory()
	})
	.catch(() => {
		 alert("Error Saving Selection")
		 document.getElementById('savebutton').innerHTML = 'Save'
	})
  }

  plot(){
	var plotx = document.getElementById('plotx').value
	var ploty = document.getElementById('ploty').value
	var plottype = document.getElementById('plottype').value
	var logscale = document.getElementById('logscale')
	//var plotmax = document.getElementById('plotmax').value
	var thumbnails = document.getElementsByClassName('thumbnail')
	
	var labels = []
    var data = []
  
	if(plottype == "bar"){
      for(var thumbnail of thumbnails){
      //  labels.push(thumbnail.getAttribute('data-' + plotx))
        labels.push('')
        if(logscale.checked){
			data.push(Math.log(thumbnail.getAttribute('data-' + ploty)))
		}else{
			data.push(thumbnail.getAttribute('data-' + ploty))
		}
      }
      var chart = new Chartist.Bar('#chart', {
          labels: labels,
          series: [ data ]
        }, {
		  showLabel: false,
          fullWidth: true,
          chartPadding: {
            right: 40
          },
          plugins: [
			Chartist.plugins.zoom({ onZoom: onZoom })
		  ]
      })
    
      chart.on("draw", (data) => {
        if(data.type === "bar"){
          var index = data.index
          var classnumber = thumbnails[index].getAttribute('data-class')
          data.element._node.setAttribute('data-plotclass', classnumber)
          data.element._node.onclick = () => {
            var thumbnail = document.querySelector('[data-class="' + classnumber + '"]')
            thumbnail.style.backgroundColor = "red"
            thumbnail.scrollIntoView()
            setTimeout(() => {
              thumbnail.style.backgroundColor = "#f3f3f3"
            }, 1000)
          }
        }
      })
    } else if (plottype == "scatter"){
	  var data = []
	  for(var thumbnail of thumbnails){
		if(logscale.checked){
			data.push({x:thumbnail.getAttribute('data-' + plotx), y:Math.log(thumbnail.getAttribute('data-' + ploty))})
		}else{
			data.push({x:thumbnail.getAttribute('data-' + plotx), y:thumbnail.getAttribute('data-' + ploty)})
		}
        
      }
      var chart = new Chartist.Line('#chart', {
          series: [ data ]
        }, {
          fullWidth: true,
          showLine: false,
          chartPadding: {
            right: 40
          },
          axisX: {
            type: Chartist.AutoScaleAxis,
            onlyInteger: true
          },
          axisY: {
            type: Chartist.AutoScaleAxis,
            onlyInteger: false
          },
          plugins: [
			Chartist.plugins.zoom({ onZoom: onZoom })
		  ]
      });
      chart.on("draw", (data) => {
        if(data.type === "point"){
			var index = data.index
			data.element._node.setAttribute('data-index', index)
			data.element._node.setAttribute('data-plotclass', index)
			data.element._node.onclick = () => {
				var thumbnail = thumbnails[index]
				thumbnail.style.backgroundColor = "red"
				thumbnail.scrollIntoView()
				setTimeout(() => {
				  thumbnail.style.backgroundColor = "#f3f3f3"
				}, 1000)
			}
		}
      })  
	}
    function onZoom(chart, reset) {
		resetFnc = reset;
	}
  }
  
  updateTarget(element){
	var selected = element.options[element.selectedIndex]

	document.getElementById('brightness').value = Number(selected.dataset.brightness)
    document.getElementById('contrast').value = Number(selected.dataset.contrast)
    document.getElementById('zoom').value = Number(selected.dataset.zoom)
  }

  updateImages() {
	var target = document.getElementById('target')
    var brightness = document.getElementById('brightness').value
    var contrast = document.getElementById('contrast').value
    var zoom = document.getElementById('zoom').value
    var showselected = document.getElementById('showselected').checked
    
	target.options[target.selectedIndex].dataset.brightness = brightness
	target.options[target.selectedIndex].dataset.contrast = contrast
	target.options[target.selectedIndex].dataset.zoom = zoom
	
    var images = document.getElementsByClassName(target.value)
    
    for(var image of images){
      image.style.filter = "contrast(" + contrast + "%) brightness(" + brightness + "%)"
      image.style.transform = "scale(" + zoom + ")"
    }
    
    for (var rule of document.styleSheets[0]['rules']){
		if(rule['selectorText'] == '.simpleview .thumbnail[data-state="0"]'){
			if(showselected){
				rule.style.display = "flex"
			}else{
				rule.style.display = "none"
			}
		}
	}
	
	for (var showtarget of document.getElementsByClassName('showtarget')){
		var targetname = showtarget.dataset.target
		var images = document.getElementsByClassName(targetname)
		for(var image of images){
			if(showtarget.checked){
				image.parentElement.style.display = "flex"
			}else{
				image.parentElement.style.display = "none"
			}
		}
	}
	
  }
  
  scaleImages() {
	  var size = document.getElementById('scale').value
	  
	  for(var thumbnailcontainerimg of document.getElementsByClassName('thumbnailcontainerimg')){
		  thumbnailcontainerimg.src = null
		  thumbnailcontainerimg.parentElement.getElementsByClassName('miniloader')[0].style.display = "block"
	  }
		  
	  for(var thumbnailcontainer of document.getElementsByClassName('thumbnailcontainer')){
		  thumbnailcontainer.style.width = size + "px"
		  thumbnailcontainer.style.height = size + "px"
	  }
	  
	  fetcher.imagewidth = Number(size)
	  fetcher.fetchImages()
  }
  
  view2D(){
	  var request = {
		  cls : "view",
		  fnc : "view2D",
		  arg : {
				filename : browser.selection
		  }
	 }
	fetcher.fetchJSON(request)
	.then(response => response.json())
	.then(json => {
		document.getElementById('popup').style.display = "block"
		document.getElementById('popupwindow').className = "simpleview"
		document.getElementById('popupwindow').innerHTML = json.html
		
		var thumbnails = document.getElementById('thumbnails')
		
		var xdim = (json.header.nx < 100) ? json.header.nx : 100;
		var ydim = (json.header.ny < 100) ? json.header.ny : 100;
		for(var thumbcount = 0; thumbcount < json.header['nz']; thumbcount++){
		  var thumbnail = document.createElement('div')
		  thumbnail.className = "thumbnail"
		  var container = document.createElement('div')
		  container.className = "thumbnailcontainer"
		  container.style.width = xdim + "px"
		  container.style.height = ydim + "px"
		  var containerloader = document.createElement('div')
		  containerloader.className = "miniloader"
		  container.appendChild(containerloader)
		  var containerimg = document.createElement('img')
		  containerimg.className = "thumbnailcontainerimg cls2D"
		  containerimg.dataset.path = window.location.href + "/image?stackfile=" + browser.selection + "&frame=" + thumbcount
		  container.appendChild(containerimg)
		  thumbnail.appendChild(container)
		  thumbnails.appendChild(thumbnail)
	  }
	  this.addImageControls(['cls2D'], true, false)
	  fetcher.imagewidth = xdim
	  fetcher.fetchImages()
	  browser.hide()
    })
  }
  
  view3D(){
	  var request = {
		  cls : "view",
		  fnc : "view2D",
		  arg : {
				filename : browser.selection
		  }
	 }
	fetcher.fetchJSON(request)
	.then(response => response.json())
	.then(json => {
		document.getElementById('popup').style.display = "block"
		document.getElementById('popupwindow').className = "simpleview"
		document.getElementById('popupwindow').innerHTML = json.html
		
		browser.hide()
		
		var request = {
		  cls : "view",
		  fnc : "prepareMRCVolume",
		  arg : {
			  volume : browser.selection
		   }
		}
		return fetcher.fetchJSON(request)
	})
	.then(response => response.json())
	.then (json => {
		  this.resetSideBar()
		  if(this.plugin){
			  this.plugin.destroy()
		  }
		  this.plugin = LiteMol.Plugin.create({ target: '#mainpane' })
		  this.mdb = json['mdb']
		  var url = window.location.href + "/DensityServer/local/"
		  var surface = this.plugin.createTransform()
		//  .add(this.plugin.root, LiteMol.Bootstrap.Entity.Transformer.Data.Download, { url: "/DensityServer/local/" + json['mdb'] + "/cell?detail=4", type: 'Binary', description: 'local Density', title: "joe2", id:"density3ddata"})
		  .add(this.plugin.root, LiteMol.Bootstrap.Entity.Transformer.Data.Download, { url: url.replace('//DensityServer', '/DensityServer') + json['mdb'] + "/cell?detail=4", type: 'Binary', description: 'local Density', title: "joe2", id:"density3ddata"})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Data.ParseBinaryCif, {id:"density3dcif"}, {})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Density.CreateFromCif, {id:"density3dd", blockIndex: 1})
		  .then(LiteMol.Bootstrap.Entity.Transformer.Density.CreateVisual, {
			style: LiteMol.Bootstrap.Visualization.Density.Style.create({
			  isoValue: 0,
			  isoValueType: LiteMol.Bootstrap.Visualization.Density.IsoValueType.Sigma,
			  color: LiteMol.Visualization.Color.fromHex(0xBB3333),
			  isWireframe: false,
			  transparency: { alpha: 1.0 }
			  })
			}, {ref : "density3d"})
		  this.plugin.applyTransform(surface);
		  LiteMol.Bootstrap.Command.Layout.SetViewportOptions.dispatch(this.plugin.context, { clearColor: LiteMol.Visualization.Color.fromRgb(255, 255, 255) });
		  LiteMol.Bootstrap.Command.Visual.ResetScene.dispatch(this.plugin.context, void 0)
		  this.addImageControls([], false, true)
	 })
  }
}

var simpleview

window.addEventListener('load', () => {simpleview = new SimpleView()});
