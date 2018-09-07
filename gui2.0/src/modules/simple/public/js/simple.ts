declare var Chartist

class Simple {

	private	plugin

	selectRecent2D(){
		var select = <HTMLSelectElement>document.getElementById('iterationselector')
		console.log(select.options.length)
		select.options[select.options.length - 1].selected = true
		select.options[select.options.length - 1].click()
		
	}

	viewCavgs(file, status){
		var request = {
			mod : "simple",
			fnc : "getCavgs",
			arg : {}
		}
		request['arg']['file'] = file
		request['arg']['status'] = status
		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				document.getElementById('cavgs').innerHTML = json.html
				viewer.sortThumbnails()
				viewer.loadImages('cavgs')
				simple.sumParticles2D('cavgs')
			})
	}
	
	viewIni3dIteration(file, status){
		var request = {
			mod : "simple",
			fnc : "viewIni3dIteration",
			arg : {}
		}
		request['arg']['file'] = file

		return postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				document.getElementById('inivol').innerHTML = json.html
				viewer.loadImages('reprojections')
				this.plugin =  LiteMol.Plugin.create({ target: '#litemol' })
				var surface = this.plugin.createTransform()
					.add(this.plugin.root, LiteMol.Bootstrap.Entity.Transformer.Data.Download, { url: "/DensityServer/local/" + json['mdb'] + "/cell?detail=4", type: 'Binary', description: 'local Density', title: "joe2", id:"density3ddata"})
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
				LiteMol.Bootstrap.Command.Layout.SetViewportOptions.dispatch(this.plugin.context, { clearColor: LiteMol.Visualization.Color.fromRgb(255, 255, 255) })
				LiteMol.Bootstrap.Command.Visual.ResetScene.dispatch(this.plugin.context, void 0);
			})
	}
	
	updateIsoSurface(sigma) {
		var density = this.plugin.context.select('density3d')[0]
		var style = LiteMol.Bootstrap.Visualization.Density.Style.create({
			isoValue: Number(sigma),
			isoValueType: LiteMol.Bootstrap.Visualization.Density.IsoValueType.Sigma,
			color: LiteMol.Visualization.Color.fromHex(0xBB3333),
			isWireframe: false,
			transparency: { alpha: 1.0 }
		}) 
		LiteMol.Bootstrap.Entity.Transformer.Density.CreateVisual.create({ style: style }, { ref: density.ref, isHidden: false }).update(this.plugin.context, density).run()	
	}
	
	save2DSelection(){
		var thumbnails = document.getElementsByClassName('thumbnail')
		var sorted = []
		for (var thumbnail of thumbnails) {
			var state = 0
			if(thumbnail.getAttribute('data-selected') == "true"){
				state = 1
			}
			sorted.push([Number(thumbnail.getAttribute('data-class')), state])
		}
		sorted.sort((x,y) => { return x[0] - y[0] })
		
		var request = {
			mod : "simple",
			fnc : "save2DSelection",
			arg : {}
		}
		request['arg']['file'] = (<HTMLInputElement>document.getElementById('savefolder')).value + "/" + (<HTMLInputElement>document.getElementById('savefilename')).value
		request['arg']['projectfile'] = (<HTMLInputElement>document.getElementById('saveprojectfile')).value
		request['arg']['selection'] = []
		for (var element of sorted){
			request['arg']['selection'].push(element[1])
		}
		return postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				if(json['status'] == "ok"){
					alert("Selection Saved")
					viewer.toggleScrollMenu('save')
				}else{
					alert("Error Saving Selection")
				}
			})
	}
	
	saveMicsSelection(){
		var thumbnails = document.getElementsByClassName('thumbnail')
		var sorted = []
		for (var thumbnail of thumbnails) {
			var state = 0
			if(thumbnail.getAttribute('data-selected') == "true"){
				state = 1
			}
			sorted.push(state)
		}
		
		var request = {
			mod : "simple",
			fnc : "saveMicsSelection",
			arg : {}
		}
		request['arg']['file'] = (<HTMLInputElement>document.getElementById('savefolder')).value + "/" + (<HTMLInputElement>document.getElementById('savefilename')).value
		request['arg']['projectfile'] = (<HTMLInputElement>document.getElementById('saveprojectfile')).value
		request['arg']['selection'] = sorted
		
		return postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				if(json['status'] == "ok"){
					alert("Selection Saved")
					viewer.toggleScrollMenu('save')
				}else{
					alert("Error Saving Selection")
				}
			})
	}
	
	viewParticles2D(classnumber, projectfile){
		var request = {
			mod : "simple",
			fnc : "viewParticles2D",
			arg : {}
		}
		request['arg']['class'] = classnumber
		request['arg']['projfile'] = projectfile
		var ptclsviewer = document.getElementById('ptclsviewer')
		ptclsviewer.style.display = "block"
		return postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				ptclsviewer.innerHTML = json.html
				viewer.loadImages('ptcls')
				return
			})
		
	}
	
	hideParticles2D(){
		var ptclsviewer = document.getElementById('ptclsviewer')
		ptclsviewer.style.display = "none"
		ptclsviewer.innerHTML = ""
	}
	
	sumParticles2D(parentid){
		var particlecounter =  document.getElementById('particlecounter')
		var thumbnails =  document.getElementById(parentid).getElementsByClassName('thumbnail')
		var total = 0
		var selected = 0
		for (var thumbnail of thumbnails){
			var pop = Number(thumbnail.getAttribute('data-pop'))
			total += pop
			if(thumbnail.getAttribute('data-selected') == "true"){
				selected += pop
			}
		}
		particlecounter.innerHTML = "Selected " + selected + " of " + total + " particles"
	}
	
	plot2D(parentid){
		var thumbnails =  document.getElementById(parentid).getElementsByClassName('thumbnail')
		var plottype = document.getElementById(parentid).querySelector('[id=plottype]') as HTMLInputElement
		var labels = []
		var data = []
		var ymax = (<HTMLInputElement>document.getElementById('ymax')).value
		
		if(plottype.value == "classvpop"){
			for(var thumbnail of thumbnails){
				labels.push(thumbnail.getAttribute('data-class'))
				data.push(thumbnail.getAttribute('data-pop'))
			}
			var chart = new Chartist.Bar('#chart', {
					//labels: labels,
					series: [ data ]
				}, {
					fullWidth: true,
					chartPadding: {
						right: 40
					},
					high:ymax,
					low:0
			});
			chart.on("draw", (data) => {
				if(data.type === "bar"){
					var index = data.index
					var classnumber = thumbnails[index].getAttribute('data-class')
					data.element._node.setAttribute('data-plotclass', classnumber)
					data.element._node.onclick = () => {
						var thumbnail = <HTMLDivElement>document.querySelector('[data-class="' + classnumber + '"]')
						thumbnail.style.backgroundColor = "red"
						thumbnail.scrollIntoView()
						setTimeout(() => { 
							thumbnail.style.backgroundColor = "grey"
						}, 1000)
					}
				}
			})
			
		}else if(plottype.value == "classvres"){
			for(var thumbnail of thumbnails){
				labels.push(thumbnail.getAttribute('data-class'))
				data.push(thumbnail.getAttribute('data-res'))
			}
			var chart = new Chartist.Bar('#chart', {
					//labels: labels,
					series: [ data ]
				}, {
					fullWidth: true,
					chartPadding: {
						right: 40
					},
					high:ymax,
					low:0
					
			});
			chart.on("draw", (data) => {
				if(data.type === "bar"){
					var index = data.index
					var classnumber = thumbnails[index].getAttribute('data-class')
					data.element._node.setAttribute('data-plotclass', classnumber)
					data.element._node.onclick = () => {
						var thumbnail = <HTMLDivElement>document.querySelector('[data-class="' + classnumber + '"]')
						thumbnail.style.backgroundColor = "red"
						thumbnail.scrollIntoView()
						setTimeout(() => { 
							thumbnail.style.backgroundColor = "grey"
						}, 1000)
					}
				}
			})
		}else if(plottype.value == "resvpop"){
			for(var thumbnail of thumbnails){
				data.push({x : thumbnail.getAttribute('data-res'), y : thumbnail.getAttribute('data-pop'), class : "joe"})
			}
			var chart = new Chartist.Line('#chart', {
					series: [ data ]
				}, {
					axisX: {
						type: Chartist.AutoScaleAxis,
						onlyInteger: false
					},
					fullWidth: true,
					showLine: false,
					chartPadding: {
						right: 40
					},
					high:ymax,
					low:0
			});
			chart.on("draw", (data) => {
				if(data.type === "point"){
					var index = data.index
					var classnumber = thumbnails[index].getAttribute('data-class')
					data.element._node.setAttribute('data-plotclass', classnumber)
					data.element._node.onclick = () => {
						var thumbnail = <HTMLDivElement>document.querySelector('[data-class="' + classnumber + '"]')
						thumbnail.style.backgroundColor = "red"
						thumbnail.scrollIntoView()
						setTimeout(() => { 
							thumbnail.style.backgroundColor = "grey"
						}, 1000)
					}
				}
			})
		}
	}
		
	plot2Dptcls(parentid){
		var thumbnails =  document.getElementById(parentid).getElementsByClassName('thumbnail')
		var plottype = document.getElementById(parentid).querySelector('[id=plottype]') as HTMLInputElement
		var plot = plottype.value.split(":")
		var labels = []
		var data = []

		for(var index = 0; index < thumbnails.length; index++){
			if(plot[2] == "scatter"){
				var x
				var y
				if(plot[0] == "ptcl" || plot[0] == "mic"){
					x = index
				}else{
					x = thumbnails[index].getAttribute(`data-${plot[0]}`)
				}
				y = (<HTMLElement>thumbnails[index]).getAttribute(`data-${plot[1]}`)
				data.push({x : x, y : y})
			}
		}
		
		if(plot[2] == "scatter"){
			var chart = new Chartist.Line('#ptcls2Dchart', {
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
					showLine: false,
					chartPadding: {
						right: 40
					}
			})
			chart.on("draw", (data) => {
				if(data.type === "point"){
					var index = data.index
					data.element._node.setAttribute('data-index', index)
					data.element._node.onclick = () => {
						var thumbnail = thumbnails[index] as HTMLElement
					//	var thumbnails = document.getElementById('thumbnails') as HTMLElement
						thumbnail.scrollIntoView()
						thumbnail.style.backgroundColor = "red"
						setTimeout(() => { 
							thumbnail.style.backgroundColor = "grey"
						}, 1000)
					}
				}
			})
		}
	}
	
	lightPlot2Dptcl(ptclid) {
		var plotpoint = document.getElementById('ptcls2Dchart').querySelector(`[data-index='${ptclid}']`)
		if(plotpoint){
			(<HTMLElement>plotpoint).style.stroke = "grey";
			(<HTMLElement>plotpoint).style.strokeWidth = "15px"
			setTimeout(() => { 
				(<HTMLElement>plotpoint).style.stroke = "#d70206";
				(<HTMLElement>plotpoint).style.strokeWidth = "10px"
			}, 1000)
		}
	}
	
	viewMicrograph(path){
		document.getElementById('micrograph').setAttribute('data-path', path)
		document.getElementById('micrograph').setAttribute('data-frame', "0")
		viewer.loadImages('mic')
	}
	
	viewParticles(path, frames){
		var thumbnails =  document.getElementById('ptcls').getElementsByClassName('thumbnails')[0]
		thumbnails.innerHTML = ""
		for(var i = 0; i < frames; i++){
			var thumbnail = document.createElement('div')
			thumbnail.className = "thumbnail"
			var dynimgcontainer = document.createElement('div')
			dynimgcontainer.className = "dynimgcontainer"
			var img = document.createElement('img')
			img.className = "dynimg"
			img.setAttribute('data-path', path)
			img.setAttribute('data-frame', i.toString())
			dynimgcontainer.appendChild(img)
			thumbnail.appendChild(dynimgcontainer)
			thumbnails.appendChild(thumbnail)
		}
		viewer.loadImages('ptcls')
	}
}

var simple
window.addEventListener('load', () => {simple = new Simple()})
