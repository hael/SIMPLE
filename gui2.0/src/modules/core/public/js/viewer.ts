///<reference path="../../../../../external/LiteMol/js/LiteMol-plugin.d.ts" />

class Viewer {
	
	private popup
	private	gauze
	private plugin
	
	constructor() {
		this.popup = document.getElementById('browserpopup')
		this.gauze = document.getElementById('gauze')
	}
	
	show(html) {
		this.gauze.style.display = "block"
		this.popup.innerHTML = html
		this.popup.className = "viewer"
	}
	
	hide() {
		this.popup.innerHTML = ""
		this.popup.className = "popup"
		this.gauze.style.display = "none"
	}
	
	loadImages(element) {
		var images = document.getElementById(element).getElementsByClassName('dynimg') 
		for (var image of images){
			(<HTMLImageElement>image).src = "/image?stackfile=" + (<HTMLImageElement>image).dataset.path + "&frame=" + (<HTMLImageElement>image).dataset.frame + "&width=" + (<HTMLImageElement>image).clientWidth
		}
	}
	
	modifyImages(element) {
		var brightness = document.getElementById(element).querySelector('[id=brightness]') as HTMLInputElement
		var contrast = document.getElementById(element).querySelector('[id=contrast]') as HTMLInputElement
		var zoom = document.getElementById(element).querySelector('[id=zoom]') as HTMLInputElement
		
		var images =  document.getElementById(element).getElementsByClassName('dynimg') 
		for (var image of images){
			(<HTMLImageElement>image).style.filter = "contrast(" + contrast.value + "%) brightness(" + brightness.value + "%)";
			(<HTMLImageElement>image).style.transform = "scale(" + zoom.value + ")";
		}
	}
	
	scaleImages(element) {
		var scale = document.getElementById(element).querySelector('[id=scale]') as HTMLInputElement
		var containers = document.getElementById(element).getElementsByClassName('dynimgcontainer') 
		for (var container of containers){
			(<HTMLDivElement>container).style.width = scale.value + "px";
			(<HTMLDivElement>container).style.height = scale.value + "px";
		}
		this.loadImages(element)
	}
	
	view2d() {
		var request = {
			mod : "core",
			fnc : "view2D",
			arg : {}
		}
		request['arg']['file'] = browser.selection
		return postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				viewer.show(json.html)
				viewer.loadImages('frames')
			})
	}
	
	view3d() {
		var request = {
			mod : "core",
			fnc : "view3D",
			arg : {}
		}
		
		request['arg']['file'] = browser.selection
		return postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				viewer.show(json.html)
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
				LiteMol.Bootstrap.Command.Layout.SetViewportOptions.dispatch(this.plugin.context, { clearColor: LiteMol.Visualization.Color.fromRgb(255, 255, 255) });
				LiteMol.Bootstrap.Command.Visual.ResetScene.dispatch(this.plugin.context, void 0)
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
	
	toggleSelect(element) {
		var images = element.getElementsByClassName('dynimg')
		var plotpoint = document.querySelector('[data-plotclass="' + element.getAttribute('data-class') + '"]')
		for (var image of images){
			if(element.dataset.selected == "true"){
				element.dataset.selected = "false"
				image.style.visibility = "hidden"
				if(plotpoint){
					(<HTMLElement>plotpoint).style.stroke = "grey"
				}
			}else{
				element.dataset.selected = "true"
				image.style.visibility = "unset"
				if(plotpoint){
					(<HTMLElement>plotpoint).style.stroke = "#d70206"
				}
			}
		}
	}
	
	sortThumbnails() {
		var attribute = (<HTMLInputElement>document.getElementById('sortattribute')).value
		var order = (<HTMLInputElement>document.getElementById('sortorder')).value
		var thumbnails = document.getElementsByClassName('thumbnail')
		var container = thumbnails[0].parentElement
		var sorted = []
		for (var thumbnail of thumbnails){
			sorted.push([thumbnail.getAttribute('data-' + attribute), thumbnail])
		}

		for (var i = thumbnails.length - 1 ; i >= 0; i--){
			thumbnails[i].remove()
		}
		if(order == "ascending"){
			sorted.sort(function(x,y){ return x[0] - y[0] })
		}else if(order == "descending"){
			sorted.sort(function(x,y){ return y[0] - x[0] })
		}
		
		for(var element of sorted){
			container.appendChild(element[1])
		}	
		
	}
	
	selectThumbnails() {
		var attribute = (<HTMLInputElement>document.getElementById('selectattribute')).value
		var operator = (<HTMLInputElement>document.getElementById('selectoperator')).value
		var value = Number((<HTMLInputElement>document.getElementById('selectvalue')).value)
		var thumbnails = document.getElementsByClassName('thumbnail')

		for (var thumbnail of thumbnails){
			var thumbvalue = Number(thumbnail.getAttribute('data-' + attribute))
			var images = thumbnail.getElementsByClassName('dynimg')
			if((<HTMLInputElement>document.getElementById('selectappend')).checked === false){
				thumbnail.setAttribute('data-selected', "true")
				for (var image of images){
					(<HTMLElement>image).style.visibility = "unset"	
				}
			}
			if(operator == "gt" && thumbvalue <= value){
				(<HTMLElement>thumbnail).dataset.selected = "false"
				for (var image of images){
					(<HTMLElement>image).style.visibility = "hidden"	
				}
			}else if (operator == "lt" && thumbvalue >= value){
				(<HTMLElement>thumbnail).dataset.selected = "false"
				for (var image of images){
					(<HTMLElement>image).style.visibility = "hidden"	
				}
			}else if (operator == "eq" && thumbvalue == value){
				(<HTMLElement>thumbnail).dataset.selected = "false"
				for (var image of images){
					(<HTMLElement>image).style.visibility = "hidden"	
				}
			}
		}

		
	}
	
	toggleScrollMenu(parentname, menuname) {
		var menus = document.getElementById(parentname).getElementsByClassName('scrollmenu')
		var selectedmenu = document.getElementById(parentname).querySelector("[id=" + menuname + "]") as HTMLElement
		var current = selectedmenu.style.display
		for (var menu of menus) {
			(<HTMLElement>menu).style.display = "none"
		}
		if(current != "block"){
			selectedmenu.style.display = "block"
		}
	}
		
}

var viewer
window.addEventListener('load', () => {viewer = new Viewer()})
