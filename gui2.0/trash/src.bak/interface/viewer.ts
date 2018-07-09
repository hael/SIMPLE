///<reference path="../../external/LiteMol/dist/js/LiteMol-plugin.d.ts" />

class Viewer {
	
	private viewpane
	private viewpanetitle
	private viewpanebody
	private viewpanesave
	private target

	constructor() {
		this.viewpane = document.getElementById('viewer')
		this.viewpanetitle = document.getElementById('viewertitle')
		this.viewpanebody = document.getElementById('viewerbody')
		this.viewpanesave = document.getElementById('viewersave')
		this.target = document.getElementById('target')
	}
	
	public show(){
		this.viewpanetitle.innerHTML = ""
		this.viewpanebody.innerHTML = ""
		fileselector.show(true)
	}
	
	public showSave(){
		this.viewpanesave.style.display = "block"
	}
	
	public hideTools(){
		document.getElementById('targetdiv').style.display = "none"
		document.getElementById('contrastdiv').style.display = "none"
		document.getElementById('brightnessdiv').style.display = "none"
		document.getElementById('isoleveldiv').style.display = "none"
		document.getElementById('zoomdiv').style.display = "none"
		document.getElementById('savebutton').style.display = "none"
		document.getElementById('scalediv').style.display = "none"
	}
	
	public updateVal(value, elementid){
		var element = document.getElementById(elementid) as HTMLInputElement
		element.value = value
	}
	
	public viewer2DInit() {
		saveselector.sourcefile = fileselector.selection
		var request = {mod : "core", fnc : "view2DInit", arg : { pth : fileselector.selection }}
		postAjax(request, this.viewer2D.bind(this))
		fileselector.hide()
		this.hideTools()
		this.viewpanetitle.innerHTML = "Loading " + fileselector.selection
	}
	
	public viewer2D(data) {
		this.viewpanebody.innerHTML = data
		this.viewpanetitle.innerHTML = "Viewing " + fileselector.selection
		document.getElementById('contrastdiv').style.display = "inline-block"
		document.getElementById('brightnessdiv').style.display = "inline-block"
		document.getElementById('savebutton').style.display = "inline-block"
		document.getElementById('zoomdiv').style.display = "inline-block"
		document.getElementById('scalediv').style.display = "inline-block"
		saveselector.action = viewer.save2D
		var contrast = document.getElementById('contrast') as HTMLInputElement
		var brightness = document.getElementById('brightness') as HTMLInputElement
		var zoom = document.getElementById('zoom') as HTMLInputElement
		var scale = document.getElementById('scale') as HTMLInputElement
		contrast.removeAttribute("onchange")
		brightness.removeAttribute("onchange")
		zoom.removeAttribute("onchange")
		scale.removeAttribute("onchange")
		contrast.onchange = function() {viewer.modifyImages()}
		brightness.onchange = function() {viewer.modifyImages()}
		zoom.onchange = function() {viewer.modifyImages()}
		scale.onchange = function() {viewer.scaleImages()}
		this.loadImages()
	}
		
	public save2D(){
		var savefolder = document.getElementById('savefolder') as HTMLInputElement
		var savefilename = document.getElementById('savefilename') as HTMLInputElement
		var request = {}
		request['mod'] = "core"
		request['fnc'] = "save2D"
		request['arg'] = {}
		request['arg']['pth'] = savefolder.value + "/" + savefilename.value
		request['arg']['src'] = saveselector.sourcefile
		request['arg']['selection'] = []
		var selection = document.querySelectorAll("[data-selected='true']") as NodeListOf<HTMLDivElement>
		for(var selected of selection){
			request['arg']['selection'].push(selected.dataset.identity)
		}
		return postAjaxPromise(request)
	}
		
	public viewer3DInit() {
		var request = {mod : "core", fnc : "createMDB", arg : { pth : fileselector.selection }}
		postAjax(request, this.viewer3DConfirm.bind(this))
		fileselector.hide()
		this.viewpanetitle.innerHTML = "Loading " + fileselector.selection
	}
	
	public viewer3DConfirm(data){
		var dataJSON = JSON.parse(data)
		var	trialcount = 0
		var waiter = setInterval(function() {
			fetch("/DensityServer/local/" + dataJSON.pth + "/").then(function(status) {
				if(status.ok){
					clearInterval(waiter)
					viewer.viewer3D(data)
				}
			})
			trialcount++
			if(trialcount == 10){
				clearInterval(waiter)
				alert("Error")
			}
		}, 500)		
	}
	
	public viewer3D(data) {
		var dataJSON = JSON.parse(data)
		
		//document.body.removeChild(document.getElementById('fileselector'))
		this.viewpanebody.innerHTML = ""
		
		this.hideTools()
		document.getElementById('isoleveldiv').style.display = "inline-block"
		
		var isolevel = document.getElementById('isolevel') as HTMLInputElement
		isolevel.removeAttribute("onchange");

		isolevel.onchange = function(){
			var density = plugin.context.select('density3d')[0]
			var style = LiteMol.Bootstrap.Visualization.Density.Style.create({
				isoValue: Number(isolevel.value),
				isoValueType: LiteMol.Bootstrap.Visualization.Density.IsoValueType.Sigma,
				color: LiteMol.Visualization.Color.fromHex(0xBB3333),
				isWireframe: false,
				transparency: { alpha: 1.0 }
			}) 
			//Required//LiteMol.Bootstrap.Entity.Transformer.Density.CreateVisual.create({ style: style }, { ref: density.ref, isHidden: false }).update(plugin.context, density).run()	
		}
		
		var plugin = LiteMol.Plugin.create({ target: '#viewerbody' })
		
		this.viewpanetitle.innerHTML = "Viewing " + fileselector.selection
		
		var surface = plugin.createTransform()
			.add(plugin.root, LiteMol.Bootstrap.Entity.Transformer.Data.Download, { url: "/DensityServer/local/" + dataJSON.pth + "/cell?detail=4", type: 'Binary', description: 'local Density', title: "joe2", id:"density3ddata"})
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
					}, {ref : "density3d"});
		plugin.applyTransform(surface);
		LiteMol.Bootstrap.Command.Layout.SetViewportOptions.dispatch(plugin.context, { clearColor: LiteMol.Visualization.Color.fromRgb(255, 255, 255) });
		LiteMol.Bootstrap.Command.Visual.ResetScene.dispatch(plugin.context, void 0)
	}
	
	public loadImages() {
		var images = document.getElementsByClassName('dynimg') as HTMLCollectionOf<HTMLImageElement>
		for (var image of images){
			image.src = "/image?mod=core&fnc=mrc2JPEG&pth=" + image.dataset.path + "&frm=" + image.dataset.frame + "&wth=" +image.clientWidth
		}
	}
	
	public modifyImages() {
		var brightness = document.getElementById('brightness') as HTMLInputElement
		var contrast = document.getElementById('contrast') as HTMLInputElement
		var zoom = document.getElementById('zoom') as HTMLInputElement
		
		var images = document.getElementsByClassName('dynimg') as HTMLCollectionOf<HTMLImageElement>
		for (var image of images){
			image.style.filter = "contrast(" + contrast.value + "%) brightness(" + brightness.value + "%)"
			image.style.transform = "scale(" + zoom.value + ")"
		}
	}
	
	public scaleImages() {
		var scale = document.getElementById('scale') as HTMLInputElement
		var containers = document.getElementsByClassName('dynimgcontainer') as HTMLCollectionOf<HTMLDivElement>
		for (var container of containers){

			container.style.width = scale.value + "px"
			container.style.height = scale.value + "px"
		}
		this.loadImages()
	}
	
	public toggleThumbnail(thumbnail){
		var images = thumbnail.getElementsByTagName('img')
		if(thumbnail.dataset.selected == "true"){
			thumbnail.dataset.selected = false
			for(var image of images){
				image.style.visibility = "hidden"
			}
		}else{
			thumbnail.dataset.selected = true
			for(var image of images){
				image.style.visibility = "visible"
			}
		}
	}
	
}

var viewer = new Viewer()
