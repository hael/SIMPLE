class AutoPickTrainingWidget {

	private widgetpopup
	private img
	
	constructor(){
		this.widgetpopup = document.getElementById('widgetpopup')
		this.img = new Image()
	}
	
	
	
	refresh(){
		var request = {
			mod : "simple",
			fnc : "getAutopickTrainingWidget",
			arg : {}
		}
		request['arg']['projfile'] = (<HTMLInputElement>document.getElementById('inputpath')).value
		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				this.widgetpopup.innerHTML = json.html
				this.widgetpopup.className = "autopicktrainingpopup"; 
				(<HTMLInputElement>document.getElementById('autopickrefs')).value = (<HTMLInputElement>document.getElementById('keyrefs')).value
			})
	}
	
	view(element){
		this.refresh()
	}
	
	hide(){
		this.widgetpopup.innerHTML = ""
		this.widgetpopup.className = ""
	}
	
	pick(){
		var request = {}
		request['fnc'] = "trainAutopick"
		request['mod'] = "simple"
		request['arg'] = {}
		request['arg']['projectfolder'] = projectselector.selectedfolder
		request['arg']['projfile'] = (<HTMLInputElement>document.getElementById('inputpath')).value
		request['arg']['micrographid'] = (<HTMLSelectElement>document.getElementById('micrographselector')).value
		request['arg']['refs'] = (<HTMLInputElement>document.getElementById('autopickrefs')).value
		request['arg']['ndev'] = (<HTMLInputElement>document.getElementById('autopickquality')).value
		request['arg']['thres'] = (<HTMLInputElement>document.getElementById('autopickthreshold')).value

		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				this.drawBoxes(json.coordinates, json.boxsize)
			})
	}

	done(){
		
	}
	
	drawBoxes(coordinates, boxsize){
		var canvas = <HTMLCanvasElement>document.getElementById('micrographcanvas')
		var ctx = canvas.getContext('2d')
		ctx.clearRect(0, 0, canvas.width, canvas.height)
		ctx.drawImage(this.img, 0, 0)
		ctx.lineWidth = 4
		ctx.lineJoin = "round"
		ctx.strokeStyle = "magenta"
		for(var coordinate of coordinates){
			ctx.strokeRect(coordinate[0]-(boxsize/2), coordinate[1]-(boxsize/2), boxsize, boxsize)
		}
		
	}
	
	loadMicrograph(path, xdim, ydim){
		var canvas = <HTMLCanvasElement>document.getElementById('micrographcanvas')
		var ctx = canvas.getContext('2d')
		this.img.addEventListener('load',() => {
			ctx.drawImage(this.img,0,0)
		}, false)
		canvas.width = 4096
		canvas.height = 4096
		this.img.src = "/image?stackfile=" + path + "&frame=0&width=" + xdim

	}
}

var autopicktrainingwidget
window.addEventListener('load', () => {autopicktrainingwidget = new AutoPickTrainingWidget()})
