declare var Chartist

class GuinierPlotWidget {

	private	plugin
	private widgetpopup
	
	constructor(){
		this.widgetpopup = document.getElementById('widgetpopup')
	}
	
	refresh(){
		var request = {
			mod : "simple",
			fnc : "getGuinierWidget",
			arg : {}
		}
		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				this.widgetpopup.innerHTML = json.html
				this.widgetpopup.className = "guinierwidgetpopup"
			})
	}
	
	view(element){
		this.refresh()
	}
	
	hide(){
		this.widgetpopup.innerHTML = ""
		this.widgetpopup.className = ""
	}
	
	calculate(){
		var request = {
			mod : "simple",
			fnc : "calculateGuinier",
			arg : {}
		}
		request['arg']['volume'] = (<HTMLInputElement>document.getElementById('guiniervolume')).value
		request['arg']['smpd'] = (<HTMLInputElement>document.getElementById('guiniersmpd')).value
		request['arg']['hp'] = (<HTMLInputElement>document.getElementById('guinierhp')).value
		request['arg']['lp'] = (<HTMLInputElement>document.getElementById('guinierlp')).value
		return postAjaxPromise(request)
			.then((response) => response.json())
			.then ((json) => {
				console.log("B Factor determined to be : ", json.bfac)
				document.getElementById('guinierbfac').innerHTML = "B Factor determined to be : " + json.bfac
				var chart = new Chartist.Line('#guinierchart', {
					series: [ json.plot ]
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
			})
	}


}

var guinierplotwidget
window.addEventListener('load', () => {guinierplotwidget = new GuinierPlotWidget()})
