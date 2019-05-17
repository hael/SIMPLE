class GuinierPlotWidget {
  
  constructor(){
  }

  refresh(){
    var request = {
      cls : "widget",
      fnc : "getGuinierWidget",
      arg : {
		file : document.getElementById('inputpath').value
	  }
    }
    return fetcher.fetchJSON(request)
    .then((response) => response.json())
    .then ((json) => {
      document.getElementById('widgetpopup').style.display = "block"
      document.getElementById('widgetwindow').className = "guinierwidgetpopup"
      document.getElementById('widgetwindow').innerHTML = json.html
    })
  }

  view(element){
    this.refresh()
  }

  hide(){
    document.getElementById('widgetpopup').style.display = "none"
    document.getElementById('widgetwindow').className = ""
    document.getElementById('widgetwindow').innerHTML = ""
  }
  
  use(){
	document.getElementById('keybfac').value = this.bfac
	this.hide()
  }

  calculate(){
	var guiniervolume = document.getElementById('guiniervolume')
    var request = {
      cls : "widget",
      fnc : "calculateGuinier",
      arg : {
        volume : guiniervolume.options[guiniervolume.selectedIndex].dataset.vol,
        smpd : guiniervolume.options[guiniervolume.selectedIndex].dataset.smpd,
        hp : document.getElementById('guinierhp').value,
        lp : document.getElementById('guinierlp').value,
        projfile:document.getElementById('inputpath').value
      }
    }
    return fetcher.fetchJSON(request)
      .then(response => response.json())
      .then (json => {
		this.bfac = json.bfac
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
