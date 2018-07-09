class Simple {
	
	private viewpane
	private viewpanetitle
	public viewpanebody

	constructor() {
		this.viewpane = document.getElementById('viewer')
		this.viewpanetitle = document.getElementById('viewertitle')
		this.viewpanebody = document.getElementById('viewerbody')
	}
	
	public show(){
		this.viewpanetitle.innerHTML = "Microscope Sync"
		this.viewpanebody.innerHTML = ""
		this.request()
	}
	
	public request () {
		var request = {mod : "sync", fnc : "display", arg : { viw : "syncdisplay" }}
		postAjax(request, this.draw)
	}
	
	public draw (data) {
		document.getElementById('viewerbody').innerHTML = data
	}
}

var simple = new Simple()
