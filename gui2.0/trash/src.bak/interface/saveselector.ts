class SaveSelector {
	
	public action
	public sourcefile
	private selector
	private gauze
	
	constructor(){
		this.selector = document.getElementById('saveselector')
		this.gauze = document.getElementById('gauze')
	}
	
	public show(){
		this.gauze.style.display = "block"
		this.selector.style.display = "block"
	}
	
	public hide(){
		this.selector.style.display = "none"
		this.gauze.style.display = "none"
	}
	
	public save(){
		var saveselection = document.getElementById('saveselection')
		saveselection.innerHTML = "Saving ..."
		this.action().then(function(){saveselection.innerHTML = "Save"; saveselector.hide(); alert("Saved")})

	}

}

var saveselector = new SaveSelector()
