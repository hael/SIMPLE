class HeaderButtons {
	
	private hbelement
	
	constructor(){
		this.hbelement = document.getElementById('headerbuttons')
		var request = {
			mod : "core",
			fnc : "headerButtons",
			arg : {}
		}
		postAjaxPromise(request)
			.then(response => response.json())
			.then ((json) => {
				this.hbelement.innerHTML = json.html
			})
		
	}
	
}

var headerbuttons
window.addEventListener('load', () => {headerbuttons = new HeaderButtons()})
