class Highlighter {

  constructor() {
	document.addEventListener('click', event => {
		if(event.shiftKey){
			var borderStyle = "5px solid red"
			if(event.target.style.border == borderStyle){
				event.target.style.border = ""
			}else{
				event.target.style.border = "5px solid red"
			}
		}
	}, false)
  }
  
}

var highlighter
window.addEventListener('load', () => {highlighter = new Highlighter()})
