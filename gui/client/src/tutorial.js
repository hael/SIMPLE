class Tutorial {

  constructor() {
  }

  getHTML() {
    for(var input of document.getElementsByTagName('input')){
      input.setAttribute("value", input.value)
    }
    for(var textarea of document.getElementsByTagName('textarea')){
     textarea.setAttribute("value", input.value)
    }
    for(var image of document.getElementsByTagName('img')){
      if(image.src.includes("/image")){
        this.replaceBase64Image(image)
      }
    }
   return(document.documentElement.innerHTML)
  }

  replaceBase64Image(img) {
    var canvas = document.createElement("canvas");
    canvas.width = img.width;
    canvas.height = img.height;
    var ctx = canvas.getContext("2d");
    ctx.drawImage(img, 0, 0);
    var dataURL = canvas.toDataURL("image/png");
    img.src = dataURL;
  }
}

var tutorial
window.addEventListener('load', () => {tutorial = new Tutorial()})
