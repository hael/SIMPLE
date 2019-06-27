class Loader {
  
  constructor() {
    setTimeout(() => {
      document.getElementById('loader').style.opacity = 0
    }, 1000)
    setTimeout(() => {
      document.getElementById('loader').style.display = 'none'
    }, 1500)
    
  }

}

let loader
window.addEventListener('load', () => {loader = new Loader()})
