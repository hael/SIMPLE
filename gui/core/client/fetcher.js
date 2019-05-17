class Fetcher{
  constructor() {}

  fetchJSON (request) {
    return fetch("/", {
      method: 'POST',
      body: JSON.stringify(request),
      headers: new Headers({
        'Content-Type': 'application/json',
      }),
      credentials: 'include'
    })
  }
  
  fetchImages() {
    var images = document.getElementsByClassName('thumbnailcontainerimg')
    var blobjs = `var loadQueue = []
        self.onmessage = (msg) => {
          var sources = msg.data;
          for(var source in sources){
            loadQueue.push([sources[source], source])
          }
          processQueue();
          //processQueue();
          //processQueue();
        };

        function processQueue(){
          var queueElement = loadQueue[0]
          loadQueue.shift()
          return fetch('` + window.location.origin  + `' + queueElement[0], {credentials: 'include'})
            .then(() => {
              postMessage({id : queueElement[1], source : queueElement[0]})
              return
            })
            .then(() => {
              processQueue()
            })
        };
    `;

    const workerblob = new Blob([blobjs]);
    this.loaderworker = new Worker(window.URL.createObjectURL(workerblob));
    this.loaderworker.onmessage = (event) =>{
      var image = images[Number(event.data['id'])]
      if(image.dataset.path){
        image.src = event.data['source']
      } else if(image.dataset.sprite){
        image.style.background = "url(" + event.data['source'] + ")" + Number(image.clientWidth) * Number(image.dataset.spriteid)  + "px 0px"
      }
     image.parentElement.getElementsByClassName('miniloader')[0].style.display = "none"
    }

    var sources = []
    for (var image of images){
	  if(image.dataset.path && image.dataset.thumb && image.dataset.boxfile){
		sources.push("/image?stackfile=" + image.dataset.path + "&boxfile=" + image.dataset.boxfile + "&width=" + image.clientWidth + "&thumb=" + image.dataset.thumb)
	  }else if(image.dataset.path){
        sources.push(image.dataset.path + "&width=" + this.imagewidth)
      }
      //}else if(image.dataset.sprite){
     //   sources.push("/image?stackfile=" + image.dataset.sprite + "&frame=0&width=" + Number(image.clientWidth) * Number(image.dataset.spritewidth))
     // }
    }

    this.loaderworker.postMessage(sources)
  }
  

  
  killWorker(){
	if(this.loaderworker != undefined){
		this.loaderworker.terminate()
		this.loaderworker = undefined
	}
  }
}

let fetcher
window.addEventListener('load', () => {fetcher = new Fetcher()})
