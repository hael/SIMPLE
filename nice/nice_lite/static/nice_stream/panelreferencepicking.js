let lastinteraction = Date.now();

scrlRight = () => {
  const micrograph_slider = document.getElementById("picking_references")
  micrograph_slider.scrollLeft += 82;
  lastinteraction = Date.now();
}

scrlLeft = () => {
  const micrograph_slider = document.getElementById("picking_references")
  micrograph_slider.scrollLeft -= 82;
  lastinteraction = Date.now();
}

scrlMicRight = () => {
  const micrograph_slider = document.getElementById("micrographs_slider")
  micrograph_slider.scrollLeft += 200;
  lastinteraction = Date.now();
}

scrlMicLeft = () => {
  const micrograph_slider = document.getElementById("micrographs_slider")
  micrograph_slider.scrollLeft -= 200;
  lastinteraction = Date.now();
}

var pick_observer = new IntersectionObserver(
  (entries, opts) => {
    entries.forEach(entry => { 
      if(entry.intersectionRatio == 1){
        const boxes_overlay = document.getElementById("boxes_overlay")
        if(boxes_overlay ==  null) return
        boxes_overlay.dataset.xdim  = entry.target.dataset.xdim
        boxes_overlay.dataset.ydim  = entry.target.dataset.ydim
        boxes_overlay.dataset.boxes = entry.target.dataset.boxes
        draw_overlay_coordinates()
      }
    })
  }, 
  {
    root: document.getElementById('micrographs_slider'), 
    threshold: .5 
  }
)

draw_overlay_coordinates = () => {
    let multipick = true
    const boxes_overlay = document.getElementById("boxes_overlay")
    if(boxes_overlay ==  null) return
    const scale = boxes_overlay.width  / Number(boxes_overlay.dataset.xdim)
    // account for rectangular images
    var yoffset = (boxes_overlay.height - (scale * Number(boxes_overlay.dataset.ydim))) / 2
    const boxes = JSON.parse(boxes_overlay.dataset.boxes.replaceAll("'", '"'))
    const ctx = boxes_overlay.getContext("2d");
    ctx.clearRect(0, 0, boxes_overlay.width, boxes_overlay.height)
    for(const box of boxes){
      ctx.strokeStyle = "yellow";
      ctx.beginPath();
      ctx.arc(box["x"] * scale, (box["y"] * scale) + yoffset, 0.5, 0, 2 * Math.PI);
      ctx.stroke();
    }
}

showMenu = (element, event) => {
  event.preventDefault()
  const selectmenu    = element.parentElement.parentElement.querySelector("[name='selectmenu']")
  const selectmenubox = selectmenu.querySelector("[name='selectmenubox']")
  selectmenubox.style.top  = event.pageY + "px"
  selectmenubox.style.left = event.pageX + "px"
  selectmenu.style.display = "flex"
  lastinteraction = Date.now();
}

hideMenu = () => {
  for(const selectmenu of document.querySelectorAll("[name='selectmenu']")){
    selectmenu.style.display = "none"
  }
  lastinteraction = Date.now();
}

updateBrightness = (element) => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--preprocess-brightness', element.value / 100);
}

updateContrast = (element) => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--preprocess-contrast', element.value / 100);
}

window.addEventListener("load", () => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--preprocess-contrast',   1.0);
  cssroot.style.setProperty('--preprocess-brightness', 1.0);
})

/* draw boxes on load */
window.addEventListener("load", () =>{
    for(const box_overlay of document.getElementsByClassName("boxes_overlay")){
      const scale = box_overlay.width  / Number(box_overlay.dataset.xdim)
      // account for rectangular images
      var yoffset = (box_overlay.height - (scale * Number(box_overlay.dataset.ydim))) / 2
      const boxes = JSON.parse(box_overlay.dataset.boxes.replaceAll("'", '"'))
      const ctx = box_overlay.getContext("2d");
      for(const box of boxes){
        ctx.strokeStyle = "yellow";
        ctx.beginPath();
        ctx.arc(box["x"] * scale, (box["y"] * scale) + yoffset, 0.5, 0, 2 * Math.PI);
        ctx.stroke();
      }
    }
},false);

window.addEventListener("load", () =>{
    for(const micrographs_pie_chart of document.getElementsByClassName("micrographs_pie_chart")){
        const ctx = micrographs_pie_chart.getContext("2d");
        const n_imported  = Number(micrographs_pie_chart.dataset.imported) 
        const n_processed = Number(micrographs_pie_chart.dataset.processed)
        const n_rejected  = Number(micrographs_pie_chart.dataset.rejected)
        new Chart(ctx, {
            type: 'doughnut',
            options:{
              maintainAspectRatio : false,
              plugins:{
                legend:{
                    position : "right",
                    labels:{
                      boxWidth: 10,
                      padding:  2,
                      font :{
                        size: 9
                      }
                    }
                }
              }
            },
            data: {
              labels: [
                  'queued',
                  'processed',
                  'rejected'
              ],
              datasets: [{
                  data: [n_imported - n_processed - n_rejected, n_processed, n_rejected],
                  backgroundColor: [
                  'rgb(255, 99, 132)',
                  'rgb(54, 162, 235)',
                  'rgb(255, 205, 86)'
                  ],
                  hoverOffset: 4
              }]
            }
        })
    }
},false);

window.addEventListener("load", () =>{
    for(const slider_micrograph of document.getElementsByClassName('slider_micrograph')){
      pick_observer.observe(slider_micrograph)
    }
},false);

window.addEventListener("load", () =>{
  document.getElementById("loadinggauze").style.opacity = "0";
  setTimeout(function () {
   document.getElementById("loadinggauze").style.display = "none";
  }, 600);
})

window.addEventListener("visibilitychange", (event) => {
  if(document.visibilityState !== "hidden"){
    location.reload();
  }
})

setInterval(function () {
  if((Date.now() - lastinteraction) > 30000 && document.visibilityState !== "hidden"){
    lastinteraction = Date.now();
    location.reload();
  }
}, 1000);