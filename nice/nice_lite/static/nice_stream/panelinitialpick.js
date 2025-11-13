let lastinteraction = Date.now();

restartProcess = (element) => {
  const confirmed = confirm("Please confirm that you wish to restart this process");
  if(confirmed){
    element.form.submit()
  }
}

stopProcess = (element) => {
  const confirmed = confirm("Please confirm that you wish to stop this process");
  if(confirmed){
    element.form.submit()
  }
}

scrlRight = (element, event) => {
  event.preventDefault()
  const micrograph_slider = document.getElementById("micrographs_slider")
  const rect = micrograph_slider.getBoundingClientRect();
  micrograph_slider.scrollLeft += rect.width;
  lastinteraction = Date.now();
}

scrlLeft = (element, event) => {
  event.preventDefault()
  const micrograph_slider = document.getElementById("micrographs_slider")
  const rect = micrograph_slider.getBoundingClientRect();
  micrograph_slider.scrollLeft -= rect.width;
  lastinteraction = Date.now();
}

refineDiameter = () => {
  const refinediameter  = document.querySelector('#refinediameter')
  const refine_diameter = refinediameter.querySelector("[name='refine_diameter']")
  const diameter_input  = document.getElementsByName("diameter")[0]
  refine_diameter.value = diameter_input.value
  refinediameter.submit()
  return false
}

var multipick_observer = new IntersectionObserver(
  (entries, opts) => {
    entries.forEach(entry => { 
      if(entry.intersectionRatio > 0.5){
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

    const closest = (arr, n) => arr.reduce((prev, curr) => Math.abs(curr - n) < Math.abs(prev - n) ? curr : prev);

    let multipick = true
    let best_diam = 0
    const boxes_overlay = document.getElementById("boxes_overlay")
    if(boxes_overlay ==  null) return
    const picking_diameters = document.getElementById("picking_diameters")
    if(picking_diameters ==  null) multipick = false
    const diameter_input = document.getElementsByName("diameter")[0]
    if(diameter_input ==  null) multipick = false
    const micrographs_slider = document.querySelector("#micrographs_slider")
    if(micrographs_slider != undefined && "best_diam" in micrographs_slider.dataset && "picking_diameters" in micrographs_slider.dataset){
      const available_diameters = JSON.parse(micrographs_slider.dataset.picking_diameters.replaceAll("'", '"'))
      multipick = true
      best_diam = closest(available_diameters, Number(micrographs_slider.dataset.best_diam))
    }
    const xdim = Number(boxes_overlay.dataset.xdim)
    const ydim = Number(boxes_overlay.dataset.ydim)
    let scale  = 1
    if(xdim > ydim){
      scale = boxes_overlay.width  / Number(boxes_overlay.dataset.xdim)
    }else{
      scale = boxes_overlay.height  / Number(boxes_overlay.dataset.ydim)
    }
    // account for rectangular images
    var yoffset = (boxes_overlay.height - (scale * Number(boxes_overlay.dataset.ydim))) / 2
    var xoffset = (boxes_overlay.width  - (scale * Number(boxes_overlay.dataset.xdim))) / 2
    const boxes = JSON.parse(boxes_overlay.dataset.boxes.replaceAll("'", '"'))
    const ctx = boxes_overlay.getContext("2d");
    ctx.clearRect(0, 0, boxes_overlay.width, boxes_overlay.height)
    for(const box of boxes){
      ctx.strokeStyle = "yellow";
      ctx.beginPath();
      ctx.arc((box["x"] * scale) + xoffset, (box["y"] * scale) + yoffset, 1, 0, 2 * Math.PI);
      ctx.stroke();
    }
    lastinteraction = Date.now();
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
  cssroot.style.setProperty('--inipick-brightness', element.value / 100);
}

updateContrast = (element) => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--inipick-contrast', element.value / 100);
}

updateScale = (element) => {
  const micrographs_slider = document.querySelector("#micrographs_slider")
  const boxes_overlay      = document.querySelector("#boxes_overlay")
  if(micrographs_slider != undefined) micrographs_slider.style.width = element.value + "px"
  if(boxes_overlay      != undefined) {
    boxes_overlay.height = element.value
    boxes_overlay.width  = element.value
  }
  for(const miccontainer of document.querySelectorAll(".miccontainer")){
    miccontainer.style.height = element.value + "px"
    miccontainer.style.width  = element.value + "px"
    const img = miccontainer.querySelector("img")
    if(img != undefined){
      img.style.height = element.value + "px"
      img.style.width  = element.value + "px"
    }
  }
  draw_overlay_coordinates()
}

window.addEventListener("load", () => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--inipick-contrast',   1.0);
  cssroot.style.setProperty('--inipick-brightness', 0.5);
})

/* draw boxes on load */
window.addEventListener("load", () =>{
    for(const slider_micrograph of document.getElementsByClassName('slider_micrograph')){
      multipick_observer.observe(slider_micrograph)
    }
},false);

window.addEventListener("load", () => {
  const logtext = document.querySelector(".logtext")
  logtext.scrollTop = logtext.scrollHeight - logtext.clientHeight
})

window.addEventListener("load", () =>{

    const closest = (arr, n) => arr.reduce((prev, curr) => Math.abs(curr - n) < Math.abs(prev - n) ? curr : prev);
    
    const micrographs_slider = document.querySelector("#micrographs_slider")
    if(micrographs_slider != undefined && "best_diam" in micrographs_slider.dataset && "picking_diameters" in micrographs_slider.dataset){
      const available_diameters = JSON.parse(micrographs_slider.dataset.picking_diameters.replaceAll("'", '"'))
      const best_diam = closest(available_diameters, Number(micrographs_slider.dataset.best_diam))
      const diameter_input = document.getElementsByName("diameter")[0]
      if(diameter_input != undefined) diameter_input.value = available_diameters.indexOf(best_diam)
    }
},false);

/* draw boxes on load */
window.addEventListener("load", () =>{
    for(const miccontainer of document.getElementsByClassName("miccontainer")){
      const box_overlay = miccontainer.querySelector('.boxes-overlay')
      const xdim = Number(miccontainer.dataset.xdim)
      const ydim = Number(miccontainer.dataset.ydim)
      let scale  = 1
      if(xdim > ydim){
        scale = box_overlay.width  / Number(miccontainer.dataset.xdim)
      }else{
        scale = box_overlay.height  / Number(miccontainer.dataset.ydim)
      }
      // account for rectangular images
      var yoffset = (box_overlay.height - (scale * Number(miccontainer.dataset.ydim))) / 2
      var xoffset = (box_overlay.width  - (scale * Number(miccontainer.dataset.xdim))) / 2
      const boxes = JSON.parse(miccontainer.dataset.boxes.replaceAll("'", '"'))
      const ctx = box_overlay.getContext("2d");
      for(const box of boxes){
        ctx.strokeStyle = "yellow";
        ctx.beginPath();
        ctx.arc((box["x"] * scale) + xoffset, (box["y"] * scale) + yoffset, 1, 0, 2 * Math.PI);
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
                    window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                    window.getComputedStyle(document.body).getPropertyValue('--color-nice4success'),
                    window.getComputedStyle(document.body).getPropertyValue('--color-nice4alert'),
                  ],
                  hoverOffset: 4,
                  borderColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble')
              }]
            }
        })
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

/*window.addEventListener("load", (event) => {
  // reload generate refs
  console.log("RELOAD GEPICKREFS")
  window.parent.document.querySelector("#genpickrefs_iframe").contentWindow.location.reload()
})*/

setInterval(function () {
  if((Date.now() - lastinteraction) > 30000 && document.visibilityState !== "hidden"){
    lastinteraction = Date.now();
    location.reload();
  }
}, 1000);