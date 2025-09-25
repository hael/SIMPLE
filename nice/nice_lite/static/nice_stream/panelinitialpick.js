let lastinteraction = Date.now();

scrlRight = () => {
  const micrograph_slider = document.getElementById("micrographs_slider")
  micrograph_slider.scrollLeft += 200;
  lastinteraction = Date.now();
}

scrlLeft = () => {
  const micrograph_slider = document.getElementById("micrographs_slider")
  micrograph_slider.scrollLeft -= 200;
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
    const picking_diameters = document.getElementById("picking_diameters")
    if(picking_diameters ==  null) multipick = false
    const diameter_input = document.getElementsByName("diameter")[0]
    if(diameter_input ==  null) multipick = false
    const scale = boxes_overlay.width  / Number(boxes_overlay.dataset.xdim)
    // account for rectangular images
    var yoffset = (boxes_overlay.height - (scale * Number(boxes_overlay.dataset.ydim))) / 2
    const boxes = JSON.parse(boxes_overlay.dataset.boxes.replaceAll("'", '"'))
    const ctx = boxes_overlay.getContext("2d");
    ctx.clearRect(0, 0, boxes_overlay.width, boxes_overlay.height)
    if(multipick){
      const available_diameters = JSON.parse(picking_diameters.dataset.diameters.replaceAll("'", '"'))
      const selected_diameter = available_diameters[picking_diameters.value - 1]
      diameter_input.value = selected_diameter
      for(const box of boxes){
        if(box["diameter"] == selected_diameter){
          ctx.strokeStyle = "yellow";
          ctx.beginPath();
          ctx.arc(box["x"] * scale, (box["y"] * scale) + yoffset, 0.5, 0, 2 * Math.PI);
          ctx.stroke();
        }
      }
    }else{
      for(const box of boxes){
        ctx.strokeStyle = "yellow";
        ctx.beginPath();
        ctx.arc(box["x"] * scale, (box["y"] * scale) + yoffset, 0.5, 0, 2 * Math.PI);
        ctx.stroke();
      }
    }
    lastinteraction = Date.now();
}

/* draw boxes on load */
window.addEventListener("load", () =>{
    for(const slider_micrograph of document.getElementsByClassName('slider_micrograph')){
      multipick_observer.observe(slider_micrograph)
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