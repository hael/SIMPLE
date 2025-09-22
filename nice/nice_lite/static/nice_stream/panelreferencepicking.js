scrlRight = () => {
  const micrograph_slider = document.getElementById("picking_references")
  micrograph_slider.scrollLeft += 82;
}

scrlLeft = () => {
  const micrograph_slider = document.getElementById("picking_references")
  micrograph_slider.scrollLeft -= 82;
}

scrlMicRight = () => {
  const micrograph_slider = document.getElementById("micrographs_slider")
  micrograph_slider.scrollLeft += 200;
}

scrlMicLeft = () => {
  const micrograph_slider = document.getElementById("micrographs_slider")
  micrograph_slider.scrollLeft -= 200;
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

setTimeout(function () {
  document.getElementById("loadinggauze").style.display = "flex";
  document.getElementById("loadinggauze").style.opacity = "1";
  setTimeout(function () {
   location.reload();
  }, 600);
}, 31000);