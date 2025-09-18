function navigation(slider) {
  let wrapper, dots, arrowLeft, arrowRight

  function markup(remove) {
    wrapperMarkup(remove)
    dotMarkup(remove)
    arrowMarkup(remove)
  }

  function removeElement(elment) {
    elment.parentNode.removeChild(elment)
  }
  function createDiv(className) {
    var div = document.createElement("div")
    var classNames = className.split(" ")
    classNames.forEach((name) => div.classList.add(name))
    return div
  }

  function arrowMarkup(remove) {
    if (remove) {
      removeElement(arrowLeft)
      removeElement(arrowRight)
      return
    }
    arrowLeft = createDiv("arrow arrow--left")
    arrowLeft.addEventListener("click", () => slider.prev())
    arrowRight = createDiv("arrow arrow--right")
    arrowRight.addEventListener("click", () => slider.next())

    wrapper.appendChild(arrowLeft)
    wrapper.appendChild(arrowRight)
  }

  function wrapperMarkup(remove) {
    if (remove) {
      var parent = wrapper.parentNode
      while (wrapper.firstChild)
        parent.insertBefore(wrapper.firstChild, wrapper)
      removeElement(wrapper)
      return
    }
    wrapper = createDiv("navigation-wrapper")
    slider.container.parentNode.appendChild(wrapper)
    wrapper.appendChild(slider.container)
  }

  function dotMarkup(remove) {
    if (remove) {
      removeElement(dots)
      return
    }
    dots = createDiv("dots")
    slider.track.details.slides.forEach((_e, idx) => {
      var dot = createDiv("dot")
      dot.addEventListener("click", () => slider.moveToIdx(idx))
      dots.appendChild(dot)
    })
    wrapper.appendChild(dots)
  }

  function updateClasses() {
    var slide = slider.track.details.rel
    slide === 0
      ? arrowLeft.classList.add("arrow--disabled")
      : arrowLeft.classList.remove("arrow--disabled")
    slide === slider.track.details.slides.length - 1
      ? arrowRight.classList.add("arrow--disabled")
      : arrowRight.classList.remove("arrow--disabled")
    Array.from(dots.children).forEach(function (dot, idx) {
      idx === slide
        ? dot.classList.add("dot--active")
        : dot.classList.remove("dot--active")
    })
  }

  slider.on("created", () => {
    markup()
    updateClasses()
  })
  slider.on("optionsChanged", () => {
    markup(true)
    markup()
    updateClasses()
  })
  slider.on("slideChanged", () => {
    updateClasses()
  })
  slider.on("destroyed", () => {
    markup(true)
  })
}

function draw_overlay_multipick_coordinates(){
    const picking_diameters = document.getElementById("picking_diameters")
    if(picking_diameters ==  null) return
    const diameter_input = document.getElementsByName("diameter")[0]
    if(diameter_input ==  null) return
    const available_diameters = JSON.parse(picking_diameters.dataset.diameters.replaceAll("'", '"'))
    const selected_diameter = available_diameters[picking_diameters.value - 1]
    diameter_input.value = selected_diameter
    for(const box_overlay of document.getElementsByClassName("boxes_overlay_multipick")){
      const scale = box_overlay.width  / Number(box_overlay.dataset.xdim)
      // account for rectangular images
      var yoffset = (box_overlay.height - (scale * Number(box_overlay.dataset.ydim))) / 2
      const boxes = JSON.parse(box_overlay.dataset.boxes.replaceAll("'", '"'))
      const ctx = box_overlay.getContext("2d");
      ctx.clearRect(0, 0, box_overlay.width, box_overlay.height)
      for(const box of boxes){
        if(box["diameter"] == selected_diameter){
          ctx.strokeStyle = "red";
          ctx.beginPath();
          ctx.arc(box["x"] * scale, (box["y"] * scale) + yoffset, 2, 0, 2 * Math.PI);
          ctx.stroke();
        }
      }
    }
}

function draw_boxes_overlay_coordinates(){
    for(const box_overlay of document.getElementsByClassName("boxes_overlay")){
      const scale = box_overlay.width  / Number(box_overlay.dataset.xdim)
      // account for rectangular images
      var yoffset = (box_overlay.height - (scale * Number(box_overlay.dataset.ydim))) / 2
      const boxes = JSON.parse(box_overlay.dataset.boxes.replaceAll("'", '"'))
      const ctx = box_overlay.getContext("2d");
      ctx.clearRect(0, 0, box_overlay.width, box_overlay.height)
      for(const box of boxes){
        ctx.strokeStyle = "red";
        ctx.beginPath();
        ctx.arc(box["x"] * scale, (box["y"] * scale) + yoffset, 2, 0, 2 * Math.PI);
        ctx.stroke();
      }
    }
}

window.addEventListener("load", () =>{
    var slider = new KeenSlider("#latest_picked_micrographs_slider", {}, [navigation])

},false);

/* draw boxes on load */
window.addEventListener("load", () =>{
    draw_overlay_multipick_coordinates()
    draw_boxes_overlay_coordinates()
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
                        size: 10
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
    if(document.body.dataset.hasOwnProperty("wide") && document.body.dataset.wide == "True"){
      window.parent.document.getElementById('inipick_iframe').width = '800px';
    }else{
      // put back to default]
      window.parent.document.getElementById('inipick_iframe').width = '400px';
    }

},false);

setTimeout(function () {
  // auto reload if not in user_input mode
  if(document.body.dataset.hasOwnProperty("wide") && document.body.dataset.wide != "True"){
    location.reload();
  }
}, 15000);