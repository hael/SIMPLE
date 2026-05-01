let lastinteraction = Date.now();

// Build a doughnut chart with stream* CSS colour variables.
const buildDonut = (canvas, dataValues, labels) => {
    const style = getComputedStyle(document.body);
    new Chart(canvas.getContext('2d'), {
        type: 'doughnut',
        options: {
            responsive: false,
            plugins: { legend: { display: false } }
        },
        data: {
            labels,
            datasets: [{
                data: dataValues,
                backgroundColor: [
                    style.getPropertyValue('--color-streamring').trim(),
                    style.getPropertyValue('--color-streamicon').trim(),
                    style.getPropertyValue('--color-streamrejected').trim()
                ],
                hoverOffset: 4
            }]
        }
    });
};

const restartProcess = (element)  => {
  const confirmed = confirm("Please confirm that you wish to restart this process");
  if(confirmed){
    element.form.submit()
  }
}

const stopProcess = (element)  => {
  const confirmed = confirm("Please confirm that you wish to stop this process");
  if(confirmed){
    element.form.submit()
  }
}

<<<<<<< Updated upstream
<<<<<<< HEAD
const scrlMicRight = ()  => {
  const micrograph_slider = document.getElementById("micrographs_slider")
  const rect = micrograph_slider.getBoundingClientRect();
  micrograph_slider.scrollLeft += rect.width;
  lastinteraction = Date.now();
}

const scrlMicLeft = ()  => {
  const micrograph_slider = document.getElementById("micrographs_slider")
  const rect = micrograph_slider.getBoundingClientRect();
  micrograph_slider.scrollLeft -= rect.width;
=======
const scrlMicRight = (element, event) => {
  event.preventDefault();
  const slider = element.previousElementSibling;
  const atEnd = slider.scrollLeft + slider.clientWidth >= slider.scrollWidth - 1;
  slider.scrollLeft = atEnd ? 0 : slider.scrollLeft + 235;
  lastinteraction = Date.now();
}

=======
const scrlMicRight = (element, event) => {
  event.preventDefault();
  const slider = element.previousElementSibling;
  const atEnd = slider.scrollLeft + slider.clientWidth >= slider.scrollWidth - 1;
  slider.scrollLeft = atEnd ? 0 : slider.scrollLeft + 235;
  lastinteraction = Date.now();
}

>>>>>>> Stashed changes
const scrlMicLeft = (element, event) => {
  event.preventDefault();
  const slider = element.nextElementSibling;
  const atStart = slider.scrollLeft <= 0;
  slider.scrollLeft = atStart ? slider.scrollWidth - slider.clientWidth : slider.scrollLeft - 235;
<<<<<<< Updated upstream
>>>>>>> a1e410fad146030f3fcbc61f288170a806ef2b04
=======
>>>>>>> Stashed changes
  lastinteraction = Date.now();
}

var pick_observer = new IntersectionObserver(
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

const draw_overlay_coordinates = ()  => {
    const boxes_overlay = document.getElementById("boxes_overlay")
    if(boxes_overlay ==  null) return
    const xdim = Number(boxes_overlay.dataset.xdim)
    const ydim = Number(boxes_overlay.dataset.ydim)
    const scale = xdim > ydim
        ? boxes_overlay.width  / xdim
        : boxes_overlay.height / ydim
    // account for rectangular images
    const xoffset = (boxes_overlay.width  - (scale * xdim)) / 2
    const yoffset = (boxes_overlay.height - (scale * ydim)) / 2
    const boxes = JSON.parse(boxes_overlay.dataset.boxes.replaceAll("'", '"'))
    const ctx = boxes_overlay.getContext("2d");
    ctx.clearRect(0, 0, boxes_overlay.width, boxes_overlay.height)
    for(const box of boxes){
      ctx.strokeStyle = getComputedStyle(document.body).getPropertyValue('--color-streamaction').trim();
      ctx.beginPath();
      ctx.arc((box["x"] * scale) + xoffset, (box["y"] * scale) + yoffset, 1, 0, 2 * Math.PI);
      ctx.stroke();
    }
}

const showMenu = (element, event)  => {
  event.preventDefault()
  const selectmenu    = element.parentElement.parentElement.querySelector("[name='selectmenu']")
  const selectmenubox = selectmenu.querySelector("[name='selectmenubox']")
  selectmenubox.style.top  = event.pageY + "px"
  selectmenubox.style.left = event.pageX + "px"
  selectmenu.style.display = "flex"
  lastinteraction = Date.now();
}

const hideMenu = ()  => {
  for(const selectmenu of document.querySelectorAll("[name='selectmenu']")){
    selectmenu.style.display = "none"
  }
  lastinteraction = Date.now();
}

const updateBrightness = (element)  => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--refpick-brightness', element.value / 100);
}

const updateContrast = (element)  => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--refpick-contrast', element.value / 100);
}

const updateScale = (element)  => {
  const micrographs_slider = document.querySelector("#micrographs_slider")
  const boxes_overlay      = document.querySelector("#boxes_overlay")
  if(micrographs_slider != undefined) micrographs_slider.style.width = element.value + "px"
  if(boxes_overlay      != undefined) {
    boxes_overlay.height = element.value
    boxes_overlay.width  = element.value
  }
  for(const miccontainer of document.querySelectorAll(".miccontainer")){
    const img = miccontainer.querySelector("img")
    if(img != undefined){
      img.style.height = element.value + "px"
      img.style.width  = element.value + "px"
    }
    const div = miccontainer.querySelector("div")
    if(div != undefined){
      div.style.height = element.value + "px"
      div.style.width  = element.value + "px"
    }
  }
  draw_overlay_coordinates()
}

window.addEventListener("load", () => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--refpick-contrast',   1.0);
  cssroot.style.setProperty('--refpick-brightness', 0.5);
})

window.addEventListener("load", () => {
  const logtext = document.querySelector(".logtext")
  if (logtext) logtext.scrollTop = logtext.scrollHeight - logtext.clientHeight
})

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
        ctx.strokeStyle = getComputedStyle(document.body).getPropertyValue('--color-streamaction').trim();
        ctx.beginPath();
        ctx.arc((box["x"] * scale) + xoffset, (box["y"] * scale) + yoffset, 1, 0, 2 * Math.PI);
        ctx.stroke();
      }
    }
},false);

window.addEventListener("load", () =>{
    for(const canvas of document.getElementsByClassName("micrographs_pie_chart")){
        const n_imported  = Number(canvas.dataset.imported)
        const n_processed = Number(canvas.dataset.processed)
        const n_rejected  = Number(canvas.dataset.rejected)
        buildDonut(canvas, [n_imported - n_processed - n_rejected, n_processed, n_rejected], ['queued', 'accepted', 'rejected'])
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
  if((Date.now() - lastinteraction) > 10_000 && document.visibilityState !== "hidden"){
    lastinteraction = Date.now();
    location.reload();
  }
}, 1000);