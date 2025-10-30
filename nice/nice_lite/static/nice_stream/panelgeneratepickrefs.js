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
  const picking_references = document.getElementById("picking_references")
  picking_references.parentElement.scrollLeft += 72;
  lastinteraction = Date.now();
}

scrlLeft = (element, event) => {
  event.preventDefault()
  const picking_references = document.getElementById("picking_references")
  picking_references.parentElement.scrollLeft -= 72;
  lastinteraction = Date.now();
}

toggleTemplate = (templ) => {
    templ.classList.toggle("disabledbutton")
    const xmark = templ.querySelector('.xmark')
    if(templ.classList.contains("disabledbutton")){
      xmark.classList.remove("hidden")
    }else{
      xmark.classList.add("hidden")
    }
    lastinteraction = Date.now();
}

selectRefs = (form) => {
    const selected = []
    var path       = ""
    for(const pick_template of document.getElementsByClassName("picktemplate")){
        const idx  = Number(pick_template.dataset.idx)
        if(!pick_template.classList.contains("disabledbutton")){
            path = pick_template.dataset.path
            selected.push(idx)
        }
    }
    document.getElementById("final_selection_source").value = path
    document.getElementById("final_selection").value = selected
}

showMenu = (element, event) => {
  event.preventDefault()
  const selectmenu    = element.parentElement.parentElement.parentElement.querySelector("[name='selectmenu']")
  const selectmenubox = selectmenu.querySelector("[name='selectmenubox']")
  const sortpop     = document.querySelector("#sortpop")
  const sortres     = document.querySelector("#sortres")
  const selectabove = document.querySelector("#selectabove")
  const selectbelow = document.querySelector("#selectbelow")
  selectmenubox.style.top  = event.pageY + "px"
  selectmenubox.style.left = event.pageX + "px"
  selectmenu.style.display = "flex"
  sortpop.onclick = () => {sortPop()}
  sortres.onclick = () => {sortRes()}
  if(selectabove != undefined) selectabove.onclick = () => {selectAbove(element)}
  if(selectbelow != undefined) selectbelow.onclick = () => {selectBelow(element)}
  lastinteraction = Date.now();
}

hideMenu = () => {
  for(const selectmenu of document.querySelectorAll("[name='selectmenu']")){
    const sortpop     = document.querySelector("#sortpop")
    const sortres     = document.querySelector("#sortres")
    const selectabove = document.querySelector("#selectabove")
    const selectbelow = document.querySelector("#selectbelow")
    selectmenu.style.display = "none"
    sortpop.onclick = null
    sortres.onclick = null
    if(selectabove != undefined) selectabove.onclick = null
    if(selectbelow != undefined) selectbelow.onclick = null
  }
  lastinteraction = Date.now();
}

sortPop = () => {
  const picking_references = document.querySelector("#picking_references")
  const picktemplates      = Array.from(picking_references.querySelectorAll(".picktemplatecontainer"))
  picktemplates.sort((a, b) => {
        return Number(b.dataset.pop) - Number(a.dataset.pop);
  });
  picktemplates.forEach((item) => {
        picking_references.appendChild(item);
  });
  hideMenu()
}

sortRes = () => {
  const picking_references = document.querySelector("#picking_references")
  const picktemplates      = Array.from(picking_references.querySelectorAll(".picktemplatecontainer"))
  picktemplates.sort((a, b) => {
        return Number(a.dataset.res) - Number(b.dataset.res);
  });
  picktemplates.forEach((item) => {
        picking_references.appendChild(item);
  });
  hideMenu()
}

selectAbove = (element) => {
  let threshold = true
  const picking_references = document.querySelector("#picking_references")
  for(const picktemplate of picking_references.querySelectorAll(".picktemplate")){
    const xmark = picktemplate.querySelector('.xmark') 
    if(threshold){
      picktemplate.classList.remove("disabledbutton")
    }else{
      picktemplate.classList.add("disabledbutton")
    }
    if(picktemplate.classList.contains("disabledbutton")){
      xmark.classList.remove("hidden")
    }else{
      xmark.classList.add("hidden")
    }
    if(picktemplate == element){
      threshold = false
    }
  }
  hideMenu()
}

selectBelow = (element) => {
  let threshold = false
  const picking_references = document.querySelector("#picking_references")
  for(const picktemplate of picking_references.querySelectorAll(".picktemplate")){
    const xmark = picktemplate.querySelector('.xmark') 
    if(picktemplate == element){
      threshold = true
    }
    if(threshold){
      picktemplate.classList.remove("disabledbutton")
    }else{
      picktemplate.classList.add("disabledbutton")
    }
    if(picktemplate.classList.contains("disabledbutton")){
      xmark.classList.remove("hidden")
    }else{
      xmark.classList.add("hidden")
    }
  }
  hideMenu()
}

const box_sizes = [32, 36, 40, 48, 52, 56, 64, 66, 70, 72, 80, 84, 88, 100, 104, 108, 112, 120, 128, 130, 132,
    140, 144, 150, 160, 162, 168, 176, 180, 182, 192, 200, 208, 216, 220, 224, 240, 256, 264, 288, 300, 308, 320, 324, 336,
    338, 352, 364, 384, 400, 420, 432, 448, 450, 462, 480, 486, 500, 504, 512, 520, 528, 546, 560, 576, 588, 600, 640, 648,
    650, 660, 672, 686, 700, 702, 704, 720, 726, 728, 750, 768, 770, 784, 800, 810, 840, 882, 896, 910, 924, 936, 972, 980,
    1008, 1014, 1020, 1024,1296, 1536, 1728, 1944,2048, 2304, 2592, 3072, 3200, 3456, 3888, 4096, 4608, 5000, 5184, 6144,
    6250, 6400, 6912, 7776, 8192, 9216, 10240, 12288, 12500]

updateBrightness = (element) => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--genpickrefs-brightness', element.value / 100);
}

updateContrast = (element) => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--genpickrefs-contrast', element.value / 100);
}

updateScale = (element) => {
  const scale = element.value
  for(const cls2D of document.querySelectorAll(".cls2D")){
    cls2D.style.width  = scale + "px"
    cls2D.style.height = scale + "px"
  }
  for(const cls2Dimg of document.querySelectorAll(".cls2Dimg")){
    cls2Dimg.style.width  = scale + "px"
    cls2Dimg.style.height = scale + "px"
  }
  for(const picktemplate of document.querySelectorAll(".picktemplate")){
    picktemplate.style.width  = scale + "px"
    picktemplate.style.height = scale + "px"
  }
  for(const picktemplateimg of document.querySelectorAll(".picktemplateimg")){
    picktemplateimg.style.width  = scale + "px"
    picktemplateimg.style.height = scale + "px"
  }
}

window.addEventListener("load", () => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--genpickrefs-contrast',   1.0);
  cssroot.style.setProperty('--genpickrefs-brightness', 1.0);
})

window.addEventListener("load", () =>{
    for(const movies_pie_chart of document.getElementsByClassName("particles_pie_chart")){
        const ctx = movies_pie_chart.getContext("2d");
        const n_imported  = Number(movies_pie_chart.dataset.imported) 
        const n_processed = Number(movies_pie_chart.dataset.processed)
        const n_rejected  = Number(movies_pie_chart.dataset.rejected)
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
                    'accepted',
                    'rejected'
                ],
                datasets: [{
                    data: [n_imported - n_processed - n_rejected, n_processed, n_rejected],
                    backgroundColor: [
                      window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                      window.getComputedStyle(document.body).getPropertyValue('--color-nice4success'),
                      window.getComputedStyle(document.body).getPropertyValue('--color-nice4alert'),
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

/*setInterval(function () {
  if((Date.now() - lastinteraction) > 30000 && document.visibilityState !== "hidden"){
    lastinteraction = Date.now();
    location.reload();
  }
}, 1000);*/