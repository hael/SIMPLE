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

function toggleTemplate(templ){
    templ.classList.toggle("disabledbutton")
}

function selectRefs(form){
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

window.addEventListener("load", () =>{
    var slider = new KeenSlider("#picking_references_slider", {}, [navigation])
},false);

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
                        size: 10
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
      window.parent.document.getElementById('genpickrefs_iframe').width = '800px';
    }else{
      // put back to default]
      window.parent.document.getElementById('genpickrefs_iframe').width = '400px';
    }

},false);

const box_sizes = [32, 36, 40, 48, 52, 56, 64, 66, 70, 72, 80, 84, 88, 100, 104, 108, 112, 120, 128, 130, 132,
    140, 144, 150, 160, 162, 168, 176, 180, 182, 192, 200, 208, 216, 220, 224, 240, 256, 264, 288, 300, 308, 320, 324, 336,
    338, 352, 364, 384, 400, 420, 432, 448, 450, 462, 480, 486, 500, 504, 512, 520, 528, 546, 560, 576, 588, 600, 640, 648,
    650, 660, 672, 686, 700, 702, 704, 720, 726, 728, 750, 768, 770, 784, 800, 810, 840, 882, 896, 910, 924, 936, 972, 980,
    1008, 1014, 1020, 1024,1296, 1536, 1728, 1944,2048, 2304, 2592, 3072, 3200, 3456, 3888, 4096, 4608, 5000, 5184, 6144,
    6250, 6400, 6912, 7776, 8192, 9216, 10240, 12288, 12500]

function updateBoxSize(){
    const box_size_selector        = document.getElementById("box_size_selector")
    const current_box_size         = document.getElementById("current_box_size")
    const final_selection_boxsize = document.getElementById("final_selection_boxsize")
    current_box_size.innerHTML = box_sizes[box_size_selector.value] + "px"
    final_selection_boxsize.value = box_sizes[box_size_selector.value]
    const scale = Number(box_size_selector.dataset.boxsize) / Number(box_sizes[box_size_selector.value])
    for(const picktemplate of document.getElementsByClassName("picktemplate")){
        picktemplate.style.transform = "scale(" + scale + ")"
    }
}

window.addEventListener("load", () =>{
    const box_size_selector = document.getElementById("box_size_selector")
    const current_box_size  = document.getElementById("current_box_size")
    if(box_size_selector == undefined || current_box_size == undefined) return
    const box_size = Number(box_size_selector.dataset.boxsize)
    box_size_selector.max = box_sizes.length - 1
    box_size_selector.step = 1
    box_size_selector.value = box_sizes.indexOf(box_size)
    updateBoxSize()
},false);

setTimeout(function () {
  // autp reload if not in user_input mode
  if(document.body.dataset.hasOwnProperty("wide") && document.body.dataset.wide != "True"){
    location.reload();
  }
}, 15000);