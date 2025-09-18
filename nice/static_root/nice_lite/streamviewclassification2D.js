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

function toggleCls(element){
    element.classList.toggle("disabledbutton")
}

function selectCls(form){
    const selected = []
    var path       = ""
    for(const cls2D of document.getElementsByClassName("cls2D")){
        const idx  = Number(cls2D.dataset.idx)
        if(!cls2D.classList.contains("disabledbutton")){
            selected.push(idx)
        }
    }
    document.getElementById("snapshot_selection").value = selected
}

function drawMask(){
  const selected_mskdiam = document.getElementById("selected_mskdiam")
  if(selected_mskdiam.value == "") return
  for(const cls2D of document.getElementsByClassName("cls2D")){
    const canvas   = cls2D.getElementsByClassName("mskcanvas")[0]
    const mskscale = Number(canvas.dataset.mskscale)
    const ctx = canvas.getContext("2d")
    ctx.strokeStyle = "yellow";
    ctx.clearRect(0, 0, canvas.width, canvas.height)
    ctx.beginPath();
    ctx.arc(canvas.width / 2, canvas.height / 2, Number(selected_mskdiam.value) * canvas.width / (mskscale * 2), 0, 2 * Math.PI);
    ctx.stroke();
  }
}

function updateMskdiam(element){
  const current_mskdiam     = document.getElementById("current_mskdiam")
  const selected_mskdiam    = document.getElementById("selected_mskdiam")
  const mskdiam             = element.value * 2 // multiply by 2 to ensure even
  current_mskdiam.innerHTML = mskdiam + "Å" 
  selected_mskdiam.value    = mskdiam
  drawMask()
}

window.addEventListener("load", () =>{
    var latest_accepted_slider = new KeenSlider("#latest_accepted_slider", {}, [navigation])
    var latest_rejected_slider = new KeenSlider("#latest_rejected_slider", {}, [navigation])
},false);

window.addEventListener("load", () =>{
    for(const movies_pie_chart of document.getElementsByClassName("particles_pie_chart")){
        const ctx = movies_pie_chart.getContext("2d");
        const n_imported = Number(movies_pie_chart.dataset.imported) 
        const n_accepted = Number(movies_pie_chart.dataset.accepted)
        const n_rejected = Number(movies_pie_chart.dataset.rejected)
        new Chart(ctx, {
            type: 'doughnut',
            data: {
            labels: [
                'particles processing',
                'particles accepted',
                'particles rejected'
            ],
            datasets: [{
                data: [n_imported - n_accepted - n_rejected, n_accepted, n_rejected],
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
  const current_mskdiam  = document.getElementById("current_mskdiam")
  const selected_mskdiam = document.getElementById("selected_mskdiam")
  const mskdiam_selector = document.getElementById("mskdiam_selector")
  for(const cls2D of document.getElementsByClassName("cls2D")){
    const canvas   = cls2D.getElementsByClassName("mskcanvas")[0]
    const mskdiam  = Number(canvas.dataset.mskdiam)
    current_mskdiam.innerHTML = mskdiam + "Å" 
    selected_mskdiam.value    = mskdiam
    const half_mskdiam = Math.round(mskdiam/2)
    mskdiam_selector.value = half_mskdiam
    break
  }
  drawMask()
},false);

window.addEventListener("load", () =>{
    if(document.body.dataset.hasOwnProperty("wide") && document.body.dataset.wide == "True"){
      window.parent.document.getElementById('cls2D_iframe').width = '800px';
    }else{
      // put back to default
      window.parent.document.getElementById('cls2D_iframe').width = '400px';
    }
},false);

setTimeout(function () {
  // auto reload if not in user_input mode
  if(document.body.dataset.hasOwnProperty("wide") && document.body.dataset.wide != "True"){
    location.reload();
  }
}, 33000);