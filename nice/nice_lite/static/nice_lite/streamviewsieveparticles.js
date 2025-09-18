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

function toggleSievecls(element){
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  const rejected_cls2D_slider = document.getElementById("rejected_cls2D_slider")
  if(element.parentElement.id == "accepted_cls2D_slider"){
    rejected_cls2D_slider.appendChild(element)
  }else{
    accepted_cls2D_slider.appendChild(element)
  }
}

function selectSievecls(form){
    const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
    const rejected_cls2D_slider = document.getElementById("rejected_cls2D_slider")
    const accepted = []
    const rejected = []
    for(const cls of accepted_cls2D_slider.getElementsByClassName("sievecls")){
        const idx  = Number(cls.dataset.idx)
        accepted.push(idx)
    }
    for(const cls of rejected_cls2D_slider.getElementsByClassName("sievecls")){
        const idx  = Number(cls.dataset.idx)
        rejected.push(idx)
    }
    document.getElementById("accepted_cls2D").value = accepted
    document.getElementById("rejected_cls2D").value = rejected
}

window.addEventListener("load", () =>{
    var latest_accepted_slider = new KeenSlider("#latest_accepted_slider", {}, [navigation])
    var latest_rejected_slider = new KeenSlider("#latest_rejected_slider", {}, [navigation])
},false);

window.addEventListener("load", () =>{
    for(const particles_pie_chart of document.getElementsByClassName("particles_pie_chart")){
        const ctx = particles_pie_chart.getContext("2d");
        const n_imported = Number(particles_pie_chart.dataset.imported) 
        const n_accepted = Number(particles_pie_chart.dataset.accepted)
        const n_rejected = Number(particles_pie_chart.dataset.rejected)
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
    if(document.body.dataset.hasOwnProperty("wide") && document.body.dataset.wide == "True"){
      window.parent.document.getElementById('sieve_iframe').width = '800px';
    }else{
      // put back to default
      window.parent.document.getElementById('sieve_iframe').width = '400px';
    }

},false);

setTimeout(function () {
  // auto reload if not in user_input mode
  if(document.body.dataset.hasOwnProperty("wide") && document.body.dataset.wide != "True"){
    location.reload();
  }
}, 33000);