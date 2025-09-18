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

window.addEventListener("load", () =>{
    var slider = new KeenSlider("#latest_micrographs_slider", {}, [navigation])
},false);

window.addEventListener("load", () =>{
    for(const movies_pie_chart of document.getElementsByClassName("movies_pie_chart")){
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
    for(const ctfres_histogram of document.getElementsByClassName("ctfres_histogram")){
        const ctx = ctfres_histogram.getContext("2d");
        const labels = JSON.parse(ctfres_histogram.dataset.labels.replaceAll("'", '"'))
        const data   = JSON.parse(ctfres_histogram.dataset.values.replaceAll("'", '"'))
  
        let savedval = Number(ctfres_histogram.dataset.savedval)
        let savedidx = labels.indexOf(savedval)
        let dragactive = false

        const dragger = {
		      id : "dragger",
          beforeEvent(chart, args) {
            switch (args.event.type) {
		          case 'mousemove':
                if(dragactive){
                  const xelements = chart.getElementsAtEventForMode(args.event.native, 'x', { intersect: false }, true)
			  	        if(xelements.length > 0){
                    const firstelement = xelements[0]
                    if(chart.config.options.plugins.annotation.annotations.cutoff.value != firstelement.index){
                      chart.config.options.plugins.annotation.annotations.cutoff.value         = firstelement.index;
                      chart.config.options.plugins.annotation.annotations.cutoff.label.content = labels[firstelement.index] + "Å"
                      chart.config.options.plugins.annotation.annotations.cutoffshadow.xMin    = firstelement.index
                      savedval = labels[firstelement.index]
                      chart.update('none')
                    }
                  }
                }
                break;
              case 'mouseout':
			          break;
		          case 'mouseup':
                dragactive = false
                const update_preprocess_ctfres = document.getElementById("update_preprocess_ctfres")
                const ctfres = document.getElementsByName("ctfres")
                ctfres[0].value = savedval
                update_preprocess_ctfres.submit()
                break
              case 'mousedown':
                dragactive = true
                break
              default:
                break
            }
          }
        }
        
        new Chart(ctx, {
            type: 'bar',
            plugins: [dragger],
            options:{
                maintainAspectRatio : false,
                events: ['mousedown', 'mouseup', 'mousemove', 'mouseout'],
                scales: {
                  x: {
                      display: false,
                  },
                  y: {
                      display: false,
                  }
                },
                plugins:{
                    legend:{
                        display: false,
                    },
                    annotation : {
                        annotations : {
                            cutoff:{
                                type: 'line',
                                borderWidth: 0,
                                label: {
                                    display: true,
                                    content: savedval + "Å",
                                    position: 'start',
                                    backgroundColor : '#d3d3d3',
                                    color : "#585858",
                                    font :{
                                        size: 9
                                    },
                                    padding: 3
                                },
                                scaleID: 'x',
                                value: savedidx,
                                z : 10,
                                enter(element){
                                    return true;
                                },
                                leave(element){
                                    return true;
                                }
                            },
                            cutoffshadow:{
                                type: 'box',
                                backgroundColor: 'rgba(211, 211, 211, 0.2)',
                                borderWidth: 0,
                                xMin: savedidx,
                                xMax :labels.length -1,
                            }
                        }
                    }
                }
            },
            data: {
              labels: labels,
              datasets: [{
                  data: data,
                  hoverOffset: 4
              }]
            }
        })
    }
},false);

window.addEventListener("load", () =>{
    for(const astig_histogram of document.getElementsByClassName("astig_histogram")){
        const ctx = astig_histogram.getContext("2d");
        const labels = JSON.parse(astig_histogram.dataset.labels.replaceAll("'", '"'))
        const data   = JSON.parse(astig_histogram.dataset.values.replaceAll("'", '"'))
  
        let savedval = Number(astig_histogram.dataset.savedval)
        let savedidx = labels.indexOf(savedval)
        let dragactive = false

        const dragger = {
		      id : "dragger",
          beforeEvent(chart, args) {
            switch (args.event.type) {
		          case 'mousemove':
                if(dragactive){
                  const xelements = chart.getElementsAtEventForMode(args.event.native, 'x', { intersect: false }, true)
			  	        if(xelements.length > 0){
                    const firstelement = xelements[0]
                    if(chart.config.options.plugins.annotation.annotations.cutoff.value != firstelement.index){
                      chart.config.options.plugins.annotation.annotations.cutoff.value         = firstelement.index;
                      chart.config.options.plugins.annotation.annotations.cutoff.label.content = labels[firstelement.index] + "%"
                      chart.config.options.plugins.annotation.annotations.cutoffshadow.xMin    = firstelement.index
                      savedval = labels[firstelement.index]
                      chart.update('none')
                    }
                  }
                }
                break;
              case 'mouseout':
			          break;
		          case 'mouseup':
                dragactive = false
                const update_preprocess_astig = document.getElementById("update_preprocess_astig")
                const astigmatism = document.getElementsByName("astigmatism")
                astigmatism[0].value = savedval
                update_preprocess_astig.submit()
                break
              case 'mousedown':
                dragactive = true
                break
              default:
                break
            }
          }
        }
        
        new Chart(ctx, {
            type: 'bar',
            plugins: [dragger],
            options:{
                maintainAspectRatio : false,
                events: ['mousedown', 'mouseup', 'mousemove', 'mouseout'],
                scales: {
                  x: {
                      display: false,
                  },
                  y: {
                      display: false,
                  }
                },
                plugins:{
                    legend:{
                        display: false,
                    },
                    annotation : {
                        annotations : {
                            cutoff:{
                                type: 'line',
                                borderWidth: 0,
                                label: {
                                    display: true,
                                    content: savedval + "%",
                                    position: 'start',
                                    backgroundColor : '#d3d3d3',
                                    color : "#585858",
                                    font :{
                                        size: 9
                                    },
                                    padding: 3
                                },
                                scaleID: 'x',
                                value: savedidx,
                                z : 10,
                                enter(element){
                                    return true;
                                },
                                leave(element){
                                    return true;
                                }
                            },
                            cutoffshadow:{
                                type: 'box',
                                backgroundColor: 'rgba(211, 211, 211, 0.2)',
                                borderWidth: 0,
                                xMin: savedidx,
                                xMax :labels.length -1,
                            }
                        }
                    }
                }
            },
            data: {
              labels: labels,
              datasets: [{
                  data: data,
                  hoverOffset: 4
              }]
            }
        })
    }
},false);

window.addEventListener("load", () =>{
    for(const icescore_histogram of document.getElementsByClassName("icescore_histogram")){
        const ctx = icescore_histogram.getContext("2d");
        const labels = JSON.parse(icescore_histogram.dataset.labels.replaceAll("'", '"'))
        const data   = JSON.parse(icescore_histogram.dataset.values.replaceAll("'", '"'))
  
        let savedval = Number(icescore_histogram.dataset.savedval)
        let savedidx = labels.indexOf(savedval)
        let dragactive = false

        const dragger = {
		      id : "dragger",
          beforeEvent(chart, args) {
            switch (args.event.type) {
		          case 'mousemove':
                if(dragactive){
                  const xelements = chart.getElementsAtEventForMode(args.event.native, 'x', { intersect: false }, true)
			  	        if(xelements.length > 0){
                    const firstelement = xelements[0]
                    if(chart.config.options.plugins.annotation.annotations.cutoff.value != firstelement.index){
                      chart.config.options.plugins.annotation.annotations.cutoff.value         = firstelement.index;
                      chart.config.options.plugins.annotation.annotations.cutoff.label.content = Number(labels[firstelement.index]).toFixed(2);
                      chart.config.options.plugins.annotation.annotations.cutoffshadow.xMin    = firstelement.index
                      savedval = labels[firstelement.index]
                      chart.update('none')
                    }
                  }
                }
                break;
              case 'mouseout':
			          break;
		          case 'mouseup':
                dragactive = false
                const update_preprocess_icescore = document.getElementById("update_preprocess_icescore")
                const icescore = document.getElementsByName("icescore")
                icescore[0].value = savedval
                update_preprocess_icescore.submit()
                break
              case 'mousedown':
                dragactive = true
                break
              default:
                break
            }
          }
        }
        
        new Chart(ctx, {
            type: 'bar',
            plugins: [dragger],
            options:{
                maintainAspectRatio : false,
                events: ['mousedown', 'mouseup', 'mousemove', 'mouseout'],
                scales: {
                  x: {
                      display: false,
                  },
                  y: {
                      display: false,
                  }
                },
                plugins:{
                    legend:{
                        display: false,
                    },
                    annotation : {
                        annotations : {
                            cutoff:{
                                type: 'line',
                                borderWidth: 0,
                                label: {
                                    display: true,
                                    content: Number(savedval).toFixed(2),
                                    position: 'start',
                                    backgroundColor : '#d3d3d3',
                                    color : "#585858",
                                    font :{
                                        size: 9
                                    },
                                    padding: 3
                                },
                                scaleID: 'x',
                                value: savedidx,
                                z : 10,
                                enter(element){
                                    return true;
                                },
                                leave(element){
                                    return true;
                                }
                            },
                            cutoffshadow:{
                                type: 'box',
                                backgroundColor: 'rgba(211, 211, 211, 0.2)',
                                borderWidth: 0,
                                xMin: savedidx,
                                xMax :labels.length -1,
                            }
                        }
                    }
                }
            },
            data: {
              labels: labels,
              datasets: [{
                  data: data,
                  hoverOffset: 4
              }]
            }
        })
    }
},false);


setTimeout(function () {
   location.reload();
}, 30000);