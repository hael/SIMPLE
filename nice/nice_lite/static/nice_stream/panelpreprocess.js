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
  const micrograph_slider = document.getElementById("micrograph_slider")
  micrograph_slider.scrollLeft += 200;
  lastinteraction = Date.now();
}

scrlLeft= (element, event) => {
  event.preventDefault()
  const micrograph_slider = document.getElementById("micrograph_slider")
  micrograph_slider.scrollLeft -= 200;
  lastinteraction = Date.now();
}

showMenu = (element, event) => {
  event.preventDefault()
  const selectmenu    = element.parentElement.parentElement.querySelector("[name='selectmenu']")
  const selectmenubox = selectmenu.querySelector("[name='selectmenubox']")
  selectmenu.style.display = "flex"
  const selectmenurect    = selectmenu.getBoundingClientRect()
  const selectmenuboxrect = selectmenubox.getBoundingClientRect()
  if(event.pageX + selectmenuboxrect.width > selectmenurect.width){
    selectmenubox.style.left = selectmenu.width - selectmenubox.width + "px"
  } else{
    selectmenubox.style.left = event.pageX + "px"
  }
  if(event.pageY + selectmenuboxrect.height > selectmenurect.height){
    selectmenubox.style.top = selectmenu.height - selectmenubox.height + "px"
  } else{
    selectmenubox.style.top  = event.pageY + "px"
  }
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
  cssroot.style.setProperty('--preprocess-brightness', element.value / 100);
}

updateContrast = (element) => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--preprocess-contrast', element.value / 100);
}

updateScale = (element) => {
  const width = element.value
  for(const miccontainer of document.querySelectorAll(".miccontainer")){
    miccontainer.style.width = width + "px"
  }
}

window.addEventListener("load", () => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--preprocess-contrast',   1.0);
  cssroot.style.setProperty('--preprocess-brightness', 1.0);
})

window.addEventListener("load", () => {
  const logtext = document.querySelector(".logtext")
  logtext.scrollTop = logtext.scrollHeight - logtext.clientHeight
})

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
                lastinteraction = Date.now();
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
                      display: true,
                      ticks: {
                        font: {
                            size: 8,
                        },
                        maxTicksLimit: 3
                      }
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
                                        size: 8
                                    },
                                    padding: 2
                                },
                                scaleID: 'x',
                                value: savedidx,
                                z : 10,
                                enter(element){
                                    ctfres_histogram.style.cursor = 'grab'; 
                                    return true;
                                },
                                leave(element){
                                    ctfres_histogram.style.cursor = 'default'; 
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
                  data:            data,
                  backgroundColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                  borderColor:     window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble'),
                  hoverOffset: 4
              }], 
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
                lastinteraction = Date.now();
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
                      display: true,
                      ticks: {
                        font: {
                            size: 8,
                        },
                        maxTicksLimit: 3
                      }
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
                                        size: 8
                                    },
                                    padding: 2
                                },
                                scaleID: 'x',
                                value: savedidx,
                                z : 10,
                                enter(element){
                                    astig_histogram.style.cursor = 'grab'; 
                                    return true;
                                },
                                leave(element){
                                    astig_histogram.style.cursor = 'default'; 
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
                  data:            data,
                  backgroundColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                  borderColor:     window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble'),
                  hoverOffset:     4
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
                lastinteraction = Date.now();
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
                      display: true,
                      ticks: {
                        font: {
                            size: 8,
                        },
                        maxTicksLimit: 3
                      }
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
                                        size: 8
                                    },
                                    padding: 2
                                },
                                scaleID: 'x',
                                value: savedidx,
                                z : 10,
                                enter(element){
                                    icescore_histogram.style.cursor = 'grab'; 
                                    return true;
                                },
                                leave(element){
                                    icescore_histogram.style.cursor = 'default'; 
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
                  data:            data,
                  backgroundColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                  borderColor:     window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble'),
                  hoverOffset:     4
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
