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
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  accepted_cls2D_slider.parentElement.scrollLeft += 77;
  lastinteraction = Date.now();
}

scrlLeft = (element, event) => {
  event.preventDefault()
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  accepted_cls2D_slider.parentElement.scrollLeft -= 77;
  lastinteraction = Date.now();
}

toggleCls = (element) => {
    element.classList.toggle("disabledbutton")
    const mskcanvas = element.querySelector('.mskcanvas')
    const xmark     = element.querySelector('.xmark')
    if(element.classList.contains("disabledbutton")){
      if(mskcanvas != undefined) mskcanvas.classList.add("hidden")
      xmark.classList.remove("hidden")
    }else{
      xmark.classList.add("hidden")
      if(mskcanvas != undefined) mskcanvas.classList.remove("hidden")
    }
    lastinteraction = Date.now();
    updateCounts()
}

selectCls = (element) => {
    const selected = []
    for(const cls2D of document.getElementsByClassName("cls2D")){
        const idx  = Number(cls2D.dataset.idx)
        if(!cls2D.classList.contains("disabledbutton")){
            selected.push(idx)
        }
    }
    document.getElementById("snapshot_selection").value = selected
    const loadinggauze = document.getElementById("loadinggauze")
    loadinggauze.innerHTML = "generating snapshot ..."
    loadinggauze.classList.remove("hidden")
    loadinggauze.style.opacity = "1";
    element.form.submit()
    timeout = setTimeout(() => { 
      // we know submit doesnt return for at least 2 seconds
      const psets_iframe = window.parent.document.querySelector("#psets_iframe")
      psets_iframe.src = psets_iframe.src
    }, 1500);
}

selectClsFinal = (element) => {
    const deselected = []
    for(const cls2D of document.getElementsByClassName("cls2D")){
        const idx  = Number(cls2D.dataset.idx)
        if(cls2D.classList.contains("disabledbutton")){
            deselected.push(idx)
        }
    }
    document.getElementById("final_deselection").value = deselected
    const loadinggauze = document.getElementById("loadinggauze")
    loadinggauze.innerHTML = "generating final selection ..."
    loadinggauze.classList.remove("hidden")
    loadinggauze.style.opacity = "1";
    element.form.submit()
    timeout = setTimeout(() => { 
      // we know submit doesnt return for at least 2 seconds
      const psets_iframe = window.parent.document.querySelector("#psets_iframe")
      psets_iframe.src = psets_iframe.src
    }, 1500);
}

drawMask = () => {
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
    lastinteraction = Date.now();
  }
}

updateMskdiam = (element) => {
  const current_mskdiam     = document.getElementById("current_mskdiam")
  const selected_mskdiam    = document.getElementById("selected_mskdiam")
  const mskdiam             = element.value * 2 // multiply by 2 to ensure even
  current_mskdiam.innerHTML = mskdiam + "Å" 
  selected_mskdiam.value    = mskdiam
  drawMask()
}

updateMskdiamSubmit = (element) => {
  const loadinggauze = document.getElementById("loadinggauze")
  loadinggauze.innerHTML = "updating mask diameter ..."
  loadinggauze.classList.remove("hidden")
  loadinggauze.style.opacity = "1";
  element.form.submit()
}

showMenu = (element, event) => {
  event.preventDefault()
  const selectmenu    = element.parentElement.parentElement.parentElement.querySelector("[name='selectmenu']")
  const selectmenubox = selectmenu.querySelector("[name='selectmenubox']")
  const sortpop       = selectmenu.querySelector("#sortpop")
  const sortres       = selectmenu.querySelector("#sortres")
  const selectabove   = selectmenu.querySelector("#selectabove")
  const selectbelow   = selectmenu.querySelector("#selectbelow")
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
    const sortpop     = selectmenu.querySelector("#sortpop")
    const sortres     = selectmenu.querySelector("#sortres")
    const selectabove = selectmenu.querySelector("#selectabove")
    const selectbelow = selectmenu.querySelector("#selectbelow")
    selectmenu.style.display = "none"
    sortpop.onclick = null
    sortres.onclick = null
    if(selectabove != undefined) selectabove.onclick = null
    if(selectbelow != undefined) selectbelow.onclick = null
  }
  lastinteraction = Date.now();
}

sortPop = () => {
  const accepted_cls2D_slider = document.querySelector("#accepted_cls2D_slider")
  const acceptedtemplates     = Array.from(accepted_cls2D_slider.querySelectorAll(".cls2dcontainer"))
  acceptedtemplates.sort((a, b) => {
        return Number(b.dataset.pop) - Number(a.dataset.pop);
  });
  acceptedtemplates.forEach((item) => {
        accepted_cls2D_slider.appendChild(item);
  });
  hideMenu()
}

sortRes = () => {
  const accepted_cls2D_slider = document.querySelector("#accepted_cls2D_slider")
  const acceptedtemplates     = Array.from(accepted_cls2D_slider.querySelectorAll(".cls2dcontainer"))
  acceptedtemplates.sort((a, b) => {
        return Number(a.dataset.res) - Number(b.dataset.res);
  });
  acceptedtemplates.forEach((item) => {
        accepted_cls2D_slider.appendChild(item);
  });
  hideMenu()
}

selectAbove = (element) => {
  let threshold = true
  for(const clscontainer of element.parentElement.parentElement.querySelectorAll(".cls2D")){
    const mskcanvas = clscontainer.querySelector('.mskcanvas')
    const xmark     = clscontainer.querySelector('.xmark') 
    if(threshold){
      clscontainer.classList.remove("disabledbutton")
    }else{
      clscontainer.classList.add("disabledbutton")
    }
    if(clscontainer.classList.contains("disabledbutton")){
      if(mskcanvas != undefined) mskcanvas.classList.add("hidden")
      xmark.classList.remove("hidden")
    }else{
      xmark.classList.add("hidden")
      if(mskcanvas != undefined) mskcanvas.classList.remove("hidden")
    }
    if(clscontainer == element){
      threshold = false
    }
  }
  updateCounts()
  hideMenu()
}

selectBelow = (element) => {
  let threshold = false
  for(const clscontainer of element.parentElement.parentElement.querySelectorAll(".cls2D")){
    const mskcanvas = clscontainer.querySelector('.mskcanvas')
    const xmark     = clscontainer.querySelector('.xmark') 
    if(clscontainer == element){
      threshold = true
    }
    if(threshold){
      clscontainer.classList.remove("disabledbutton")
    }else{
      clscontainer.classList.add("disabledbutton")
    }
    if(clscontainer.classList.contains("disabledbutton")){
      if(mskcanvas != undefined) mskcanvas.classList.add("hidden")
      xmark.classList.remove("hidden")
    }else{
      xmark.classList.add("hidden")
      if(mskcanvas != undefined) mskcanvas.classList.remove("hidden")
    }
  }
  updateCounts()
  hideMenu()
}

updateCounts = () => {
  const clustercount          = document.querySelector("#clustercount")
  const particlecount         = document.querySelector("#particlecount")
  const final_selection_ptcls = document.querySelector("#final_selection_ptcls")
  const cls2dcontainers       = document.querySelectorAll(".cls2dcontainer")
  ncls       = 0
  nptcls     = 0
  nptcls_tot = 0
  for(const cls2dcontainer of cls2dcontainers){
    const deselected = cls2dcontainer.querySelector(".disabledbutton")
    const pop = Number(cls2dcontainer.dataset.pop)
    const idx = Number(cls2dcontainer.dataset.idx)
    nptcls_tot = nptcls_tot + pop
    if(deselected == undefined){
      nptcls = nptcls + pop
      ncls   = ncls + 1
    }
  }
  clustercount.innerHTML = ncls.toLocaleString() + " / " + cls2dcontainers.length.toLocaleString() 
  particlecount.innerHTML = nptcls.toLocaleString() + " / " + nptcls_tot.toLocaleString()
  if(final_selection_ptcls != undefined) final_selection_ptcls.value = nptcls
}

updateBrightness = (element) => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--classification2D-brightness', element.value / 100);
}

updateContrast = (element) => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--classification2D-contrast', element.value / 100);
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
  for(const mskcanvas of document.querySelectorAll(".mskcanvas")){
    mskcanvas.style.width  = scale + "px"
    mskcanvas.style.height = scale + "px"
  }
}

window.addEventListener("load", () => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--classification2D-contrast',   1.0);
  cssroot.style.setProperty('--classification2D-brightness', 1.0);
})

window.addEventListener("load", () => {
  const logtext = document.querySelector(".logtext")
  logtext.scrollTop = logtext.scrollHeight - logtext.clientHeight
})

window.addEventListener("load", () =>{
    for(const movies_pie_chart of document.getElementsByClassName("particles_pie_chart")){
        const ctx = movies_pie_chart.getContext("2d");
        const n_imported = Number(movies_pie_chart.dataset.imported) 
        const n_accepted = Number(movies_pie_chart.dataset.accepted)
        const n_rejected = Number(movies_pie_chart.dataset.rejected)
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
         //       'queued',
                'accepted',
                'rejected'
            ],
            datasets: [{
             //   data: [n_imported - n_accepted - n_rejected, n_accepted, n_rejected],
                data: [n_accepted, n_rejected],
                backgroundColor: [
               //     window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
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
  const current_mskdiam  = document.getElementById("current_mskdiam")
  const selected_mskdiam = document.getElementById("selected_mskdiam")
  const mskdiam_selector = document.getElementById("mskdiam_selector")
  if(current_mskdiam == undefined || selected_mskdiam == undefined || mskdiam_selector == undefined) return
  for(const cls2D of document.getElementsByClassName("cls2D")){
    const canvas   = cls2D.getElementsByClassName("mskcanvas")[0]
    const mskdiam  = Number(canvas.dataset.mskdiam)
    current_mskdiam.innerHTML = mskdiam + "Å" 
    selected_mskdiam.value    = mskdiam
    const half_mskdiam = Math.round(mskdiam/2)
    mskdiam_selector.value = half_mskdiam
    break
  }
  updateCounts()
  drawMask()
},false);

window.addEventListener("load", () =>{
  document.getElementById("loadinggauze").style.opacity = "0";
  setTimeout(function () {
   document.getElementById("loadinggauze").classList.add("hidden");
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