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

const scrlRight = (element, event)  => {
  event.preventDefault()
  const slider = document.getElementById("accepted_cls2D_slider").parentElement;
  const atEnd = slider.scrollLeft + slider.clientWidth >= slider.scrollWidth - 1;
  slider.scrollLeft = atEnd ? 0 : slider.scrollLeft + 72;
  lastinteraction = Date.now();
}

const scrlLeft = (element, event)  => {
  event.preventDefault()
  const slider = document.getElementById("accepted_cls2D_slider").parentElement;
  const atStart = slider.scrollLeft <= 0;
  slider.scrollLeft = atStart ? slider.scrollWidth - slider.clientWidth : slider.scrollLeft - 72;
  lastinteraction = Date.now();
}

const scrlRejectedRight = (element, event)  => {
  event.preventDefault()
  const slider = document.getElementById("rejected_cls2D_slider").parentElement;
  const atEnd = slider.scrollLeft + slider.clientWidth >= slider.scrollWidth - 1;
  slider.scrollLeft = atEnd ? 0 : slider.scrollLeft + 72;
  lastinteraction = Date.now();
}

const scrlRejectedLeft = (element, event)  => {
  event.preventDefault()
  const slider = document.getElementById("rejected_cls2D_slider").parentElement;
  const atStart = slider.scrollLeft <= 0;
  slider.scrollLeft = atStart ? slider.scrollWidth - slider.clientWidth : slider.scrollLeft - 72;
  lastinteraction = Date.now();
}

const toggleSievecls = (element)  => {
  if(element.classList.contains("sieveclsselected")){
    for(const xmark of element.querySelectorAll(".xmark")){
      xmark.classList.remove("hidden")
    }
  }else{
    for(const xmark of element.querySelectorAll(".xmark")){
      xmark.classList.add("hidden")
    }
  }
  element.classList.toggle("sieveclsselected")
  updateCounts()
  lastinteraction = Date.now();
}

const selectSievecls = (form)  => {
    const accepted = []
    for(const cls of document.getElementsByClassName("sieveclsselected")){
        const idx  = Number(cls.dataset.idx)
        accepted.push(idx)
    }
    document.getElementById("accepted_cls2D").value = accepted
    form.submit()
}

const showMenu = (element, event)  => {
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

const hideMenu = ()  => {
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

const sortPop = ()  => {
  const accepted_cls2D_slider = document.querySelector("#accepted_cls2D_slider")
  const rejected_cls2D_slider = document.querySelector("#rejected_cls2D_slider")
  const acceptedtemplates     = Array.from(accepted_cls2D_slider.querySelectorAll(".sieveclscontainer"))
  const rejectedtemplates     = Array.from(rejected_cls2D_slider.querySelectorAll(".sieveclscontainer"))
  acceptedtemplates.sort((a, b) => {
        return Number(b.dataset.pop) - Number(a.dataset.pop);
  });
  acceptedtemplates.forEach((item) => {
        accepted_cls2D_slider.appendChild(item);
  });
  rejectedtemplates.sort((a, b) => {
        return Number(b.dataset.pop) - Number(a.dataset.pop);
  });
  rejectedtemplates.forEach((item) => {
        rejected_cls2D_slider.appendChild(item);
  });
  hideMenu()
}

const sortRes = ()  => {
  const accepted_cls2D_slider = document.querySelector("#accepted_cls2D_slider")
  const rejected_cls2D_slider = document.querySelector("#rejected_cls2D_slider")
  const acceptedtemplates     = Array.from(accepted_cls2D_slider.querySelectorAll(".sieveclscontainer"))
  const rejectedtemplates     = Array.from(rejected_cls2D_slider.querySelectorAll(".sieveclscontainer"))
  acceptedtemplates.sort((a, b) => {
        return Number(a.dataset.res) - Number(b.dataset.res);
  });
  acceptedtemplates.forEach((item) => {
        accepted_cls2D_slider.appendChild(item);
  });
  rejectedtemplates.sort((a, b) => {
        return Number(a.dataset.res) - Number(b.dataset.res);
  });
  rejectedtemplates.forEach((item) => {
        rejected_cls2D_slider.appendChild(item);
  });
  hideMenu()
}

const selectAbove = (element)  => {
  let threshold = true
  const accepted_cls2D_slider = document.querySelector("#accepted_cls2D_slider")
  const rejected_cls2D_slider = document.querySelector("#rejected_cls2D_slider")
  for(const clscontainer of element.parentElement.parentElement.querySelectorAll(".sieveclscontainer")){
    if(threshold){
      accepted_cls2D_slider.appendChild(clscontainer)
    }else{
      rejected_cls2D_slider.appendChild(clscontainer)
    }
    if(clscontainer == element.parentElement){
      threshold = false
    }
  }
  updateCounts()
  hideMenu()
}

const selectBelow = (element)  => {
  let threshold = false
  const accepted_cls2D_slider = document.querySelector("#accepted_cls2D_slider")
  const rejected_cls2D_slider = document.querySelector("#rejected_cls2D_slider")
  for(const clscontainer of element.parentElement.parentElement.querySelectorAll(".sieveclscontainer")){
    if(clscontainer == element.parentElement){
      threshold = true
    }
    if(threshold){
      accepted_cls2D_slider.appendChild(clscontainer)
    }else{
      rejected_cls2D_slider.appendChild(clscontainer)
    }
  }
  updateCounts()
  hideMenu()
}

const updateCounts = ()  => {
  const clustercount    = document.querySelector("#clustercount")
  const particlecount   = document.querySelector("#particlecount")
  if(clustercount == undefined || particlecount == undefined) return
  const sieveclscontainers = document.querySelectorAll(".sieveclscontainer")
  const sieveclscontainers_accepted = document.querySelectorAll("#accepted_cls2D_slider .sieveclscontainer")
  let ncls       = 0
  let nptcls     = 0
  let nptcls_tot = 0
  for(const sieveclscontainer of sieveclscontainers){
    const pop = Number(sieveclscontainer.dataset.pop)
    nptcls_tot = nptcls_tot + pop
  }
  for(const sieveclscontainer of sieveclscontainers_accepted){
    const pop = Number(sieveclscontainer.dataset.pop)
    nptcls = nptcls + pop
    ncls = ncls + 1
  }
  clustercount.innerHTML = ncls.toLocaleString() + " / " + sieveclscontainers.length.toLocaleString() 
  particlecount.innerHTML = nptcls.toLocaleString() + " / " + nptcls_tot.toLocaleString()
}

const updateBrightness = (element)  => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--sieve-brightness', element.value / 100);
}

const updateContrast = (element)  => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--sieve-contrast', element.value / 100);
}

const updateScale = (element)  => {
  const scale = element.value
  for(const cls2D of document.querySelectorAll(".sievecls")){
    cls2D.style.width  = scale + "px"
    cls2D.style.height = scale + "px"
  }
  for(const cls2Dimg of document.querySelectorAll(".sieveclsimg")){
    cls2Dimg.style.width  = scale + "px"
    cls2Dimg.style.height = scale + "px"
  }
}

window.addEventListener("load", () => {
  var cssroot = document.querySelector(':root');
  cssroot.style.setProperty('--sieve-contrast',   1.0);
  cssroot.style.setProperty('--sieve-brightness', 1.0);
})

window.addEventListener("load", () => {
  const logtext = document.querySelector(".logtext")
  if (logtext) logtext.scrollTop = logtext.scrollHeight - logtext.clientHeight
})

window.addEventListener("load", () =>{
    for(const canvas of document.getElementsByClassName("particles_pie_chart")){
        const n_imported = Number(canvas.dataset.imported)
        const n_accepted = Number(canvas.dataset.accepted)
        const n_rejected = Number(canvas.dataset.rejected)
        buildDonut(canvas, [n_imported - n_accepted - n_rejected, n_accepted, n_rejected], ['unclassified', 'accepted', 'rejected'])
    }
},false);

window.addEventListener("load", () =>{
  const loadinggauze = document.getElementById("loadinggauze")
  if(loadinggauze != undefined){
    loadinggauze.style.opacity = "0";
    setTimeout(function () {
    document.getElementById("loadinggauze").style.display = "none";
    }, 600);
  }
  updateCounts()
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