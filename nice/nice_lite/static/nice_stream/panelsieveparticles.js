let lastinteraction = Date.now();

scrlRight = () => {
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  accepted_cls2D_slider.parentElement.scrollLeft += 72;
  lastinteraction = Date.now();
}

scrlLeft = () => {
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  accepted_cls2D_slider.parentElement.scrollLeft -= 72;
  lastinteraction = Date.now();
}

scrlRejectedRight = () => {
  const rejected_cls2D_slider = document.getElementById("rejected_cls2D_slider")
  rejected_cls2D_slider.parentElement.scrollLeft += 72;
  lastinteraction = Date.now();
}

scrlRejectedLeft = () => {
  const rejected_cls2D_slider = document.getElementById("rejected_cls2D_slider")
  rejected_cls2D_slider.parentElement.scrollLeft -= 72;
  lastinteraction = Date.now();
}

toggleSievecls = (element) => {
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  const rejected_cls2D_slider = document.getElementById("rejected_cls2D_slider")
  if(element.parentElement.parentElement.id == "accepted_cls2D_slider"){
    rejected_cls2D_slider.appendChild(element.parentElement)
  }else{
    accepted_cls2D_slider.appendChild(element.parentElement)
  }
  lastinteraction = Date.now();
}

selectSievecls = (form) => {
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

sortRes = () => {
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

selectAbove = (element) => {
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
  hideMenu()
}

selectBelow = (element) => {
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
  hideMenu()
}

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