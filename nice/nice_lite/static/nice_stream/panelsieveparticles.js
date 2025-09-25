let lastinteraction = Date.now();

scrlRight = () => {
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  accepted_cls2D_slider.scrollLeft += 72;
  lastinteraction = Date.now();
}

scrlLeft = () => {
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  accepted_cls2D_slider.scrollLeft -= 72;
  lastinteraction = Date.now();
}

scrlRejectedRight = () => {
  const rejected_cls2D_slider = document.getElementById("rejected_cls2D_slider")
  rejected_cls2D_slider.scrollLeft += 72;
  lastinteraction = Date.now();
}

scrlRejectedLeft = () => {
  const rejected_cls2D_slider = document.getElementById("rejected_cls2D_slider")
  rejected_cls2D_slider.scrollLeft -= 72;
  lastinteraction = Date.now();
}

toggleSievecls = (element) => {
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  const rejected_cls2D_slider = document.getElementById("rejected_cls2D_slider")
  if(element.parentElement.parentElement.parentElement.id == "accepted_cls2D_slider"){
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