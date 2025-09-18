scrlRight = () => {
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  accepted_cls2D_slider.scrollLeft += 77;
}

scrlLeft = () => {
  const accepted_cls2D_slider = document.getElementById("accepted_cls2D_slider")
  accepted_cls2D_slider.scrollLeft -= 77;
}

toggleCls = (element) => {
    element.classList.toggle("disabledbutton")
}

selectCls = (form) => {
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
                        size: 10
                      }
                    }
                }
              }
            },
            data: {
            labels: [
                'processing',
                'accepted',
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
  document.getElementById("loadinggauze").style.opacity = "0";
  setTimeout(function () {
   document.getElementById("loadinggauze").style.display = "none";
  }, 600);
})

setTimeout(function () {
  document.getElementById("loadinggauze").style.display = "flex";
  document.getElementById("loadinggauze").style.opacity = "1";
  setTimeout(function () {
   location.reload();
  }, 600);
}, 31000);