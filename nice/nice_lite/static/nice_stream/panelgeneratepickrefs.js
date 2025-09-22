let lastinteraction = Date.now();

scrlRight = () => {
  const micrograph_slider = document.getElementById("picking_references")
  micrograph_slider.scrollLeft += 72;
  lastinteraction = Date.now();
}

scrlLeft = () => {
  const micrograph_slider = document.getElementById("picking_references")
  micrograph_slider.scrollLeft -= 72;
  lastinteraction = Date.now();
}

toggleTemplate = (templ) => {
    templ.classList.toggle("disabledbutton")
    lastinteraction = Date.now();
}

selectRefs = (form) => {
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
                        size: 9
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

const box_sizes = [32, 36, 40, 48, 52, 56, 64, 66, 70, 72, 80, 84, 88, 100, 104, 108, 112, 120, 128, 130, 132,
    140, 144, 150, 160, 162, 168, 176, 180, 182, 192, 200, 208, 216, 220, 224, 240, 256, 264, 288, 300, 308, 320, 324, 336,
    338, 352, 364, 384, 400, 420, 432, 448, 450, 462, 480, 486, 500, 504, 512, 520, 528, 546, 560, 576, 588, 600, 640, 648,
    650, 660, 672, 686, 700, 702, 704, 720, 726, 728, 750, 768, 770, 784, 800, 810, 840, 882, 896, 910, 924, 936, 972, 980,
    1008, 1014, 1020, 1024,1296, 1536, 1728, 1944,2048, 2304, 2592, 3072, 3200, 3456, 3888, 4096, 4608, 5000, 5184, 6144,
    6250, 6400, 6912, 7776, 8192, 9216, 10240, 12288, 12500]

updateBoxSize = () => {
    const box_size_selector        = document.getElementById("box_size_selector")
    const current_box_size         = document.getElementById("current_box_size")
    const final_selection_boxsize = document.getElementById("final_selection_boxsize")
    current_box_size.innerHTML = box_sizes[box_size_selector.value] + "px"
    final_selection_boxsize.value = box_sizes[box_size_selector.value]
    const scale = Number(box_size_selector.dataset.boxsize) / Number(box_sizes[box_size_selector.value])
    for(const picktemplate of document.getElementsByClassName("picktemplate")){
        picktemplate.style.transform = "scale(" + scale + ")"
    }
    lastinteraction = Date.now();
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

window.addEventListener("load", () =>{
  document.getElementById("loadinggauze").style.opacity = "0";
  setTimeout(function () {
   document.getElementById("loadinggauze").style.display = "none";
  }, 600);
})

setInterval(() => {
  if((Date.now() - lastinteraction) > 30000){
    document.getElementById("loadinggauze").style.display = "flex";
    document.getElementById("loadinggauze").style.opacity = "1";
    setTimeout(() => {
      location.reload();
    }, 600);
  }
}, 1000);