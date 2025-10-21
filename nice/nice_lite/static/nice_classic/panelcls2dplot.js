window.addEventListener("load", () =>{
    for(const micrographs_plot of document.getElementsByClassName("cls2d_plot")){
        const ctx = micrographs_plot.getContext("2d");
        const data   = JSON.parse(micrographs_plot.dataset.values.replaceAll("'", '"'))
        new Chart(ctx, {
            type: 'scatter',
            options:{
                scales: {
                    x: {
                        type: 'linear',
                        position: 'bottom',
                        ticks: {
                          stepSize: 1,
                        }
                    },
                    y:{
                        type: 'linear',
                        position: 'bottom',
                        ticks: {
                          stepSize: 1,
                        }
                    }
                },
                plugins: {
                    legend: {
                        display : false
                    },
                },
                maintainAspectRatio : false,
            },
            data: {
                datasets: [{
                  data: data,
                  backgroundColor: [
                    window.getComputedStyle(document.body).getPropertyValue('--color-nice4header')
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