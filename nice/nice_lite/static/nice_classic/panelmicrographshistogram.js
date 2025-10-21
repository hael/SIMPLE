window.addEventListener("load", () =>{
    for(const micrographs_histogram of document.getElementsByClassName("micrographs_histogram")){
        const ctx = micrographs_histogram.getContext("2d");
        const labels = JSON.parse(micrographs_histogram.dataset.labels.replaceAll("'", '"'))
        while(labels.length < 24){
          labels.push("")
        }
        const data   = JSON.parse(micrographs_histogram.dataset.values.replaceAll("'", '"'))
  
        new Chart(ctx, {
            type: 'bar',
            options:{
                maintainAspectRatio : false,
                maxBarThickness : 4,
                scales: {
                  x: {
                      display: false,
                  },
                  y: {
                      display: true,
                  }
                },
                plugins:{
                    legend:{
                        display: false,
                    },
                }
            },
            data: {
              labels: labels,
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