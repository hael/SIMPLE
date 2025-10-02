window.addEventListener("load", () =>{
    for(const cls2d_histogram of document.getElementsByClassName("cls2d_histogram")){
        const ctx = cls2d_histogram.getContext("2d");
        const labels = JSON.parse(cls2d_histogram.dataset.labels.replaceAll("'", '"'))
        while(labels.length < 24){
          labels.push("")
        }
        const data = JSON.parse(cls2d_histogram.dataset.values.replaceAll("'", '"'))
  
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