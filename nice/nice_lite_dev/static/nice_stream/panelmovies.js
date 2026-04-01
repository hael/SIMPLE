let lastinteraction = Date.now();

window.addEventListener("load", () =>{
    for(const astig_histogram of document.getElementsByClassName("astig_histogram")){
        const ctx = astig_histogram.getContext("2d");
        const labels = JSON.parse(astig_histogram.dataset.labels.replaceAll("'", '"'))
        while(labels.length < 24){
          labels.push("")
        }
        const data   = JSON.parse(astig_histogram.dataset.values.replaceAll("'", '"'))
        const data2  = JSON.parse(astig_histogram.dataset.values2.replaceAll("'", '"'))
        new Chart(ctx, {
            type: 'bar',
            options:{
                maintainAspectRatio : false,
                maxBarThickness : 4,
                scales: {
                  x: {
                      display: false,
                      stacked: true
                  },
                  y: {
                      display: true,
                      stacked: true,
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
                }
            },
            data: {
              labels: labels,
              datasets: [{
                  data: data,
                  hoverOffset: 4,
                  backgroundColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                  borderColor:     window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble'),
                },{
                  data: data2,
                  hoverOffset: 4,
                  backgroundColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4bg'),
                  borderColor:     window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble'),
                }
              ]
            }
        })
    }
},false);

window.addEventListener("load", () =>{
    for(const ctfres_histogram of document.getElementsByClassName("ctfres_histogram")){
        const ctx = ctfres_histogram.getContext("2d");
        const labels = JSON.parse(ctfres_histogram.dataset.labels.replaceAll("'", '"'))
        while(labels.length < 24){
          labels.push("")
        }
        const data   = JSON.parse(ctfres_histogram.dataset.values.replaceAll("'", '"'))
        const data2  = JSON.parse(ctfres_histogram.dataset.values2.replaceAll("'", '"'))
        new Chart(ctx, {
            type: 'bar',
            options:{
                maintainAspectRatio : false,
                maxBarThickness : 4,
                scales: {
                  x: {
                      display: false,
                      stacked: true
                  },
                  y: {
                      display: true,
                      stacked: true,
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
                }
            },
            data: {
              labels: labels,
              datasets: [{
                  data: data,
                  hoverOffset: 4,
                  backgroundColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                  borderColor:     window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble'),
                },{
                  data: data2,
                  hoverOffset: 4,
                  backgroundColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4bg'),
                  borderColor:     window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble'),
                }
              ]
            }
        })
    }
},false);

window.addEventListener("load", () =>{
    for(const df_histogram of document.getElementsByClassName("df_histogram")){
        const ctx = df_histogram.getContext("2d");
        const labels = JSON.parse(df_histogram.dataset.labels.replaceAll("'", '"'))
        while(labels.length < 24){
          labels.push("")
        }
        const data   = JSON.parse(df_histogram.dataset.values.replaceAll("'", '"'))
        const data2  = JSON.parse(df_histogram.dataset.values2.replaceAll("'", '"'))
        new Chart(ctx, {
            type: 'bar',
            options:{
                maintainAspectRatio : false,
                maxBarThickness : 4,
                scales: {
                  x: {
                      display: false,
                      stacked: true
                  },
                  y: {
                      display: true,
                      stacked: true,
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
                }
            },
            data: {
              labels: labels,
              datasets: [{
                  data: data,
                  hoverOffset: 4,
                  backgroundColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                  borderColor:     window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble'),
                },{
                  data: data2,
                  hoverOffset: 4,
                  backgroundColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4bg'),
                  borderColor:     window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble'),
                }
              ]
            }
        })
    }
},false);

window.addEventListener("load", () =>{
    for(const rate_histogram of document.getElementsByClassName("rate_histogram")){
        const ctx = rate_histogram.getContext("2d");
        const labels = JSON.parse(rate_histogram.dataset.labels.replaceAll("'", '"'))
        while(labels.length < 24){
          labels.push("")
        }
        const data   = JSON.parse(rate_histogram.dataset.values.replaceAll("'", '"'))
  
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
                      ticks: {
                        font: {
                            size: 8,
                        },
                        maxTicksLimit: 3
                      }
                  },
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
                  hoverOffset: 4,
                  backgroundColor: window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                  borderColor:     window.getComputedStyle(document.body).getPropertyValue('--color-nice4bubble'),
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