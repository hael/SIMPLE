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

window.addEventListener("load", () =>{
    for(const optics_scatter_chart of document.getElementsByClassName("optics_scatter_chart")){
        const ctx = optics_scatter_chart.getContext("2d");
        const assignments = JSON.parse(optics_scatter_chart.dataset.assignments.replaceAll("'", '"'))
        const datasets = []
        for(const group of assignments){
            const dataset = {}
            dataset["label"] = group["id"]
            dataset["data"]  = group["coordinates"]
            datasets.push(dataset)
        }
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
                datasets:datasets 
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

window.addEventListener("load", () => {
  const logtext = document.querySelector(".logtext")
  logtext.scrollTop = logtext.scrollHeight - logtext.clientHeight
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
