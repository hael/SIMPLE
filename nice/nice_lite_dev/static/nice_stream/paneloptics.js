let lastinteraction = Date.now();

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

window.addEventListener("load", () =>{
    for(const optics_scatter_chart of document.getElementsByClassName("optics_scatter_chart")){
        const ctx = optics_scatter_chart.getContext("2d");
        const assignments = JSON.parse(optics_scatter_chart.dataset.assignments.replaceAll("'", '"'))
        const datasets = []
    const points = []
        for(const group of assignments){
            const dataset = {}
            dataset["label"] = group["id"]
            dataset["data"]  = group["coordinates"]
            datasets.push(dataset)
      for(const coordinate of group["coordinates"]){
        points.push(coordinate)
      }
        }

    const xs = points.map((point) => point.x)
    const ys = points.map((point) => point.y)
    const minX = xs.length > 0 ? Math.min(...xs) : 0
    const maxX = xs.length > 0 ? Math.max(...xs) : 1
    const minY = ys.length > 0 ? Math.min(...ys) : 0
    const maxY = ys.length > 0 ? Math.max(...ys) : 1
    const globalMin = Math.floor(Math.min(minX, minY) * 10) / 10
    const globalMax = Math.ceil(Math.max(maxX, maxY) * 10) / 10
    
    // Use identical limits on both axes so the scatter is square in data space.
    const maxAxis = Math.max(Math.abs(globalMin), Math.abs(globalMax))
    const minAxis = -maxAxis
    

        new Chart(ctx, {
            type: 'scatter',
            options:{
                scales: {
                    x: {
                        type: 'linear',
                        position: 'bottom',
            min: minAxis,
            max: maxAxis,
                        ticks: {
                          stepSize: 1,
              font: { size: 8 },
              maxTicksLimit: 3,
                        }
                    },
                    y:{
                        type: 'linear',
                        position: 'bottom',
            min: minAxis,
            max: maxAxis,
                        ticks: {
                          stepSize: 1,
              font: { size: 8 },
              maxTicksLimit: 3,
                        }
                    }
                },
                plugins: {
                    legend: {
                        display : false
                    },
                },
                responsive: false,
        maintainAspectRatio: true,
        aspectRatio: 1,
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
  if((Date.now() - lastinteraction) > 10_000 && document.visibilityState !== "hidden"){
    lastinteraction = Date.now();
    location.reload();
  }
}, 1000);
