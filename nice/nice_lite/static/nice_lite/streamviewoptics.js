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

setTimeout(function () {
    location.reload();
}, 31000);

const data = {
  datasets: [{
    label: 'Scatter Dataset',
    data: [{
      x: -10,
      y: 0
    }, {
      x: 0,
      y: 10
    }, {
      x: 10,
      y: 5
    }, {
      x: 0.5,
      y: 5.5
    }],
    backgroundColor: 'rgb(255, 99, 132)'
  }],
};

setTimeout(function () {
    location.reload();
}, 31000);