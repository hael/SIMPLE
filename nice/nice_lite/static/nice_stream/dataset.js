let lastinteraction = Date.now();

submitUpdate = (element, event) => {
  if(event.key === 'Enter') {
    element.form.submit()       
  }
}

enableDatasetRename = () => {
    const datasetrename = document.querySelector("#datasetrename")
    datasetrename.disabled = false
    datasetrename.focus()
    lastinteraction = Date.now()
}

enableDatasetDescription = () => {
    const datasetdescription = document.querySelector("#datasetdescription")
    datasetdescription.disabled = false
    datasetdescription.focus()
    lastinteraction = Date.now()
}

enableStreamDescription = (element) => {
    const new_stream_description = element.parentElement.querySelector('input[name="new_stream_description"]')
    new_stream_description.disabled = false
    new_stream_description.focus()
    lastinteraction = Date.now()
}

deleteStream = (element, jobid) => {
  const confirmed = confirm("Please confirm that you wish to delete stream " + jobid);
  if(confirmed){
    element.form.submit()
  }
}

deleteDataset = (element, datasetid) => {
  const confirmed = confirm("Please confirm that you wish to delete dataset " + datasetid + " and all jobs within. If this is the last workspace or dataset in a project, the project will also be deleted");
  if(confirmed){
    element.form.submit()
  }
}

terminateStream = (element, jobid) => {
  const confirmed = confirm("Please confirm that you wish to terminate all processes in stream " + jobid);
  if(confirmed){
    element.form.submit()
  }
}

window.addEventListener("load", () =>{
    /* loads the current dataset on page load */
    for(const movies_pie_chart of document.getElementsByClassName("movies_pie_chart")){
        const ctx = movies_pie_chart.getContext("2d");
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
                  'imported',
                  'processed',
                  'rejected'
              ],
              datasets: [{
                  data: [Number(movies_pie_chart.dataset.imported), Number(movies_pie_chart.dataset.processed), Number(movies_pie_chart.dataset.rejected)],
                  backgroundColor: [
                    window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                    window.getComputedStyle(document.body).getPropertyValue('--color-nice4success'),
                    window.getComputedStyle(document.body).getPropertyValue('--color-nice4alert'),
                  ],
                  hoverOffset: 4
              }]
            }
        })
    }
},false);

window.addEventListener("load", () =>{
    /* loads the current dataset on page load */
    for(const particles_pie_chart of document.getElementsByClassName("particles_pie_chart")){
        const ctx = particles_pie_chart.getContext("2d");
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
                  'imported',
                  'accepted',
                  'rejected'
              ],
              datasets: [{
                  data: [Number(particles_pie_chart.dataset.imported) - Number(particles_pie_chart.dataset.accepted) - Number(particles_pie_chart.dataset.rejected), Number(particles_pie_chart.dataset.accepted), Number(particles_pie_chart.dataset.rejected)],
                  backgroundColor: [
                    window.getComputedStyle(document.body).getPropertyValue('--color-nice4header'),
                    window.getComputedStyle(document.body).getPropertyValue('--color-nice4success'),
                    window.getComputedStyle(document.body).getPropertyValue('--color-nice4alert'),
                  ],
                  hoverOffset: 4
              }]
            }
        })
    }
},false);

window.addEventListener("visibilitychange", (event) => {
  if(document.visibilityState !== "hidden"){
    location.reload();
  }
})

setInterval(function () {
  if((Date.now() - lastinteraction) > 10000 && document.visibilityState !== "hidden"){
    lastinteraction = Date.now();
    location.reload();
  }
}, 1000);