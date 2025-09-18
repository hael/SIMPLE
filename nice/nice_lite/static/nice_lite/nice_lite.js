changeProject = (form) => {
    selected_project_id = form.elements.selected_project_id.value
    if(selected_project_id == "new"){
        workspace_iframe = document.getElementById("workspace_iframe")
        workspace_iframe.src = "newproject"
    }else{
        form.submit()
    }
}

enableDatasetRename = () => {
    const datasetrename = document.querySelector("#datasetrename")
    datasetrename.disabled = false
    datasetrename.focus()
}

enableDatasetDescription = () => {
    const datasetdescription = document.querySelector("#datasetdescription")
    datasetdescription.disabled = false
    datasetdescription.focus()
}

enableStreamDescription = (element) => {
    const new_stream_description = element.parentElement.querySelector('input[name="new_stream_description"]')
    new_stream_description.disabled = false
    new_stream_description.focus()
}

changeDataset = (form) => {
    selected_dataset_id = form.elements.selected_dataset_id.value
    if(selected_dataset_id == "new"){
        form.submit()
    }else{
        document.cookie = "selected_dataset_id=" + selected_dataset_id
        workspace_iframe = document.getElementById("workspace_iframe")
        workspace_iframe.src = "dataset"
    }
}

window.addEventListener("load", () =>{
    /* loads the current dataset on page load */
    for(const element of document.getElementsByName("selected_dataset_id")){
        document.cookie = "selected_dataset_id=" + element.value
        workspace_iframe = document.getElementById("workspace_iframe")
        workspace_iframe.src = "dataset"
    }
},false);

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
                        size: 10
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
                        size: 10
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