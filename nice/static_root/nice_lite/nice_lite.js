changeProject = (form) => {
    selected_project_id = form.elements.selected_project_id.value
    if(selected_project_id == "new"){
        workspace_iframe = document.getElementById("workspace_iframe")
        workspace_iframe.src = "newproject"
    }else{
        form.submit()
    }
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
    for(const micrographs_pie_chart of document.getElementsByClassName("micrographs_pie_chart")){
        const ctx = micrographs_pie_chart.getContext("2d");
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
                  'queued',
                  'picked',
                  'rejected'
              ],
              datasets: [{
                  data: [Number(micrographs_pie_chart.dataset.imported) - Number(micrographs_pie_chart.dataset.processed) - Number(micrographs_pie_chart.dataset.rejected), Number(micrographs_pie_chart.dataset.processed), Number(micrographs_pie_chart.dataset.rejected)],
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