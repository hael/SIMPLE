let lastinteraction = Date.now();

changeProject = (form) => {
    selected_project_id = form.elements.selected_project_id.value
    if(selected_project_id == "new"){
        workspace_iframe = document.getElementById("workspace_iframe")
        workspace_iframe.src = form.dataset.newurl
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
