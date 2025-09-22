changeProject = (form) => {
    selected_project_id = form.elements.selected_project_id.value
    if(selected_project_id == "new"){
        workspace_iframe = document.getElementById("workspace_iframe")
        workspace_iframe.src = form.dataset.newurl
    }else{
        form.submit()
    }
}

changeWorkspace = (form) => {
    selected_workspace_id = form.elements.selected_workspace_id.value
    if(selected_workspace_id == "new"){
        form.submit()
    }else{
        document.cookie = "selected_workspace_id=" + selected_workspace_id
        workspace_iframe = document.getElementById("workspace_iframe")
        workspace_iframe.src = "dataset"
    }
}

window.addEventListener("load", () =>{
    /* loads the current dataset on page load */
    for(const element of document.getElementsByName("selected_workspace_id")){
        document.cookie = "selected_workspace_id=" + element.value
        workspace_iframe = document.getElementById("workspace_iframe")
        workspace_iframe.src = "workspace"
    }
},false);
