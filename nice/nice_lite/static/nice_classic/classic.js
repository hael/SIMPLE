changeProject = (element) => {
    selected_project_id = element.form.elements.selected_project_id.value
    if(selected_project_id == "new"){
        workspace_iframe = document.getElementById("workspace_iframe")
        workspace_iframe.src = element.form.dataset.newurl
        let selectselect = true
        for(const option of element.options){
            if(option.hasAttribute("selected")){
                option.selected = true
                selectselect = false
                break
            }
        }
        if(selectselect) element.options[0].selected = true
    }else{
        element.form.submit()
    }
}

changeWorkspace = (element) => {
    selected_workspace_id = element.form.elements.selected_workspace_id.value
    if(selected_workspace_id == "new"){
        element.form.submit()
    }else{
        document.cookie = "selected_workspace_id=" + selected_workspace_id
        workspace_iframe = document.getElementById("workspace_iframe")
        workspace_iframe.src = "workspace"
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
