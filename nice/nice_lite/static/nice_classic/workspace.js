let lastinteraction = Date.now();

submitUpdate = (element, event) => {
  if(event.key === 'Enter') {
    element.form.submit()       
  }
}

enableWorkspaceRename = () => {
    const workspacerename = document.querySelector("#workspacerename")
    workspacerename.disabled = false
    workspacerename.focus()
}

enableWorkspaceDescription = () => {
    const workspacedescription = document.querySelector("#workspacedescription")
    workspacedescription.disabled = false
    workspacedescription.focus()
}

enableJobDescription = (element) => {
    const new_job_description = element.parentElement.querySelector('input[name="new_job_description"]')
    new_job_description.disabled = false
    new_job_description.focus()
}

deleteWorkspace = (element, workspaceid) => {
  const confirmed = confirm("Please confirm that you wish to delete workspace " + workspaceid + " and all jobs within. If this is the last workspace or dataset in a project, the project will also be deleted");
  if(confirmed){
    element.form.submit()
  }
}

window.addEventListener("load", () =>{

    const chart_config = {
        chart: {
            container    : "#tree-simple",
            connectors   : { 
                type  : "step",
                style : {
                    "stroke"       : "#022b3a",
                    "stroke-width" : 1
                }
            },
            hideRootNode : true,
            nodeAlign    : "BOTTOM",

        },
        nodeStructure: nodestructure
    };
    
    var my_chart = new Treant(chart_config);
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