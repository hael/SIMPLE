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

