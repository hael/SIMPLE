selectSegment = (element, key, jobid) => {
    const upperrow              = document.querySelector("#upperrow")
    const jobdisplaypanel       = document.querySelector("#jobdisplaypanel")
    const histogramdisplaypanel = document.querySelector("#histogramdisplaypanel")
    const plotdisplaypanel      = document.querySelector("#plotdisplaypanel")
    for(const jobdatasection of element.parentElement.querySelectorAll(".jobdatasection")){
        jobdatasection.classList.remove("text-nice4dark")
        jobdatasection.classList.add("text-nice4header")
    }
    element.classList.remove("text-nice4header")
    element.classList.add("text-nice4dark")
    if(key == "mic"){
        upperrow.classList.remove("hidden")
        jobdisplaypanel.src       = "/viewjobmicrographs/" + jobid
        histogramdisplaypanel.src = "/viewjobmicrographshistogram/" + jobid
        plotdisplaypanel.src      = "/viewjobmicrographsplot/" + jobid
    }else if(key == "ptcl2D"){
        upperrow.classList.remove("hidden")
        jobdisplaypanel.src       = "/viewjobcls2d/" + jobid
        histogramdisplaypanel.src = "/viewjobcls2dhistogram/" + jobid
    }else if(key == "cls2D"){
        upperrow.classList.remove("hidden")
        jobdisplaypanel.src       = "/viewjobcls2d/" + jobid
        histogramdisplaypanel.src = "/viewjobcls2dhistogram/" + jobid
        plotdisplaypanel.src      = "/viewjobcls2dplot/" + jobid
    }else if(key == "logs"){
        upperrow.classList.add("hidden")
        jobdisplaypanel.src       = "/viewjoblogs/" + jobid
        histogramdisplaypanel.src = "about:blank"
    }
}

window.addEventListener("load", () =>{
    sessionStorage.removeItem("micrographsdeselection")
    sessionStorage.removeItem("cls2ddeselection")
    if(document.querySelector(".jobdatasection[data-key='mic']") != undefined){
        document.querySelector(".jobdatasection[data-key='mic']").click()
    }else if(document.querySelector(".jobdatasection[data-key='ptcl2D']") != undefined){
        document.querySelector(".jobdatasection[data-key='ptcl2D']").click()
    }else if(document.querySelector(".jobdatasection[data-key='cls2D']") != undefined){
        document.querySelector(".jobdatasection[data-key='cls2D']").click()
    }else if(document.querySelector(".jobdatasection[data-key='logs']") != undefined){
        document.querySelector(".jobdatasection[data-key='logs']").click()
    }
})
