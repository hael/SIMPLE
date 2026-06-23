toggleSingleMenu = () => {
    document.querySelector("#simpleprograms").style.display = "none"
    document.querySelector("#singleprograms").style.display = "flex"
}

toggleSimpleMenu = () => {
    document.querySelector("#singleprograms").style.display = "none"
    document.querySelector("#simpleprograms").style.display = "flex"
}

filterPrograms = (filterelement) => {
    const filtertext = filterelement.value.toLowerCase();
    for(const program of document.querySelectorAll(".program")){
        let display = false
        const disp = program.workspace.disp.toLowerCase()
        const desc = program.workspace.desc.toLowerCase()
        if(disp.includes(filtertext)) display = true
        if(desc.includes(filtertext)) display = true
        if(display){
            program.classList.remove("hidden")
        }else{
            program.classList.add("hidden")
        }
    }
}