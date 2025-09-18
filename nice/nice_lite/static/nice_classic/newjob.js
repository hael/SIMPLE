toggleAdvanced = (element) => {
    element.classList.toggle("advanced")
    if(element.classList.contains("advanced")){
        element.querySelector("#advancedhidden").style.display  = "none"
        element.querySelector("#advancedvisible").style.display = "flex"
        for(var input of element.parentElement.parentElement.querySelectorAll(".advancedinput")){
            input.style.display = "flex"
        }
    }else{
        element.querySelector("#advancedvisible").style.display = "none"
        element.querySelector("#advancedhidden").style.display  = "flex"
        for(const input of element.parentElement.parentElement.querySelectorAll(".advancedinput")){
            input.style.display = "none"
        }
    }
}