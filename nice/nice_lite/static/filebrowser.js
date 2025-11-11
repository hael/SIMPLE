let watchDouble = 0

const selectFile = (element, path, file, type) => {
    for(const button of document.querySelectorAll("button")){
        button.classList.remove("bg-nice4success")
    }
    element.classList.add("bg-nice4success")
    watchDouble += 1;
    setTimeout(()=>{
        if (watchDouble === 2) {
            element.form.submit()
        } else if (watchDouble === 1) {
            document.querySelector('#selectedpath').value = path + "/" + file
        }
        watchDouble = 0
    },200);
}

const selectDir = (element, path, dir, type) => {
    if(type == "dir"){
        for(const button of document.querySelectorAll("button")){
            button.classList.remove("bg-nice4success")
        }
        element.classList.add("bg-nice4success")
    }
    watchDouble += 1;
    setTimeout(()=>{
        if (watchDouble === 2) {
            element.form.submit()
        } else if (watchDouble === 1) {
            if(type == "dir"){ 
                document.querySelector('#selectedpath').value = path + "/" + dir
            }
        }
        watchDouble = 0
    },200);
}

const hideWindow = () => {
    parent.document.getElementsByName(window.name)[0].classList.add('hidden')
}

const select = () => {
    const selectedpath = document.querySelector('#selectedpath')
    const iframe = parent.document.getElementsByName(window.name)[0]
    const iframe_parent = iframe.parentElement
    iframe_parent.querySelector("input").value = selectedpath.value
    hideWindow()
}
