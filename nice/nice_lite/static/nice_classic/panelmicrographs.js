toggleMicrograph = (element, id) => {
    const deselected = element.querySelector(".deselected")
    let micrographsdeselectiontext = sessionStorage.getItem("micrographsdeselection")
    let micrographsdeselection = []
    if(micrographsdeselectiontext != null) micrographsdeselection = JSON.parse(micrographsdeselectiontext)
    deselected.classList.toggle("hidden")
    if(deselected.classList.contains("hidden")){
        micrographsdeselection = micrographsdeselection.filter(item => (item != id));
    }else{
        micrographsdeselection.push(id)
    }
    console.log(micrographsdeselection)
    sessionStorage.setItem("micrographsdeselection", JSON.stringify(micrographsdeselection));
}

selectAbove = (element, id) => {
  let micrographsdeselection = indices_post
  let threshold = true
  for(const miccontainer of element.parentElement.querySelectorAll(".micrographcontainer")){
    const deselected = miccontainer.querySelector(".deselected")
    if(threshold){
      deselected.classList.add("hidden")
    }else{
      deselected.classList.remove("hidden")
      micrographsdeselection.push(Number(miccontainer.dataset.idx))
    }
    if(miccontainer == element){
      threshold = false
    }
  }
  hideMenu()
  sessionStorage.setItem("micrographsdeselection", JSON.stringify(micrographsdeselection));
}

selectBelow = (element, id) => {
  let micrographsdeselection = indices_pre
  let threshold = false
  for(const miccontainer of element.parentElement.querySelectorAll(".micrographcontainer")){
    if(miccontainer == element){
      threshold = true
    }
    const deselected = miccontainer.querySelector(".deselected")
    if(threshold){
      deselected.classList.add("hidden")
    }else{
      deselected.classList.remove("hidden")
      micrographsdeselection.push(Number(miccontainer.dataset.idx))
    }

  }
  hideMenu()
  sessionStorage.setItem("micrographsdeselection", JSON.stringify(micrographsdeselection));
}

showMenu = (element, event) => {
  event.preventDefault()
  const selectmenu    = element.parentElement.querySelector("[name='selectmenu']")
  const selectmenubox = selectmenu.querySelector("[name='selectmenubox']")
  const selectabove   = selectmenu.querySelector("#selectabove")
  const selectbelow   = selectmenu.querySelector("#selectbelow")
  selectmenubox.style.top  = event.pageY + "px"
  selectmenubox.style.left = event.pageX + "px"
  selectmenu.style.display = "flex"
  if(selectabove != undefined) selectabove.onclick = () => {selectAbove(element)}
  if(selectbelow != undefined) selectbelow.onclick = () => {selectBelow(element)}
}

hideMenu = () => {
  for(const selectmenu of document.querySelectorAll("[name='selectmenu']")){
    const selectabove = selectmenu.querySelector("#selectabove")
    const selectbelow = selectmenu.querySelector("#selectbelow")
    selectmenu.style.display = "none"
    if(selectabove != undefined) selectabove.onclick = null
    if(selectbelow != undefined) selectbelow.onclick = null
  }
}

changePage = (element) => {
    const fromp = element.parentElement.parentElement.querySelector("[name='fromp']")
    const top   = element.parentElement.parentElement.querySelector("[name='top']")
    fromp.value = Number(element.options[element.selectedIndex].dataset.fromp)
    top.value   = Number(element.options[element.selectedIndex].dataset.top)
    element.form.submit()
}

setSortDesc = (element) => {
    const sort_micrographs_asc = element.parentElement.parentElement.querySelector("[name='sort_micrographs_asc']")
    const page  = element.parentElement.parentElement.querySelector("[name='page']")
    const fromp = element.parentElement.parentElement.querySelector("[name='fromp']")
    const top   = element.parentElement.parentElement.querySelector("[name='top']")
    page.value  = 1
    fromp.value = Number(page.options[0].dataset.fromp)
    top.value   = Number(page.options[0].dataset.top)
    sort_micrographs_asc.value = false
    sort_micrographs_asc.form.submit()
}

setSortAsc = (element) => {
    const sort_micrographs_asc = element.parentElement.parentElement.querySelector("[name='sort_micrographs_asc']")
    const page  = element.parentElement.parentElement.querySelector("[name='page']")
    const fromp = element.parentElement.parentElement.querySelector("[name='fromp']")
    const top   = element.parentElement.parentElement.querySelector("[name='top']")
    page.value  = 1
    fromp.value = Number(page.options[0].dataset.fromp)
    top.value   = Number(page.options[0].dataset.top)
    sort_micrographs_asc.value = true
    sort_micrographs_asc.form.submit()
}

saveSelectionMic = (element) => {
  element.innerHTML = "saving ..."
  element.disabled  = true
  let micrographsdeselection = []
  let micrographsdeselectiontext = sessionStorage.getItem("micrographsdeselection")
  if(micrographsdeselectiontext != null) micrographsdeselection = JSON.parse(micrographsdeselectiontext)
  const deselected_mic = element.parentElement.querySelector("[name='deselected_mic']")
  deselected_mic.value = micrographsdeselection
  element.form.submit()
}

window.addEventListener("load", () =>{
  let micrographsdeselectiontext = sessionStorage.getItem("micrographsdeselection")
  let micrographsdeselection = []
  if(micrographsdeselectiontext != null) micrographsdeselection = JSON.parse(micrographsdeselectiontext)
  for(const miccontainer of  document.querySelectorAll(".micrographcontainer")){
    const idx = Number(miccontainer.dataset.idx)
    if(micrographsdeselection.includes(idx)){
      const deselected = miccontainer.querySelector(".deselected")
      deselected.classList.remove("hidden")
    }   
  }
})

window.addEventListener("load", () =>{
  document.getElementById("loadinggauze").style.opacity = "0";
  setTimeout(function () {
   document.getElementById("loadinggauze").style.display = "none";
  }, 600);
})
