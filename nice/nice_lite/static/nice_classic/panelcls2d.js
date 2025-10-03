toggleCls2d = (element, id) => {
    const deselected = element.querySelector(".deselected")
    let cls2ddeselectiontext = sessionStorage.getItem("cls2ddeselection")
    let cls2ddeselection = []
    if(cls2ddeselectiontext != null) cls2ddeselection = JSON.parse(cls2ddeselectiontext)
    deselected.classList.toggle("hidden")
    if(deselected.classList.contains("hidden")){
        cls2ddeselection = cls2ddeselection.filter(item => (item != id));
    }else{
        cls2ddeselection.push(id)
    }
    sessionStorage.setItem("cls2ddeselection", JSON.stringify(cls2ddeselection));
    updateCounts()
}

selectAbove = (element, id) => {
  let cls2ddeselection = []
  let threshold = true
  for(const cls2dcontainer of element.parentElement.querySelectorAll(".cls2dcontainer")){
    const deselected = cls2dcontainer.querySelector(".deselected")
    if(threshold){
      deselected.classList.add("hidden")
    }else{
      deselected.classList.remove("hidden")
      cls2ddeselection.push(Number(cls2dcontainer.dataset.idx))
    }
    if(cls2dcontainer == element){
      threshold = false
    }
  }
  hideMenu()
  sessionStorage.setItem("cls2ddeselection", JSON.stringify(cls2ddeselection));
  updateCounts()
}

selectBelow = (element, id) => {
  let cls2ddeselection = []
  let threshold = false
  for(const cls2dcontainer of element.parentElement.querySelectorAll(".cls2dcontainer")){
    if(cls2dcontainer == element){
      threshold = true
    }
    const deselected = cls2dcontainer.querySelector(".deselected")
    if(threshold){
      deselected.classList.add("hidden")
    }else{
      deselected.classList.remove("hidden")
      cls2ddeselection.push(Number(cls2dcontainer.dataset.idx))
    }
  }
  hideMenu()
  sessionStorage.setItem("cls2ddeselection", JSON.stringify(cls2ddeselection));
  updateCounts()
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

setSortDesc = (element) => {
    const sort_cls2d_asc = element.parentElement.parentElement.querySelector("[name='sort_cls2d_asc']")
    sort_cls2d_asc.value = false
    sort_cls2d_asc.form.submit()
}

setSortAsc = (element) => {
    const sort_cls2d_asc = element.parentElement.parentElement.querySelector("[name='sort_cls2d_asc']")
    sort_cls2d_asc.value = true
    sort_cls2d_asc.form.submit()
}

saveSelectionCls2D = (element) => {
  element.innerHTML = "saving ..."
  element.disabled  = true
  let cls2ddeselection = []
  let cls2ddeselectiontext = sessionStorage.getItem("cls2ddeselection")
  if(cls2ddeselectiontext != null) cls2ddeselection = JSON.parse(cls2ddeselectiontext)
  const deselected_cls2D = element.parentElement.querySelector("[name='deselected_cls2D']")
  deselected_cls2D.value = cls2ddeselection
  element.form.submit()
}

updateCounts = () => {
  const clustercount    = document.querySelector("#clustercount")
  const particlecount   = document.querySelector("#particlecount")
  const cls2dcontainers = document.querySelectorAll(".cls2dcontainer")
  let cls2ddeselectiontext = sessionStorage.getItem("cls2ddeselection")
  let cls2ddeselection = []
  if(cls2ddeselectiontext != null) cls2ddeselection = JSON.parse(cls2ddeselectiontext)
  clustercount.innerHTML = (cls2dcontainers.length - cls2ddeselection.length) + "/" + cls2dcontainers.length 
  nptcls     = 0
  nptcls_tot = 0
  for(const cls2dcontainer of cls2dcontainers){
    const deselected = cls2dcontainer.querySelector(".deselected")
    const pop = Number(cls2dcontainer.dataset.pop)
    const idx = Number(cls2dcontainer.dataset.idx)
    nptcls_tot = nptcls_tot + pop
    if(!cls2ddeselection.includes(idx)){
      nptcls = nptcls + pop
    }
  }
  particlecount.innerHTML = nptcls + "/" + nptcls_tot
}

window.addEventListener("load", () =>{
  let cls2ddeselectiontext = sessionStorage.getItem("cls2ddeselection")
  let cls2ddeselection = []
  if(cls2ddeselectiontext != null) cls2ddeselection = JSON.parse(cls2ddeselectiontext)
  for(const cls2dcontainer of  document.querySelectorAll(".cls2dcontainer")){
    const idx = Number(cls2dcontainer.dataset.idx)
    if(cls2ddeselection.includes(idx)){
      const deselected = cls2dcontainer.querySelector(".deselected")
      deselected.classList.remove("hidden")
    }   
  }
  updateCounts()
})

window.addEventListener("load", () =>{
  document.getElementById("loadinggauze").style.opacity = "0";
  setTimeout(function () {
   document.getElementById("loadinggauze").style.display = "none";
  }, 600);
})
