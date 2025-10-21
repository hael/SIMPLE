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
    sessionStorage.setItem("micrographsdeselection", JSON.stringify(micrographsdeselection));
    updateCounts()
}

selectAbove = (element, id) => {
  let micrographsdeselection = [...indices_post]
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
  updateCounts()
}

selectBelow = (element, id) => {
  let micrographsdeselection = [...indices_pre]
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

updateCounts = () => {
  const micscount            = document.querySelector("#micscount")
  const micrographcontainers = document.querySelectorAll(".micrographcontainer")
  let micrographsdeselection = []
  let micrographsdeselectiontext = sessionStorage.getItem("micrographsdeselection")
  if(micrographsdeselectiontext != null) micrographsdeselection = JSON.parse(micrographsdeselectiontext)
  micscount.innerHTML = (indices_pre.length + indices_post.length + micrographcontainers.length - micrographsdeselection.length) + "/" + (indices_pre.length + indices_post.length + micrographcontainers.length)
}

drawBoxes = () => {
  const boxestoggle = document.querySelector('#boxestoggle')
  if(boxestoggle == undefined) return
  const display = boxestoggle.checked
  for(const micrographcontainer of document.querySelectorAll(".micrographcontainer")){
    const canvas = micrographcontainer.querySelector('.boxes_overlay')
    if(canvas == undefined) continue
    const micimg = micrographcontainer.querySelector('.micimg')
    if(micimg == undefined) continue
    if(!micrographcontainer.dataset.hasOwnProperty('boxes')) continue
    if(!micrographcontainer.dataset.hasOwnProperty('xdim'))  continue
    if(!micrographcontainer.dataset.hasOwnProperty('ydim'))  continue
    const boxes = JSON.parse(micrographcontainer.dataset.boxes.replaceAll("'", '"'))
    const xdim  = Number(micrographcontainer.dataset.xdim)
    const ydim  = Number(micrographcontainer.dataset.ydim)
    const imgrect   = micimg.getBoundingClientRect()
    const imgwidth  = imgrect.width / 2 // divide by 2 as has PS
    const imgheight = imgrect.height
    canvas.width  = imgwidth
    canvas.height = imgheight
    canvas.style.marginLeft = imgwidth + "px"
    const scale = imgwidth  / xdim
    // account for rectangular images
    var yoffset = (canvas.height - (scale * ydim)) / 2
    const ctx = canvas.getContext("2d");
    ctx.clearRect(0, 0, canvas.width, canvas.height)
    if(display){
      for(const box of boxes){
        ctx.strokeStyle = "yellow";
        ctx.beginPath();
        ctx.arc(box["x"] * scale, (box["y"] * scale) + yoffset, 0.5, 0, 2 * Math.PI);
        ctx.stroke();
      }
    }
  }
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
  updateCounts()
})

window.addEventListener("load", () =>{
  document.getElementById("loadinggauze").style.opacity = "0";
  setTimeout(function () {
   document.getElementById("loadinggauze").style.display = "none";
  }, 600);
})

window.addEventListener("load", () =>{
  drawBoxes()
})