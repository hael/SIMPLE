linkParticleSet = (element) => {
  const loadinggauze = document.getElementById("loadinggauze")
  loadinggauze.innerHTML = "creating snapshot ..."
  loadinggauze.style.display = "flex";
  loadinggauze.style.opacity = 1;
  setTimeout(function () {
   document.getElementById("loadinggauze").style.display = "none";
   element.form.submit()
  }, 2000);
  return false
}

window.addEventListener("load", () =>{
  document.getElementById("loadinggauze").style.opacity = "0";
  setTimeout(function () {
   document.getElementById("loadinggauze").style.display = "none";
  }, 600);
})

setTimeout(function () {
  document.getElementById("loadinggauze").style.display = "flex";
  document.getElementById("loadinggauze").style.opacity = "1";
  setTimeout(function () {
   location.reload();
  }, 600);
}, 31000);