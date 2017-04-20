
function statusPopup(data){
	var status = getArgument(data, 'status');
	if(status == "complete"){
		console.log("complete")
	}else{
		console.log("Error");
		console.log(data);
	}
}








function filesViewBack(){
	var filesviewtable = document.getElementById('filesviewtable');
	filesviewtable.style.display = "block";
	var filesviewdiv = document.getElementById('filesviewdiv');
	filesviewdiv.innerHTML = "";
}





function Null(data){
	console.log(data);
}
