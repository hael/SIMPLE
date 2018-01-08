function getAjax (url, success) {
    var xhr = new XMLHttpRequest();
    xhr.open('GET', url, true);
    xhr.onreadystatechange = function() {
        if (xhr.readyState>3 && xhr.status==200){
			success(xhr.responseText);
		}
    };
    xhr.send();
}

function get3DViewData () {
 getAjax('JSONhandler?function=3dview&' + window.location.toString().split('?')[1], function(data){view3DInit(data)});
}

function view3DInit (data) {
	var JSONdata = JSON.parse(data);
	var iterations = JSONdata.iterations;
	var snapshots = document.getElementById('snapshots');

	for(var state = 1; state < 99; state++){
		for(var i = 0; i < iterations.length; i++){
				if(iterations[i].state == state){	
					if(!document.getElementById('state'+state)){
						var p = document.createElement('pre');
						p.innerHTML = "State " + state + "\n";
						p.id = 'state'+state;
						snapshots.appendChild(p);
					}
				}
			
		}
	}
	for(var iteration = 1; iteration < 99; iteration++){
		for(var i = 0; i < iterations.length; i++){
				if(iterations[i].iteration == iteration){	
					var p = document.getElementById('state'+iterations[i].state)
					p.innerHTML += "\tIteration " + iteration + "\n";
					p.innerHTML += "\t\tRESOLUTION AT FSC=0.500 " + iterations[i].fsc50  + "\n";
					p.innerHTML += "\t\tRESOLUTION AT FSC=0.143 " + iterations[i].fsc143 + "\n";
				}
		}
	}
}

function showProjectHistory(){
	var mainpaneiframe = parent.document.getElementById('mainpaneiframe');
	mainpaneiframe.src = "projecthistory.html"
}

get3DViewData();
