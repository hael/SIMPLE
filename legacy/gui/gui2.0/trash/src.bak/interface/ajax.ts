function postAjax (request, success) {
	
	var xhr = new XMLHttpRequest()
	xhr.open('POST', '/', true)
	xhr.setRequestHeader("Content-type", "application/json")
	xhr.onreadystatechange = function() {
		if (xhr.readyState>3 && xhr.status==200){
			success(xhr.responseText)
		}
	};
	xhr.send(JSON.stringify(request))
}

function postAjaxPromise (request) {
	return fetch("/", {
		method: 'POST',
		body: JSON.stringify(request),
		headers: new Headers({
			'Content-Type': 'application/json',
		}),
		credentials: 'include'
	})
}
