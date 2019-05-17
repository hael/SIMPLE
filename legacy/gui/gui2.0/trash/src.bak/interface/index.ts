fetch("/landing",{credentials: 'include'}).then(function(response) {
	if(response.status !== 200){
		console.log("Error")
		return
	}
	response.text().then(function(data){
		document.open()
		document.write(data)
		document.close()
	})
})
