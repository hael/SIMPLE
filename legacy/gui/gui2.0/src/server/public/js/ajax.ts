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
