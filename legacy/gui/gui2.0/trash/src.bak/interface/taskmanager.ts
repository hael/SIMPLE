class TaskManager {
	
	private viewpane
	private viewpanetitle
	private viewpanebody
	private target

	constructor() {
		this.viewpane = document.getElementById('viewer')
		this.viewpanetitle = document.getElementById('viewertitle')
		this.viewpanebody = document.getElementById('viewerbody')
		this.target = document.getElementById('target')
	}
	
	public show(){
		this.viewpanetitle.innerHTML = ""
		this.viewpanebody.innerHTML = ""
		taskselector.show()
	}
	
}

var taskmanager = new TaskManager()
