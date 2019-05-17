const pug = require('pug')
const sqlite = require('./sqlite')
const simpleexec = require('./simpleExec')

class Project{

  constructor(){
    this.projectselector = pug.compileFile(global.appPath + '/core/client/projectselector.pug')
    this.projecthistory = pug.compileFile(global.appPath + '/core/client/projecthistory.pug')
    this.createproject = pug.compileFile(global.appPath + '/core/client/createproject.pug')
    var query = "CREATE TABLE IF NOT EXISTS projects (id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, user TEXT, description TEXT, history TEXT UNIQUE, folder TEXT)"
    sqlite.sqlQuery(query)
  }

  getProjectSelector(arg){
    var query = "SELECT * FROM projects WHERE user='" + arg['user'] + "'"
    return sqlite.sqlQuery(query)
      .then((rows) => {
        var projects = {}
        projects['projects'] = rows
        return ({html : this.projectselector(projects)})
      })
  }

  getProjectHistory(arg){
    var query = "SELECT * FROM " + arg['history']
    return sqlite.sqlQuery(query)
    .then(rows => {
      var history = {}
      history['history'] = rows
      history['name'] = arg['name']
      history['description'] = arg['description']
      return ({html : this.projecthistory(history), data: rows})
    })
  }
  
  getCreateProject(arg){
	return new Promise((resolve, reject) => {
		resolve({html : this.createproject({})})
	})
  }

  createProject(arg){
    var projecttable = "t" + Math.floor(Math.random() * Math.floor(1000000))
    var query = "CREATE TABLE IF NOT EXISTS " + projecttable + " (id INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT, description TEXT, arguments TEXT, folder TEXT, pid TEXT, status TEXT, type TEXT, view TEXT, parent TEXT)"
    return sqlite.sqlQuery(query)
    .then(rows => {
       // var query = "INSERT INTO projects (name, user, description, history, folder) VALUES ('" + arg['keys']['keyname'] + "','" + arg['user'] + "','" + arg['keys']['keydescription'] + "','" + projecttable + "','" + arg['keys']['keyfolder'] + "')"
		var query = "INSERT INTO projects (name, user, description, history, folder) VALUES ('" + arg['keys']['keyname'] + "','" + arg['user'] + "','" + arg['keys']['keydescription'] + "','" + projecttable + "','" + arg['keys']['keyfolder'] + "/" +  arg['keys']['keyname'] + "')"
        return sqlite.sqlQuery(query)
    })
    .then(rows => {
        return simpleexec.createProject(arg['keys']['keyname'], arg['keys']['keyfolder'])
    })
    .then(output => {
		return({status:output.code})
	})
  }
  
  deleteProject(arg){
	  var query = "DROP TABLE " + arg['history']
	  return sqlite.sqlQuery(query)
	  .then(() =>{
			query = "DELETE FROM projects WHERE history='" + arg['history'] + "'"
			return sqlite.sqlQuery(query)
	  })
  }

}

module.exports = Project
