const fs = require('fs-extra')
const pug = require('pug')

class Browser{

  constructor() {
    this.browser = pug.compileFile(global.appPath + '/core/client/browser.pug')
  }

  getAll(arg) {
	var filter = (arg['filter']) ? arg['filter'] : ""
	var newfolder = (arg['newfolder']) ? arg['newfolder'] : false
    var view = {
      files : [],
      folders : [],
      path : arg['path'] ? arg['path'] : global.simplepath,
      filter : filter,
      newfolder : newfolder,
      buttons : arg['buttons']
    }
    var files
    return fs.readdir(view['path'])
    .then(contents => {
      files = contents
      var promises = []
      for (var object of contents){
        promises.push(fs.stat(view['path'] + "/" + object))
      }
      return Promise.all(promises.map(p => p.catch(e => e)))
    })
    .then(statarray => {
      for(var i in statarray){
		  if (statarray[i] instanceof Error == false){
				if(statarray[i].isFile() && files[i].charAt(0) != "."){
				  if(view['filter'] != "" && files[i].includes(view['filter'])){
					view['files'].push(files[i])
				  }else if(view['filter'] == ""){
					view['files'].push(files[i])
				  }
				}else if (statarray[i].isDirectory() && files[i].charAt(0) != "."){
				  view['folders'].push(files[i])
				}
		}
		}
      return({html : this.browser(view)})
    })
  }

  createFolder(arg){
	  return fs.mkdir(arg['path'] + "/" + arg['name'])
	  .then(() => {
		return ({})
	  })
  }

}

module.exports = Browser
