const fs = require(global.base + '/node_modules/fs-extra')

class FileSystem{

  constructor() {}

  getFolderContents(folder) {
	var contents
	var filestmp = []
	var files = []
	var folders = []
	
    return fs.readdir(folder)
    .then(contents => {
      filestmp = contents
      var promises = []
      for (var object of contents){
        promises.push(fs.stat(folder + "/" + object))
      }
     // return Promise.all(promises)
      return Promise.all(promises.map(p => p.catch(e => e)))
    })
    .then(statarray => {
      for(var i in statarray){
		 if (statarray[i] instanceof Error == false){
			if(statarray[i].isFile() && filestmp[i].charAt(0) != "."){
			  files.push(filestmp[i])
			}else if (statarray[i].isDirectory() && filestmp[i].charAt(0) != "."){
			  folders.push(filestmp[i])
			}
		}
      }
      return({files:files, folders:folders})
    })
  }
  
}

module.exports = new FileSystem()
