const pug = require('pug')
const fs = require('fs-extra')

class Tutorial {

  constructor(){
    this.listtutorialview = pug.compileFile(global.appPath + '/core/client/tutoriallist.pug')
    this.tutorialview = pug.compileFile(global.appPath + '/core/client/tutorialview.pug')
    this.tutorials = []
	fs.readdir(global.simplepath + '/tutorials/config')
	.then(contents => {
      for (var object of contents){
        if (object.includes(".tut")) {
          this.tutorials.push(object.replace(".tut",""))
        }
      }
    })
  }

  listTutorials(){
    return new Promise((resolve, reject) => {
      resolve(this.listtutorialview({tutorials:this.tutorials}))
    })
  }
  
  showTutorial(name){
    return fs.readJson(global.simplepath +'/tutorials/config/' + name + '.tut')
    .then(json => {
      return(this.tutorialview(json))
    })
  }

}

module.exports = Tutorial


