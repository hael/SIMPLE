const pug = require(global.base + '/node_modules/pug')
const fs = require(global.base + '/node_modules/fs-extra')

class Tutorial {

  constructor(){
    this.listtutorialview = pug.compileFile(global.simplepath + '/gui_data/client/src/tutoriallist.pug')
    this.tutorialview = pug.compileFile(global.simplepath + '/gui_data/client/src/tutorialview.pug')
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


