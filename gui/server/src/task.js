const pug = require(global.base + '/node_modules/pug')
const running = require(global.base + '/node_modules/is-running')
const path = require('path')
const spawn = require(global.base + '/node_modules/child-process-promise').spawn
const fs = require(global.base + '/node_modules/fs-extra')

const filesystem = require('./fileSystem')
const sqlite = require('./sqlite')
const simpleexec = require('./simpleExec')

class Task{

  constructor() {
    this.taskselector = pug.compileFile(global.simplepath + '/gui_data/client/src/taskselector.pug')
    this.taskinput = pug.compileFile(global.simplepath + '/gui_data/client/src/taskinput.pug')
    var simplejson = require(global.userdata + '/simple_user_interface.json')
    this.simpleTasks = {}
    for (var command of simplejson){
      var commandkeys = Object.keys(command)
      var task = command[commandkeys[0]]
      task['pages'] = []
      for(var i = 1; i < commandkeys.length; i++){
        var page = {}
        page['title'] = commandkeys[i]
        page['keys'] = []
        for(var key of command[commandkeys[i]]){
          page['keys'].push(key)
        }
        if(page['keys'].length > 0){
          task['pages'].push(page)
        }
      }
      for (var page of task['pages']){
		 if(page['title'] == "computer controls"){
			  if(global.slurm && global.pbs){
				page['keys'].push({"key": "envqsys_name", "keytype": "multi", "options" : ["local", "slurm", "pbs"], "descr_short": "Submission system", "descr_long": "Queue system kind; (local|slurm|pbs)", "descr_placeholder": "", "required": true})
			  } else if (global.slurm){
				page['keys'].push({"key": "envqsys_name", "keytype": "multi", "options" : ["local", "slurm"], "descr_short": "Submission system", "descr_long": "Queue system kind; (local|slurm|pbs)", "descr_placeholder": "", "required": true})
			  } else if (global.pbs){
				page['keys'].push({"key": "envqsys_name", "keytype": "multi", "options" : ["local", "pbs"], "descr_short": "Submission system", "descr_long": "Queue system kind; (local|slurm|pbs)", "descr_placeholder": "", "required": true})
			  } else {
				page['keys'].push({"key": "envqsys_name", "keytype": "multi", "options" : ["local"], "descr_short": "Submission system", "descr_long": "Queue system kind; (local|slurm|pbs)", "descr_placeholder": "", "required": true})
			  }
			  if(global.slurm || global.pbs){
				page['keys'].push({"key": "envjob_memory_per_task", "keytype": "num", "descr_short": "Memory per part", "descr_long": "Memory per part; MB per part{1600}", "descr_placeholder": "", "required": false})
				page['keys'].push({"key": "envqsys_partition", "keytype": "str", "descr_short": " Name of SLURM/PBS partition", "descr_long": " Name of SLURM/PBS partition; give part name", "descr_placeholder": "", "required": false})
				page['keys'].push({"key": "envtime_per_image", "keytype": "num", "descr_short": "Time per image", "descr_long": "Time per image; in seconds{100}", "descr_placeholder": "", "required": false})
				page['keys'].push({"key": "envuser_email", "keytype": "str", "descr_short": "Your e-mail address", "descr_long": "Your e-mail address; e.g. myname@uni.edu", "descr_placeholder": "", "required": false})
				page['keys'].push({"key": "envqsys_qos", "keytype": "str", "descr_short": "Schedule priority", "descr_long": "Schedule priority; give priority", "descr_placeholder": "", "required": false})
				page['keys'].push({"key": "envqsys_reservation", "keytype": "str", "descr_short": "Name of reserved partition", "descr_long": "Name of reserved partition; give yourpart", "descr_placeholder": "", "required": false})
				page['keys'].push({"key": "envuser_account", "keytype": "str", "descr_short": "User account name in SLURM/PBS", "descr_long": "User account name in SLURM/PBS; e.g. Account084", "descr_placeholder": "", "required": false})
				page['keys'].push({"key": "envuser_project", "keytype": "str", "descr_short": "User project name in SLURM/PBS", "descr_long": "User project name in SLURM/PBS; e.g. Project001", "descr_placeholder": "", "required": false})
			  }
		}
	  }
      this.simpleTasks[commandkeys[0]] = task
	}
	this.createManualpickTask()
	this.attachFiletabGeneratorWidget()
	this.attachGuinierPlotWidget()
  }
  
  createManualpickTask(){
	  this.simpleTasks['manualpick'] = {
			"name": "manualpick",
			"descr_short": "Manual Particle Picking",
			"descr_long": "manually pick particles from micrographs",
			"executable": null,
			"pages" : []
		}
  }
  
  attachFiletabGeneratorWidget(){
    for(var widgetpage of this.simpleTasks['import_movies']['pages']){
      for(var widgetkeys of widgetpage['keys']){
        if(widgetkeys['key'] == "filetab"){
          widgetkeys['keytype'] = "widget"
          widgetkeys['widget'] = "filetabgeneratorwidget.view('keyfiletab')"
        }
      }
    }
  }
  
  attachGuinierPlotWidget(){
    for(var widgetpage of this.simpleTasks['postprocess']['pages']){
      for(var widgetkeys of widgetpage['keys']){
        if(widgetkeys['key'] == "bfac"){
          widgetkeys['keytype'] = "widget"
          widgetkeys['widget'] = "guinierplotwidget.view('keybfac')"
        }
      }
    }
  }

  getTaskSelector(arg) {
    return new Promise((resolve, reject) => {
      resolve({html : this.taskselector({tasks : this.simpleTasks})})
    })
  }

  getTaskSetup(arg){
	return (filesystem.getFolderContents(arg['projectfolder']))
	.then(contents =>{
		var inputs = []
		for(var folder of contents.folders){
			if(folder.match(/^\d/) && !folder.includes('_RELION')){
				inputs.push(folder)
			}
		}
		inputs.sort((a, b) => {
			var anum = a.substring(0, a.indexOf('_'))
			var bnum = b.substring(0, b.indexOf('_'))
			return bnum - anum
		});
		return({html : this.taskinput({task : this.simpleTasks[arg['task']], inputs:inputs, folder:arg['projectfolder'], name:arg['projectname']})})
	})
  }
  
  start(arg){
        arg['userdata'] = global.userdata
        if(global.exe.includes('simple_multiuser')){
		var promise = spawn("export", ["PKG_EXECPATH=''", '&&', 'simple_multiuser', 'execute', "'" + JSON.stringify(arg) + "'"], {detached: true, shell: true, cwd:arg['projectfolder']})
	}else{
		var promise = spawn(global.exe, ['execute', JSON.stringify(arg)], {detached: true, cwd:global.appPath})
	}
        var executeprocess = promise.childProcess
	promise.catch(error => {
		console.log('ERROR', error)
	})
        executeprocess.stdout.on('data', data => {
			console.log('DATA', data.toString())
		})
	executeprocess.stderr.on('data', data => {
                        console.log('DATA ERROR', data.toString())
                })
	//	return new Promise((resolve, reject) => {
	//		resolve({status:'running', jobid:})
	//	})
	return sqlite.sqlQuery("SELECT seq FROM sqlite_sequence WHERE name='" + arg['projecttable'] + "'")
	.then(rows => {
		if(rows != undefined){
			return({status:'running', jobid: 1})
		}else{
			return({status:'running', jobid: rows[0]['seq'] + 1})
		}
	})
  }
  
  kill(arg){
	  console.log("killing")
	  return sqlite.sqlQuery("UPDATE " + arg['projecttable'] + " SET pid='null',status='Killed' WHERE id=" + arg['id'])
	  .then(() => {
		  process.kill(arg['pid'], 'SIGKILL')
		  return
	  })
  }
  
  delete(arg){
 //   var query = "DELETE FROM " + arg['history'] + " WHERE id=" + arg['jobid']
	var query = "UPDATE " + arg['history'] + " SET pid='null',status='Deleted' WHERE id=" + arg['jobid']
    return sqlite.sqlQuery(query)
    .then(() => {
		if(arg['removefolder']){
			return fs.ensureDir(arg['projectfolder'] + '/Trash')
			.then(() => {
				return fs.move(arg['folder'], arg['projectfolder'] + '/Trash/' + path.basename(arg['folder']))
			})
			.then(() => {
				return fs.ensureDir(arg['folder'])
			})
			.then(() => {
				return({})
			})
		}else{
			return({})
		}
	})
  }
  
/*
  createNew(arg){
    var query = "INSERT into " + arg['projecttable'] + " (name, description, arguments, status, view, type) VALUES ('" + arg['name'] + "','" + arg['description'] + "','" + JSON.stringify(arg) + "','running', '" + JSON.stringify(arg['view']) + "', '" + arg['type'] + "')"
    var jobfolder
    var jobid
    return this.sqlQuery(query)
      .then((rows) => {
        var query = "SELECT seq FROM sqlite_sequence WHERE name='" + arg['projecttable'] + "'"
        return this.sqlQuery(query)
      })
      .then((rows) => {
        jobid = rows[0]['seq']
        jobfolder = arg['projectfolder'] + "/" + jobid + "_" + arg['type']
        fs.mkdirSync(jobfolder)
        query = "UPDATE " + arg['projecttable'] + " SET folder='" + jobfolder + "' WHERE id=" + jobid
        return this.sqlQuery(query)
      })
      .then((rows) => {
        return({jobfolder: jobfolder, jobid : jobid})
      })
  }

  updatePid(table, jobid, jobpid){
    var query = "UPDATE " + table + " SET pid='" + jobpid + "' WHERE id=" + jobid
    return this.sqlQuery(query)
  }

  updateStatus(table, jobid, status){
    var query = "UPDATE " + table + " SET status='" + status + "' WHERE id=" + jobid
    return this.sqlQuery(query)
  }

  delete(arg){
    var query = "DELETE FROM " + arg['history'] + " WHERE id=" + arg['jobid']
    return this.sqlQuery(query)
  }

  kill(arg){
    process.kill(arg['jobpid'])
    return this.updateStatus(arg['history'], arg['jobid'], "Killed")
  }

  updateJobs(){
    var query = "SELECT * FROM projects"
    var histories
		return this.sqlQuery(query)
    .then(rows => {
      histories = rows
      var promises = []
      for(var row of rows){
        var query = "SELECT * FROM " + row['history']
        promises.push(this.sqlQuery(query))
      }
      return Promise.all(promises)
    })
    .then(lines => {
      var promises = []
      for(var history in histories){
        for(var row of lines[history]){
          if(row['status'] == "running"  && row['pid'] != null && Number.isInteger(Number(row['pid'])) && !running(row['pid'])){
            promises.push(this.updateStatus(histories[history]['history'], row['id'], 'unknown'))
          } else if (row['status'] == "running"  && row['pid'] == null){
			promises.push(this.updateStatus(histories[history]['history'], row['id'], 'unknown'))
		  }
        }
      }
      return Promise.all(promises)
    })
  }*/
}

module.exports = Task
