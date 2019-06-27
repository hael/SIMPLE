const sqlite3 = require(global.base + '/node_modules/sqlite3').verbose()

class Sqlite{

  constructor(){}

  sqlQuery(query){
    return new Promise((resolve,reject) => {
      var db = new sqlite3.Database(global.userdata +'/simple.sqlite')
      db.configure("busyTimeout", 10000)
      db.all(query, (err, rows) => {
        if(err == null){
          resolve(rows)
        } else {
          reject(err)
        }
      })
      db.close()
    })
  }

}

module.exports = new Sqlite()
