"use strict";
/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
Object.defineProperty(exports, "__esModule", { value: true });
var express = require("express");
var compression = require("compression");
var web_api_1 = require("./server/web-api");
var version_1 = require("./server/version");
var server_config_1 = require("./server-config");
var Logger = require("./server/utils/logger");
var state_1 = require("./server/state");
function setupShutdown() {
    if (server_config_1.default.shutdownParams.timeoutVarianceMinutes > server_config_1.default.shutdownParams.timeoutMinutes) {
        Logger.logPlain('Server', 'Shutdown timeout variance is greater than the timer itself, ignoring.');
    }
    else {
        var tVar = 0;
        if (server_config_1.default.shutdownParams.timeoutVarianceMinutes > 0) {
            tVar = 2 * (Math.random() - 0.5) * server_config_1.default.shutdownParams.timeoutVarianceMinutes;
        }
        var tMs = (server_config_1.default.shutdownParams.timeoutMinutes + tVar) * 60 * 1000;
        console.log("----------------------------------------------------------------------------");
        console.log("  The server will shut down in " + Logger.formatTime(tMs) + " to prevent slow performance.");
        console.log("  Please make sure a daemon is running that will automatically restart it.");
        console.log("----------------------------------------------------------------------------");
        console.log();
        setTimeout(function () {
            if (state_1.State.pendingQueries > 0) {
                state_1.State.shutdownOnZeroPending = true;
            }
            else {
                Logger.logPlain('Server', "Shut down due to timeout.");
                process.exit(0);
            }
        }, tMs);
    }
}
var port = process.env.port || server_config_1.default.defaultPort;
var app = express();
app.use(compression({ level: 6, memLevel: 9, chunkSize: 16 * 16384, filter: function () { return true; } }));
web_api_1.default(app);
app.listen(port);
console.log("DensityServer " + version_1.default + ", (c) 2016 - now, David Sehnal");
console.log("");
console.log("The server is running on port " + port + ".");
console.log("");
if (server_config_1.default.shutdownParams && server_config_1.default.shutdownParams.timeoutMinutes > 0) {
    setupShutdown();
}
