"use strict";
/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : new P(function (resolve) { resolve(result.value); }).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
var __generator = (this && this.__generator) || function (thisArg, body) {
    var _ = { label: 0, sent: function() { if (t[0] & 1) throw t[1]; return t[1]; }, trys: [], ops: [] }, f, y, t, g;
    return g = { next: verb(0), "throw": verb(1), "return": verb(2) }, typeof Symbol === "function" && (g[Symbol.iterator] = function() { return this; }), g;
    function verb(n) { return function (v) { return step([n, v]); }; }
    function step(op) {
        if (f) throw new TypeError("Generator is already executing.");
        while (_) try {
            if (f = 1, y && (t = y[op[0] & 2 ? "return" : op[0] ? "throw" : "next"]) && !(t = t.call(y, op[1])).done) return t;
            if (y = 0, t) op = [0, t.value];
            switch (op[0]) {
                case 0: case 1: t = op; break;
                case 4: _.label++; return { value: op[1], done: false };
                case 5: _.label++; y = op[1]; op = [0]; continue;
                case 7: op = _.ops.pop(); _.trys.pop(); continue;
                default:
                    if (!(t = _.trys, t = t.length > 0 && t[t.length - 1]) && (op[0] === 6 || op[0] === 2)) { _ = 0; continue; }
                    if (op[0] === 3 && (!t || (op[1] > t[0] && op[1] < t[3]))) { _.label = op[1]; break; }
                    if (op[0] === 6 && _.label < t[1]) { _.label = t[1]; t = op; break; }
                    if (t && _.label < t[2]) { _.label = t[2]; _.ops.push(op); break; }
                    if (t[2]) _.ops.pop();
                    _.trys.pop(); continue;
            }
            op = body.call(thisArg, _);
        } catch (e) { op = [6, e]; y = 0; } finally { f = t = 0; }
        if (op[0] & 5) throw op[1]; return { value: op[0] ? op[1] : void 0, done: true };
    }
};
Object.defineProperty(exports, "__esModule", { value: true });
var LocalApi = require("./server/local-api");
var version_1 = require("./server/version");
var fs = require("fs");
console.log("DensityServer " + version_1.default + ", (c) 2016 - now, David Sehnal");
console.log();
function help() {
    var exampleJobs = [{
            source: {
                filename: "g:/test/mdb/xray-1tqn.mdb",
                name: 'xray',
                id: '1tqn',
            },
            query: {
                kind: 'box',
                space: 'cartesian',
                bottomLeft: [-42.996, -64.169, -45.335],
                topRight: [8.768, 15.316, 21.599]
            },
            params: {
                forcedSamplingLevel: 2,
                asBinary: true
            },
            outputFolder: 'g:/test/local-test'
        }, {
            source: {
                filename: "g:/test/mdb/emd-8116.mdb",
                name: 'em',
                id: '8116',
            },
            query: {
                kind: 'cell'
            },
            params: {
                detail: 4,
                asBinary: true
            },
            outputFolder: 'g:/test/local-test'
        }];
    console.log('Usage: node local jobs.json');
    console.log();
    console.log('Example jobs.json:');
    console.log(JSON.stringify(exampleJobs, null, 2));
}
function run() {
    return __awaiter(this, void 0, void 0, function () {
        var jobs;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    if (process.argv.length !== 3) {
                        help();
                        return [2 /*return*/];
                    }
                    try {
                        jobs = JSON.parse(fs.readFileSync(process.argv[2], 'utf-8'));
                    }
                    catch (e) {
                        console.log('Error:');
                        console.error(e);
                        return [2 /*return*/];
                    }
                    return [4 /*yield*/, LocalApi.run(jobs)];
                case 1:
                    _a.sent();
                    return [2 /*return*/];
            }
        });
    });
}
run();
