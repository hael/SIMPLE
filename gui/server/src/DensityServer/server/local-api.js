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
var Api = require("./api");
var Coordinate = require("./algebra/coordinate");
var fs = require("fs");
var path = require("path");
function run(jobs) {
    return __awaiter(this, void 0, void 0, function () {
        var progress, started, _i, jobs_1, job, e_1, elapsed;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    progress = 0;
                    started = getTime();
                    _i = 0, jobs_1 = jobs;
                    _a.label = 1;
                case 1:
                    if (!(_i < jobs_1.length)) return [3 /*break*/, 7];
                    job = jobs_1[_i];
                    _a.label = 2;
                case 2:
                    _a.trys.push([2, 4, , 5]);
                    return [4 /*yield*/, query(job)];
                case 3:
                    _a.sent();
                    return [3 /*break*/, 5];
                case 4:
                    e_1 = _a.sent();
                    console.error(e_1);
                    return [3 /*break*/, 5];
                case 5:
                    progress++;
                    elapsed = (getTime() - started) / 1000;
                    console.log("[Progress] " + progress + "/" + jobs.length + " in " + elapsed.toFixed(2) + "s");
                    _a.label = 6;
                case 6:
                    _i++;
                    return [3 /*break*/, 1];
                case 7: return [2 /*return*/];
            }
        });
    });
}
exports.run = run;
function getTime() {
    var t = process.hrtime();
    return t[0] * 1000 + t[1] / 1000000;
}
function query(job) {
    return __awaiter(this, void 0, void 0, function () {
        var box, params, filename, res;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    if (job.query.kind.toLocaleLowerCase() === 'cell') {
                        box = { kind: 'Cell' };
                    }
                    else if (job.query.space === 'fractional') {
                        box = {
                            kind: 'Fractional',
                            a: Coordinate.fractional(job.query.bottomLeft[0], job.query.bottomLeft[1], job.query.bottomLeft[2]),
                            b: Coordinate.fractional(job.query.topRight[0], job.query.topRight[1], job.query.topRight[2]),
                        };
                    }
                    else {
                        box = {
                            kind: 'Cartesian',
                            a: Coordinate.cartesian(job.query.bottomLeft[0], job.query.bottomLeft[1], job.query.bottomLeft[2]),
                            b: Coordinate.cartesian(job.query.topRight[0], job.query.topRight[1], job.query.topRight[2]),
                        };
                    }
                    params = {
                        sourceFilename: job.source.filename,
                        sourceId: job.source.id,
                        asBinary: job.params.asBinary,
                        box: box,
                        detail: !job.params.detail ? 0 : job.params.detail,
                        forcedSamplingLevel: job.params.forcedSamplingLevel
                    };
                    if (!fs.existsSync(job.outputFolder)) {
                        makeDir(job.outputFolder);
                    }
                    filename = path.join(job.outputFolder, Api.getOutputFilename(job.source.name, job.source.id, params));
                    res = function () { return wrapFile(filename); };
                    return [4 /*yield*/, Api.queryBox(params, res)];
                case 1:
                    _a.sent();
                    return [2 /*return*/];
            }
        });
    });
}
function makeDir(path, root) {
    var dirs = path.split(/\/|\\/g), dir = dirs.shift();
    root = (root || '') + dir + '/';
    try {
        fs.mkdirSync(root);
    }
    catch (e) {
        if (!fs.statSync(root).isDirectory())
            throw new Error(e);
    }
    return !dirs.length || makeDir(dirs.join('/'), root);
}
function wrapFile(fn) {
    var w = {
        open: function () {
            if (this.opened)
                return;
            this.file = fs.openSync(fn, 'w');
            this.opened = true;
        },
        writeBinary: function (data) {
            this.open();
            fs.writeSync(this.file, new Buffer(data));
            return true;
        },
        writeString: function (data) {
            this.open();
            fs.writeSync(this.file, data);
            return true;
        },
        end: function () {
            if (!this.opened || this.ended)
                return;
            fs.close(this.file, function () { });
            this.ended = true;
        },
        file: 0,
        ended: false,
        opened: false
    };
    return w;
}
