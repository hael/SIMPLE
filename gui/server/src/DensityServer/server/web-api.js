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
var Coords = require("./algebra/coordinate");
var documentation_1 = require("./documentation");
var server_config_1 = require("../server-config");
var Logger = require("./utils/logger");
var state_1 = require("./state");
function init(app) {
    function makePath(p) {
        return server_config_1.default.apiPrefix + '/' + p;
    }
    // Header
    app.get(makePath(':source/:id/?$'), function (req, res) { return getHeader(req, res); });
    // Box /:src/:id/box/:a1,:a2,:a3/:b1,:b2,:b3?text=0|1&space=cartesian|fractional
    app.get(makePath(':source/:id/box/:a1,:a2,:a3/:b1,:b2,:b3/?'), function (req, res) { return queryBox(req, res, getQueryParams(req, false)); });
    // Cell /:src/:id/cell/?text=0|1&space=cartesian|fractional
    app.get(makePath(':source/:id/cell/?'), function (req, res) { return queryBox(req, res, getQueryParams(req, true)); });
    app.get('*', function (req, res) {
        res.writeHead(200, { 'Content-Type': 'text/html; charset=utf-8' });
        res.end(documentation_1.default);
    });
}
exports.default = init;
function mapFile(type, id) {
    return server_config_1.default.mapFile(type || '', id || '');
}
function wrapResponse(fn, res) {
    var w = {
        do404: function () {
            if (!this.headerWritten) {
                res.writeHead(404);
                this.headerWritten = true;
            }
            this.end();
        },
        writeHeader: function (binary) {
            if (this.headerWritten)
                return;
            res.writeHead(200, {
                'Content-Type': binary ? 'application/octet-stream' : 'text/plain; charset=utf-8',
                'Access-Control-Allow-Origin': '*',
                'Access-Control-Allow-Headers': 'X-Requested-With',
                'Content-Disposition': "inline; filename=\"" + fn + "\""
            });
            this.headerWritten = true;
        },
        writeBinary: function (data) {
            if (!this.headerWritten)
                this.writeHeader(true);
            return res.write(new Buffer(data.buffer));
        },
        writeString: function (data) {
            if (!this.headerWritten)
                this.writeHeader(false);
            return res.write(data);
        },
        end: function () {
            if (this.ended)
                return;
            res.end();
            this.ended = true;
        },
        ended: false,
        headerWritten: false
    };
    return w;
}
function getSourceInfo(req) {
    return {
        filename: mapFile(req.params.source, req.params.id),
        id: req.params.source + "/" + req.params.id
    };
}
function validateSourndAndId(req, res) {
    if (!req.params.source || req.params.source.length > 32 || !req.params.id || req.params.source.id > 32) {
        res.writeHead(404);
        res.end();
        Logger.errorPlain("Query Box", 'Invalid source and/or id');
        return true;
    }
    return false;
}
function getHeader(req, res) {
    return __awaiter(this, void 0, void 0, function () {
        var headerWritten, _a, filename, id, header, e_1;
        return __generator(this, function (_b) {
            switch (_b.label) {
                case 0:
                    if (validateSourndAndId(req, res)) {
                        return [2 /*return*/];
                    }
                    headerWritten = false;
                    _b.label = 1;
                case 1:
                    _b.trys.push([1, 3, 4, 5]);
                    _a = getSourceInfo(req), filename = _a.filename, id = _a.id;
                    return [4 /*yield*/, Api.getHeaderJson(filename, id)];
                case 2:
                    header = _b.sent();
                    if (!header) {
                        res.writeHead(404);
                        return [2 /*return*/];
                    }
                    res.writeHead(200, {
                        'Content-Type': 'application/json; charset=utf-8',
                        'Access-Control-Allow-Origin': '*',
                        'Access-Control-Allow-Headers': 'X-Requested-With'
                    });
                    headerWritten = true;
                    res.write(header);
                    return [3 /*break*/, 5];
                case 3:
                    e_1 = _b.sent();
                    Logger.errorPlain("Header " + req.params.source + "/" + req.params.id, e_1);
                    if (!headerWritten) {
                        res.writeHead(404);
                    }
                    return [3 /*break*/, 5];
                case 4:
                    res.end();
                    return [7 /*endfinally*/];
                case 5: return [2 /*return*/];
            }
        });
    });
}
function getQueryParams(req, isCell) {
    var a = [+req.params.a1, +req.params.a2, +req.params.a3];
    var b = [+req.params.b1, +req.params.b2, +req.params.b3];
    var detail = Math.min(Math.max(0, (+req.query.detail) | 0), server_config_1.default.limits.maxOutputSizeInVoxelCountByPrecisionLevel.length - 1);
    var isCartesian = (req.query.space || '').toLowerCase() !== 'fractional';
    var box = isCell
        ? { kind: 'Cell' }
        : (isCartesian
            ? { kind: 'Cartesian', a: Coords.cartesian(a[0], a[1], a[2]), b: Coords.cartesian(b[0], b[1], b[2]) }
            : { kind: 'Fractional', a: Coords.fractional(a[0], a[1], a[2]), b: Coords.fractional(b[0], b[1], b[2]) });
    var asBinary = (req.query.encoding || '').toLowerCase() !== 'cif';
    var sourceFilename = mapFile(req.params.source, req.params.id);
    return {
        sourceFilename: sourceFilename,
        sourceId: req.params.source + "/" + req.params.id,
        asBinary: asBinary,
        box: box,
        detail: detail
    };
}
function queryBox(req, res, params) {
    return __awaiter(this, void 0, void 0, function () {
        var outputFilename, response, ok, e_2;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    if (validateSourndAndId(req, res)) {
                        return [2 /*return*/];
                    }
                    outputFilename = Api.getOutputFilename(req.params.source, req.params.id, params);
                    response = wrapResponse(outputFilename, res);
                    _a.label = 1;
                case 1:
                    _a.trys.push([1, 3, 4, 5]);
                    if (!params.sourceFilename) {
                        response.do404();
                        return [2 /*return*/];
                    }
                    return [4 /*yield*/, Api.queryBox(params, function () { return response; })];
                case 2:
                    ok = _a.sent();
                    if (!ok) {
                        response.do404();
                        return [2 /*return*/];
                    }
                    return [3 /*break*/, 5];
                case 3:
                    e_2 = _a.sent();
                    Logger.errorPlain("Query Box " + JSON.stringify(req.params || {}) + " | " + JSON.stringify(req.query || {}), e_2);
                    response.do404();
                    return [3 /*break*/, 5];
                case 4:
                    response.end();
                    queryDone();
                    return [7 /*endfinally*/];
                case 5: return [2 /*return*/];
            }
        });
    });
}
function queryDone() {
    if (state_1.State.shutdownOnZeroPending) {
        process.exit(0);
    }
}
