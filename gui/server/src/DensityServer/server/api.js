"use strict";
/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
var __assign = (this && this.__assign) || Object.assign || function(t) {
    for (var s, i = 1, n = arguments.length; i < n; i++) {
        s = arguments[i];
        for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
            t[p] = s[p];
    }
    return t;
};
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
var File = require("../common/file");
var execute_1 = require("./query/execute");
var Logger = require("./utils/logger");
var DataFormat = require("../common/data-format");
var server_config_1 = require("../server-config");
function getOutputFilename(source, id, _a) {
    var asBinary = _a.asBinary, box = _a.box, detail = _a.detail, forcedSamplingLevel = _a.forcedSamplingLevel;
    function n(s) { return (s || '').replace(/[ \n\t]/g, '').toLowerCase(); }
    function r(v) { return Math.round(10 * v) / 10; }
    var det = forcedSamplingLevel !== void 0
        ? "l" + forcedSamplingLevel
        : "d" + Math.min(Math.max(0, detail | 0), server_config_1.default.limits.maxOutputSizeInVoxelCountByPrecisionLevel.length - 1);
    var boxInfo = box.kind === 'Cell'
        ? 'cell'
        : (box.kind === 'Cartesian' ? 'cartn' : 'frac') + "_" + r(box.a[0]) + "_" + r(box.a[1]) + "_" + r(box.a[2]) + "_" + r(box.b[0]) + "_" + r(box.b[1]) + "_" + r(box.b[2]);
    return n(source) + "_" + n(id) + "-" + boxInfo + "_" + det + "." + (asBinary ? 'bcif' : 'cif');
}
exports.getOutputFilename = getOutputFilename;
/** Reads the header and includes information about available detail levels */
function getHeaderJson(filename, sourceId) {
    return __awaiter(this, void 0, void 0, function () {
        var header, _a, sampleCount, maxVoxelCount, precisions, availablePrecisions, _i, precisions_1, p, e_1;
        return __generator(this, function (_b) {
            switch (_b.label) {
                case 0:
                    Logger.logPlain('Header', sourceId);
                    _b.label = 1;
                case 1:
                    _b.trys.push([1, 3, , 4]);
                    if (!filename || !File.exists(filename)) {
                        Logger.errorPlain("Header " + sourceId, 'File not found.');
                        return [2 /*return*/, void 0];
                    }
                    _a = [{}];
                    return [4 /*yield*/, readHeader(filename, sourceId)];
                case 2:
                    header = __assign.apply(void 0, _a.concat([_b.sent()]));
                    sampleCount = header.sampling[0].sampleCount;
                    maxVoxelCount = sampleCount[0] * sampleCount[1] * sampleCount[2];
                    precisions = server_config_1.default.limits.maxOutputSizeInVoxelCountByPrecisionLevel
                        .map(function (maxVoxels, precision) { return ({ precision: precision, maxVoxels: maxVoxels }); });
                    availablePrecisions = [];
                    for (_i = 0, precisions_1 = precisions; _i < precisions_1.length; _i++) {
                        p = precisions_1[_i];
                        availablePrecisions.push(p);
                        if (p.maxVoxels > maxVoxelCount)
                            break;
                    }
                    header.availablePrecisions = availablePrecisions;
                    header.isAvailable = true;
                    return [2 /*return*/, JSON.stringify(header, null, 2)];
                case 3:
                    e_1 = _b.sent();
                    Logger.errorPlain("Header " + sourceId, e_1);
                    return [2 /*return*/, void 0];
                case 4: return [2 /*return*/];
            }
        });
    });
}
exports.getHeaderJson = getHeaderJson;
function queryBox(params, outputProvider) {
    return __awaiter(this, void 0, void 0, function () {
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0: return [4 /*yield*/, execute_1.default(params, outputProvider)];
                case 1: return [2 /*return*/, _a.sent()];
            }
        });
    });
}
exports.queryBox = queryBox;
function readHeader(filename, sourceId) {
    return __awaiter(this, void 0, void 0, function () {
        var file, header, e_2;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    file = void 0;
                    _a.label = 1;
                case 1:
                    _a.trys.push([1, 4, 5, 6]);
                    if (!filename)
                        return [2 /*return*/, void 0];
                    return [4 /*yield*/, File.openRead(filename)];
                case 2:
                    file = _a.sent();
                    return [4 /*yield*/, DataFormat.readHeader(file)];
                case 3:
                    header = _a.sent();
                    return [2 /*return*/, header.header];
                case 4:
                    e_2 = _a.sent();
                    Logger.errorPlain("Info " + sourceId, e_2);
                    return [2 /*return*/, void 0];
                case 5:
                    File.close(file);
                    return [7 /*endfinally*/];
                case 6: return [2 /*return*/];
            }
        });
    });
}
