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
var DataFormat = require("../../common/data-format");
var File = require("../../common/file");
var Coords = require("../algebra/coordinate");
var Box = require("../algebra/box");
var Logger = require("../utils/logger");
var state_1 = require("../state");
var server_config_1 = require("../../server-config");
var identify_1 = require("./identify");
var compose_1 = require("./compose");
var encode_1 = require("./encode");
function execute(params, outputProvider) {
    return __awaiter(this, void 0, void 0, function () {
        var start, guid, sourceFile, e_1;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    start = getTime();
                    state_1.State.pendingQueries++;
                    guid = generateUUID();
                    params.detail = Math.min(Math.max(0, params.detail | 0), server_config_1.default.limits.maxOutputSizeInVoxelCountByPrecisionLevel.length - 1);
                    Logger.log(guid, 'Info', "id=" + params.sourceId + ",encoding=" + (params.asBinary ? 'binary' : 'text') + ",detail=" + params.detail + "," + queryBoxToString(params.box));
                    sourceFile = void 0;
                    _a.label = 1;
                case 1:
                    _a.trys.push([1, 4, 5, 6]);
                    return [4 /*yield*/, File.openRead(params.sourceFilename)];
                case 2:
                    sourceFile = _a.sent();
                    return [4 /*yield*/, _execute(sourceFile, params, guid, outputProvider)];
                case 3:
                    _a.sent();
                    return [2 /*return*/, true];
                case 4:
                    e_1 = _a.sent();
                    Logger.error(guid, e_1);
                    return [2 /*return*/, false];
                case 5:
                    File.close(sourceFile);
                    Logger.log(guid, 'Time', Math.round(getTime() - start) + "ms");
                    state_1.State.pendingQueries--;
                    return [7 /*endfinally*/];
                case 6: return [2 /*return*/];
            }
        });
    });
}
exports.default = execute;
function getTime() {
    var t = process.hrtime();
    return t[0] * 1000 + t[1] / 1000000;
}
function generateUUID() {
    var d = getTime();
    var uuid = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function (c) {
        var r = (d + Math.random() * 16) % 16 | 0;
        d = Math.floor(d / 16);
        return (c === 'x' ? r : (r & 0x3 | 0x8)).toString(16);
    });
    return uuid;
}
function blockDomain(domain, blockSize) {
    var delta = Coords.fractional(blockSize * domain.delta[0], blockSize * domain.delta[1], blockSize * domain.delta[2]);
    return Coords.domain('Block', {
        origin: domain.origin,
        dimensions: domain.dimensions,
        delta: delta,
        sampleCount: Coords.sampleCounts(domain.dimensions, delta)
    });
}
function createSampling(header, index, dataOffset) {
    var sampling = header.sampling[index];
    var dataDomain = Coords.domain('Data', {
        origin: Coords.fractional(header.origin[0], header.origin[1], header.origin[2]),
        dimensions: Coords.fractional(header.dimensions[0], header.dimensions[1], header.dimensions[2]),
        delta: Coords.fractional(header.dimensions[0] / sampling.sampleCount[0], header.dimensions[1] / sampling.sampleCount[1], header.dimensions[2] / sampling.sampleCount[2]),
        sampleCount: sampling.sampleCount
    });
    return {
        index: index,
        rate: sampling.rate,
        byteOffset: sampling.byteOffset + dataOffset,
        dataDomain: dataDomain,
        blockDomain: blockDomain(dataDomain, header.blockSize)
    };
}
function createDataContext(file) {
    return __awaiter(this, void 0, void 0, function () {
        var _a, header, dataOffset, origin, dimensions;
        return __generator(this, function (_b) {
            switch (_b.label) {
                case 0: return [4 /*yield*/, DataFormat.readHeader(file)];
                case 1:
                    _a = _b.sent(), header = _a.header, dataOffset = _a.dataOffset;
                    origin = Coords.fractional(header.origin[0], header.origin[1], header.origin[2]);
                    dimensions = Coords.fractional(header.dimensions[0], header.dimensions[1], header.dimensions[2]);
                    return [2 /*return*/, {
                            file: file,
                            header: header,
                            spacegroup: Coords.spacegroup(header.spacegroup),
                            dataBox: { a: origin, b: Coords.add(origin, dimensions) },
                            sampling: header.sampling.map(function (s, i) { return createSampling(header, i, dataOffset); })
                        }];
            }
        });
    });
}
function createQuerySampling(data, sampling, queryBox) {
    var fractionalBox = Box.gridToFractional(Box.expandGridBox(Box.fractionalToGrid(queryBox, sampling.dataDomain), 1));
    var blocks = identify_1.default(data, sampling, fractionalBox);
    var ret = {
        sampling: sampling,
        fractionalBox: fractionalBox,
        gridDomain: Box.fractionalToDomain(fractionalBox, 'Query', sampling.dataDomain.delta),
        blocks: blocks
    };
    return ret;
}
function pickSampling(data, queryBox, forcedLevel, precision) {
    if (forcedLevel > 0) {
        return createQuerySampling(data, data.sampling[Math.min(data.sampling.length, forcedLevel) - 1], queryBox);
    }
    var sizeLimit = server_config_1.default.limits.maxOutputSizeInVoxelCountByPrecisionLevel[precision] || (2 * 1024 * 1024);
    for (var _i = 0, _a = data.sampling; _i < _a.length; _i++) {
        var s = _a[_i];
        var gridBox = Box.fractionalToGrid(queryBox, s.dataDomain);
        var approxSize = Box.volume(gridBox);
        if (approxSize <= sizeLimit) {
            var sampling = createQuerySampling(data, s, queryBox);
            if (sampling.blocks.length <= server_config_1.default.limits.maxRequestBlockCount) {
                return sampling;
            }
        }
    }
    return createQuerySampling(data, data.sampling[data.sampling.length - 1], queryBox);
}
function emptyQueryContext(data, params, guid) {
    return { kind: 'Empty', guid: guid, params: params, data: data };
}
function getQueryBox(data, queryBox) {
    switch (queryBox.kind) {
        case 'Cartesian': return Box.fractionalBoxReorderAxes(Box.cartesianToFractional(queryBox, data.spacegroup), data.header.axisOrder);
        case 'Fractional': return Box.fractionalBoxReorderAxes(queryBox, data.header.axisOrder);
        default: return data.dataBox;
    }
}
function allocateValues(domain, numChannels, valueType) {
    var values = [];
    for (var i = 0; i < numChannels; i++) {
        values[values.length] = DataFormat.createValueArray(valueType, domain.sampleVolume);
    }
    return values;
}
function createQueryContext(data, params, guid) {
    var inputQueryBox = getQueryBox(data, params.box);
    var queryBox;
    if (!data.header.spacegroup.isPeriodic) {
        if (!Box.areIntersecting(data.dataBox, inputQueryBox)) {
            return emptyQueryContext(data, params, guid);
        }
        queryBox = Box.intersect(data.dataBox, inputQueryBox);
    }
    else {
        queryBox = inputQueryBox;
    }
    var dimensions = Box.dimensions(queryBox);
    if (dimensions.some(function (d) { return isNaN(d); })) {
        throw "The query box is not defined.";
    }
    if (dimensions[0] * dimensions[1] * dimensions[2] > server_config_1.default.limits.maxFractionalBoxVolume) {
        throw "The query box volume is too big.";
    }
    var samplingInfo = pickSampling(data, queryBox, params.forcedSamplingLevel !== void 0 ? params.forcedSamplingLevel : 0, params.detail);
    if (samplingInfo.blocks.length === 0)
        return emptyQueryContext(data, params, guid);
    return {
        kind: 'Data',
        guid: guid,
        data: data,
        params: params,
        samplingInfo: samplingInfo,
        values: allocateValues(samplingInfo.gridDomain, data.header.channels.length, data.header.valueType)
    };
}
function _execute(file, params, guid, outputProvider) {
    return __awaiter(this, void 0, void 0, function () {
        var output, data, query, e_2, query;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    output = void 0;
                    _a.label = 1;
                case 1:
                    _a.trys.push([1, 5, 6, 7]);
                    return [4 /*yield*/, createDataContext(file)];
                case 2:
                    data = _a.sent();
                    query = createQueryContext(data, params, guid);
                    if (!(query.kind === 'Data')) return [3 /*break*/, 4];
                    // Step 3b: Compose the result data
                    return [4 /*yield*/, compose_1.default(query)];
                case 3:
                    // Step 3b: Compose the result data
                    _a.sent();
                    _a.label = 4;
                case 4:
                    // Step 4: Encode the result
                    output = outputProvider();
                    encode_1.default(query, output);
                    output.end();
                    return [3 /*break*/, 7];
                case 5:
                    e_2 = _a.sent();
                    query = { kind: 'Error', guid: guid, params: params, message: "" + e_2 };
                    try {
                        if (!output)
                            output = outputProvider();
                        encode_1.default(query, output);
                    }
                    catch (f) {
                        throw f;
                    }
                    throw e_2;
                case 6:
                    if (output)
                        output.end();
                    return [7 /*endfinally*/];
                case 7: return [2 /*return*/];
            }
        });
    });
}
function roundCoord(c) {
    return Math.round(100000 * c) / 100000;
}
function queryBoxToString(queryBox) {
    switch (queryBox.kind) {
        case 'Cartesian':
        case 'Fractional':
            var a = queryBox.a, b = queryBox.b;
            var r = roundCoord;
            return "box-type=" + queryBox.kind + ",box-a=(" + r(a[0]) + "," + r(a[1]) + "," + r(a[2]) + "),box-b=(" + r(b[0]) + "," + r(b[1]) + "," + r(b[2]) + ")";
        default:
            return "box-type=" + queryBox.kind;
    }
}
