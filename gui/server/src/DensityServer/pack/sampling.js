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
var CCP4 = require("./ccp4");
var Data = require("./data-model");
var File = require("../common/file");
var Downsampling = require("./downsampling");
var Writer = require("./writer");
var DataFormat = require("../common/data-format");
function createContext(filename, channels, blockSize, isPeriodic) {
    return __awaiter(this, void 0, void 0, function () {
        var header, samplingCounts, valueType, cubeBuffer, litteEndianCubeBuffer, ctx, _a, byteOffset, _i, _b, s;
        return __generator(this, function (_c) {
            switch (_c.label) {
                case 0:
                    header = channels[0].header;
                    samplingCounts = getSamplingCounts(channels[0].header.extent, blockSize);
                    valueType = CCP4.getValueType(header);
                    cubeBuffer = new Buffer(new ArrayBuffer(channels.length * blockSize * blockSize * blockSize * DataFormat.getValueByteSize(valueType)));
                    litteEndianCubeBuffer = File.IsNativeEndianLittle
                        ? cubeBuffer
                        : new Buffer(new ArrayBuffer(channels.length * blockSize * blockSize * blockSize * DataFormat.getValueByteSize(valueType)));
                    // The data can be periodic iff the extent is the same as the grid and origin is 0.
                    if (header.grid.some(function (v, i) { return v !== header.extent[i]; }) || header.origin.some(function (v) { return v !== 0; })) {
                        isPeriodic = false;
                    }
                    _a = {};
                    return [4 /*yield*/, File.createFile(filename)];
                case 1:
                    ctx = (_a.file = _c.sent(),
                        _a.isPeriodic = isPeriodic,
                        _a.channels = channels,
                        _a.valueType = valueType,
                        _a.blockSize = blockSize,
                        _a.cubeBuffer = cubeBuffer,
                        _a.litteEndianCubeBuffer = litteEndianCubeBuffer,
                        _a.kernel = { size: 5, coefficients: [1, 4, 6, 4, 1], coefficientSum: 16 },
                        _a.sampling = samplingCounts.map(function (__, i) { return createSampling(i, valueType, channels.length, samplingCounts, blockSize); }),
                        _a.dataByteOffset = 0,
                        _a.totalByteSize = 0,
                        _a.progress = { current: 0, max: 0 },
                        _a);
                    byteOffset = 0;
                    for (_i = 0, _b = ctx.sampling; _i < _b.length; _i++) {
                        s = _b[_i];
                        // Max progress = total number of blocks that need to be written.
                        ctx.progress.max += Data.samplingBlockCount(s, blockSize);
                        s.byteOffset = byteOffset;
                        byteOffset += s.byteSize;
                    }
                    ctx.dataByteOffset = 4 + DataFormat.encodeHeader(Data.createHeader(ctx)).byteLength;
                    ctx.totalByteSize = ctx.dataByteOffset + byteOffset;
                    return [2 /*return*/, ctx];
            }
        });
    });
}
exports.createContext = createContext;
function processData(ctx) {
    return __awaiter(this, void 0, void 0, function () {
        var channel, _i, _a, src;
        return __generator(this, function (_b) {
            switch (_b.label) {
                case 0:
                    channel = ctx.channels[0];
                    _b.label = 1;
                case 1:
                    if (!!channel.slices.isFinished) return [3 /*break*/, 7];
                    _i = 0, _a = ctx.channels;
                    _b.label = 2;
                case 2:
                    if (!(_i < _a.length)) return [3 /*break*/, 5];
                    src = _a[_i];
                    return [4 /*yield*/, CCP4.readSlices(src)];
                case 3:
                    _b.sent();
                    _b.label = 4;
                case 4:
                    _i++;
                    return [3 /*break*/, 2];
                case 5: return [4 /*yield*/, processSlices(ctx)];
                case 6:
                    _b.sent();
                    return [3 /*break*/, 1];
                case 7: return [2 /*return*/];
            }
        });
    });
}
exports.processData = processData;
/** Determine the suitable sampling rates for the input data */
function getSamplingCounts(baseSampleCount, blockSize) {
    var ret = [baseSampleCount];
    var prev = baseSampleCount;
    var hasSingleBoxSampling = false;
    while (true) {
        var next = [0, 0, 0];
        var max = 0;
        for (var i = 0; i < 3; i++) {
            var s = Math.floor((prev[i] + 1) / 2);
            if (s < 2)
                return ret;
            if (s > max)
                max = s;
            next[i] = s;
        }
        // no point in downsampling below the block size.
        if (max < blockSize) {
            if (hasSingleBoxSampling)
                return ret;
            hasSingleBoxSampling = true;
        }
        ret.push(next);
        prev = next;
    }
}
function createBlockBuffer(sampleCount, blockSize, valueType, numChannels) {
    var values = [];
    for (var i = 0; i < numChannels; i++)
        values[i] = DataFormat.createValueArray(valueType, sampleCount[0] * sampleCount[1] * blockSize);
    return {
        values: values,
        buffers: values.map(function (xs) { return new Buffer(xs.buffer); }),
        slicesWritten: 0
    };
}
function createDownsamplingBuffer(valueType, sourceSampleCount, targetSampleCount, numChannels) {
    var ret = [];
    for (var i = 0; i < numChannels; i++) {
        ret[ret.length] = {
            downsampleH: DataFormat.createValueArray(valueType, sourceSampleCount[1] * targetSampleCount[0]),
            downsampleHK: DataFormat.createValueArray(valueType, 5 * targetSampleCount[0] * targetSampleCount[1]),
            slicesWritten: 0,
            startSliceIndex: 0
        };
    }
    return ret;
}
function createSampling(index, valueType, numChannels, sampleCounts, blockSize) {
    var sampleCount = sampleCounts[index];
    var valuesInfo = [];
    for (var i = 0; i < numChannels; i++) {
        valuesInfo[valuesInfo.length] = {
            sum: 0.0,
            sqSum: 0.0,
            max: Number.NEGATIVE_INFINITY,
            min: Number.POSITIVE_INFINITY
        };
    }
    return {
        rate: 1 << index,
        sampleCount: sampleCount,
        blocks: createBlockBuffer(sampleCount, blockSize, valueType, numChannels),
        valuesInfo: valuesInfo,
        downsampling: index < sampleCounts.length - 1 ? createDownsamplingBuffer(valueType, sampleCount, sampleCounts[index + 1], numChannels) : void 0,
        byteOffset: 0,
        byteSize: numChannels * sampleCount[0] * sampleCount[1] * sampleCount[2] * DataFormat.getValueByteSize(valueType),
        writeByteOffset: 0
    };
}
function copyLayer(ctx, sliceIndex) {
    var channels = ctx.channels;
    var _a = ctx.sampling[0], blocks = _a.blocks, sampleCount = _a.sampleCount;
    var size = sampleCount[0] * sampleCount[1];
    var srcOffset = sliceIndex * size;
    var targetOffset = blocks.slicesWritten * size;
    for (var channelIndex = 0; channelIndex < channels.length; channelIndex++) {
        var src = channels[channelIndex].slices.values;
        var target = blocks.values[channelIndex];
        for (var i = 0; i < size; i++) {
            var v = src[srcOffset + i];
            target[targetOffset + i] = v;
        }
    }
    blocks.slicesWritten++;
}
function updateValuesInfo(sampling) {
    var blocks = sampling.blocks, sampleCount = sampling.sampleCount;
    var size = blocks.slicesWritten * sampleCount[0] * sampleCount[1];
    for (var channelIndex = 0; channelIndex < blocks.values.length; channelIndex++) {
        var values = blocks.values[channelIndex];
        var valuesInfo = sampling.valuesInfo[channelIndex];
        var sum = valuesInfo.sum, sqSum = valuesInfo.sqSum, max = valuesInfo.max, min = valuesInfo.min;
        for (var i = 0; i < size; i++) {
            var v = values[i];
            sum += v;
            sqSum += v * v;
            if (v > max)
                max = v;
            else if (v < min)
                min = v;
        }
        valuesInfo.sum = sum;
        valuesInfo.sqSum = sqSum;
        valuesInfo.max = max;
        valuesInfo.min = min;
    }
}
function shouldSamplingBeWritten(sampling, blockSize, isDataFinished) {
    if (isDataFinished)
        return sampling.blocks.slicesWritten > 0;
    return sampling.blocks.slicesWritten >= blockSize;
}
function writeBlocks(ctx, isDataFinished) {
    return __awaiter(this, void 0, void 0, function () {
        var _i, _a, s;
        return __generator(this, function (_b) {
            switch (_b.label) {
                case 0:
                    _i = 0, _a = ctx.sampling;
                    _b.label = 1;
                case 1:
                    if (!(_i < _a.length)) return [3 /*break*/, 4];
                    s = _a[_i];
                    if (!shouldSamplingBeWritten(s, ctx.blockSize, isDataFinished)) return [3 /*break*/, 3];
                    updateValuesInfo(s);
                    return [4 /*yield*/, Writer.writeBlockLayer(ctx, s)];
                case 2:
                    _b.sent();
                    _b.label = 3;
                case 3:
                    _i++;
                    return [3 /*break*/, 1];
                case 4: return [2 /*return*/];
            }
        });
    });
}
function processSlices(ctx) {
    return __awaiter(this, void 0, void 0, function () {
        var channel, sliceCount, i, isDataFinished;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    channel = ctx.channels[0];
                    sliceCount = channel.slices.sliceCount;
                    i = 0;
                    _a.label = 1;
                case 1:
                    if (!(i < sliceCount)) return [3 /*break*/, 5];
                    copyLayer(ctx, i);
                    Downsampling.downsampleLayer(ctx);
                    return [4 /*yield*/, writeBlocks(ctx, false)];
                case 2:
                    _a.sent();
                    isDataFinished = i === sliceCount - 1 && channel.slices.isFinished;
                    if (!isDataFinished) return [3 /*break*/, 4];
                    Downsampling.finalize(ctx);
                    return [4 /*yield*/, writeBlocks(ctx, true)];
                case 3:
                    _a.sent();
                    _a.label = 4;
                case 4:
                    i++;
                    return [3 /*break*/, 1];
                case 5: return [2 /*return*/];
            }
        });
    });
}
