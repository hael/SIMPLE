"use strict";
/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
Object.defineProperty(exports, "__esModule", { value: true });
var CCP4 = require("./ccp4");
var FORMAT_VERSION = '1.0.0';
function createHeader(ctx) {
    var header = ctx.channels[0].header;
    var grid = header.grid;
    function normalize(data) {
        return [data[0] / grid[0], data[1] / grid[1], data[2] / grid[2]];
    }
    return {
        formatVersion: FORMAT_VERSION,
        valueType: CCP4.getValueType(header),
        blockSize: ctx.blockSize,
        axisOrder: header.axisOrder,
        origin: normalize(header.origin),
        dimensions: normalize(header.extent),
        spacegroup: { number: header.spacegroupNumber, size: header.cellSize, angles: header.cellAngles, isPeriodic: ctx.isPeriodic },
        channels: ctx.channels.map(function (c) { return c.header.name; }),
        sampling: ctx.sampling.map(function (s) {
            var N = s.sampleCount[0] * s.sampleCount[1] * s.sampleCount[2];
            var valuesInfo = [];
            for (var _i = 0, _a = s.valuesInfo; _i < _a.length; _i++) {
                var _b = _a[_i], sum = _b.sum, sqSum = _b.sqSum, min = _b.min, max = _b.max;
                var mean = sum / N;
                var sigma = Math.sqrt(Math.max(0, sqSum / N - mean * mean));
                valuesInfo.push({ mean: mean, sigma: sigma, min: min, max: max });
            }
            return {
                byteOffset: s.byteOffset,
                rate: s.rate,
                valuesInfo: valuesInfo,
                sampleCount: s.sampleCount,
            };
        })
    };
}
exports.createHeader = createHeader;
function samplingBlockCount(sampling, blockSize) {
    return sampling.sampleCount.map(function (c) { return Math.ceil(c / blockSize); }).reduce(function (c, v) { return c * v; }, 1);
}
exports.samplingBlockCount = samplingBlockCount;
