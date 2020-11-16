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
var File = require("../common/file");
var DataFormat = require("../common/data-format");
/** Converts a layer to blocks and writes them to the output file. */
function writeBlockLayer(ctx, sampling) {
    return __awaiter(this, void 0, void 0, function () {
        var nU, nV, startOffset, v, u, size;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    nU = Math.ceil(sampling.sampleCount[0] / ctx.blockSize);
                    nV = Math.ceil(sampling.sampleCount[1] / ctx.blockSize);
                    startOffset = ctx.dataByteOffset + sampling.byteOffset;
                    v = 0;
                    _a.label = 1;
                case 1:
                    if (!(v < nV)) return [3 /*break*/, 6];
                    u = 0;
                    _a.label = 2;
                case 2:
                    if (!(u < nU)) return [3 /*break*/, 5];
                    size = fillCubeBuffer(ctx, sampling, u, v);
                    return [4 /*yield*/, File.writeBuffer(ctx.file, startOffset + sampling.writeByteOffset, ctx.litteEndianCubeBuffer, size)];
                case 3:
                    _a.sent();
                    sampling.writeByteOffset += size;
                    updateProgress(ctx.progress, 1);
                    _a.label = 4;
                case 4:
                    u++;
                    return [3 /*break*/, 2];
                case 5:
                    v++;
                    return [3 /*break*/, 1];
                case 6:
                    sampling.blocks.slicesWritten = 0;
                    return [2 /*return*/];
            }
        });
    });
}
exports.writeBlockLayer = writeBlockLayer;
/** Fill a cube at position (u,v) with values from each of the channel */
function fillCubeBuffer(ctx, sampling, u, v) {
    var blockSize = ctx.blockSize, cubeBuffer = ctx.cubeBuffer;
    var sampleCount = sampling.sampleCount;
    var _a = sampling.blocks, buffers = _a.buffers, slicesWritten = _a.slicesWritten;
    var elementSize = DataFormat.getValueByteSize(ctx.valueType);
    var sizeH = sampleCount[0], sizeHK = sampleCount[0] * sampleCount[1];
    var offsetH = u * blockSize, offsetK = v * blockSize;
    var copyH = Math.min(blockSize, sampleCount[0] - offsetH) * elementSize, maxK = offsetK + Math.min(blockSize, sampleCount[1] - offsetK), maxL = slicesWritten;
    var writeOffset = 0;
    for (var _i = 0, buffers_1 = buffers; _i < buffers_1.length; _i++) {
        var src = buffers_1[_i];
        for (var l = 0; l < maxL; l++) {
            for (var k = offsetK; k < maxK; k++) {
                // copying the bytes direct is faster than using buffer.write* functions.
                var start = (l * sizeHK + k * sizeH + offsetH) * elementSize;
                src.copy(cubeBuffer, writeOffset, start, start + copyH);
                writeOffset += copyH;
            }
        }
    }
    // flip the byte order if needed.
    File.ensureLittleEndian(ctx.cubeBuffer, ctx.litteEndianCubeBuffer, writeOffset, elementSize, 0);
    return writeOffset;
}
function updateProgress(progress, progressDone) {
    var old = (100 * progress.current / progress.max).toFixed(0);
    progress.current += progressDone;
    var $new = (100 * progress.current / progress.max).toFixed(0);
    if (old !== $new) {
        process.stdout.write("\rWriting data...    " + $new + "%");
    }
}
