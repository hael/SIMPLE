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
var File = require("../common/file");
var Data = require("./data-model");
var Sampling = require("./sampling");
var DataFormat = require("../common/data-format");
var fs = require("fs");
function pack(input, blockSize, isPeriodic, outputFilename) {
    return __awaiter(this, void 0, void 0, function () {
        var e_1;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    _a.trys.push([0, 2, , 3]);
                    return [4 /*yield*/, create(outputFilename, input, blockSize, isPeriodic)];
                case 1:
                    _a.sent();
                    return [3 /*break*/, 3];
                case 2:
                    e_1 = _a.sent();
                    console.error('[Error] ' + e_1);
                    return [3 /*break*/, 3];
                case 3: return [2 /*return*/];
            }
        });
    });
}
exports.default = pack;
function getTime() {
    var t = process.hrtime();
    return t[0] * 1000 + t[1] / 1000000;
}
function updateAllocationProgress(progress, progressDone) {
    var old = (100 * progress.current / progress.max).toFixed(0);
    progress.current += progressDone;
    var $new = (100 * progress.current / progress.max).toFixed(0);
    if (old !== $new) {
        process.stdout.write("\rAllocating...      " + $new + "%");
    }
}
/**
 * Pre allocate the disk space to be able to do "random" writes into the entire file.
 */
function allocateFile(ctx) {
    return __awaiter(this, void 0, void 0, function () {
        var totalByteSize, file, buffer, progress, written;
        return __generator(this, function (_a) {
            totalByteSize = ctx.totalByteSize, file = ctx.file;
            buffer = new Buffer(Math.min(totalByteSize, 8 * 1024 * 1024));
            progress = { current: 0, max: Math.ceil(totalByteSize / buffer.byteLength) };
            written = 0;
            while (written < totalByteSize) {
                written += fs.writeSync(file, buffer, 0, Math.min(totalByteSize - written, buffer.byteLength));
                updateAllocationProgress(progress, 1);
            }
            return [2 /*return*/];
        });
    });
}
function determineBlockSize(data, blockSize) {
    var extent = data.header.extent;
    var maxLayerSize = 1024 * 1024 * 1024;
    var valueCount = extent[0] * extent[1];
    if (valueCount * blockSize <= maxLayerSize)
        return blockSize;
    while (blockSize > 0) {
        blockSize -= 4;
        if (valueCount * blockSize <= maxLayerSize)
            return blockSize;
    }
    throw new Error('Could not determine a valid block size.');
}
function writeHeader(ctx) {
    return __awaiter(this, void 0, void 0, function () {
        var header;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    header = DataFormat.encodeHeader(Data.createHeader(ctx));
                    return [4 /*yield*/, File.writeInt(ctx.file, header.byteLength, 0)];
                case 1:
                    _a.sent();
                    return [4 /*yield*/, File.writeBuffer(ctx.file, 4, header)];
                case 2:
                    _a.sent();
                    return [2 /*return*/];
            }
        });
    });
}
function create(filename, sourceDensities, sourceBlockSize, isPeriodic) {
    return __awaiter(this, void 0, void 0, function () {
        var startedTime, files, channels_1, _i, sourceDensities_1, s, _a, _b, isOk, blockSize, _c, channels_2, ch, context, _d, channels_3, s, time, _e, files_1, f;
        return __generator(this, function (_f) {
            switch (_f.label) {
                case 0:
                    startedTime = getTime();
                    if (sourceBlockSize % 4 !== 0 || sourceBlockSize < 4) {
                        throw Error('Block size must be a positive number divisible by 4.');
                    }
                    if (!sourceDensities.length) {
                        throw Error('Specify at least one source density.');
                    }
                    process.stdout.write('Initializing... ');
                    files = [];
                    _f.label = 1;
                case 1:
                    _f.trys.push([1, , 10, 11]);
                    channels_1 = [];
                    _i = 0, sourceDensities_1 = sourceDensities;
                    _f.label = 2;
                case 2:
                    if (!(_i < sourceDensities_1.length)) return [3 /*break*/, 5];
                    s = sourceDensities_1[_i];
                    _b = (_a = channels_1).push;
                    return [4 /*yield*/, CCP4.open(s.name, s.filename)];
                case 3:
                    _b.apply(_a, [_f.sent()]);
                    _f.label = 4;
                case 4:
                    _i++;
                    return [3 /*break*/, 2];
                case 5:
                    isOk = channels_1.reduce(function (ok, s) { return ok && CCP4.compareHeaders(channels_1[0].header, s.header); }, true);
                    if (!isOk) {
                        throw new Error('Input file headers are not compatible (different grid, etc.).');
                    }
                    blockSize = determineBlockSize(channels_1[0], sourceBlockSize);
                    for (_c = 0, channels_2 = channels_1; _c < channels_2.length; _c++) {
                        ch = channels_2[_c];
                        CCP4.assignSliceBuffer(ch, blockSize);
                    }
                    return [4 /*yield*/, Sampling.createContext(filename, channels_1, blockSize, isPeriodic)];
                case 6:
                    context = _f.sent();
                    for (_d = 0, channels_3 = channels_1; _d < channels_3.length; _d++) {
                        s = channels_3[_d];
                        files.push(s.file);
                    }
                    files.push(context.file);
                    process.stdout.write('   done.\n');
                    console.log("Block size: " + blockSize);
                    // Step 2: Allocate disk space.        
                    process.stdout.write('Allocating...      0%');
                    return [4 /*yield*/, allocateFile(context)];
                case 7:
                    _f.sent();
                    process.stdout.write('\rAllocating...      done.\n');
                    // Step 3: Process and write the data 
                    process.stdout.write('Writing data...    0%');
                    return [4 /*yield*/, Sampling.processData(context)];
                case 8:
                    _f.sent();
                    process.stdout.write('\rWriting data...    done.\n');
                    // Step 4: Write the header at the start of the file.
                    // The header is written last because the sigma/min/max values are computed 
                    // during step 3.
                    process.stdout.write('Writing header...  ');
                    return [4 /*yield*/, writeHeader(context)];
                case 9:
                    _f.sent();
                    process.stdout.write('done.\n');
                    time = getTime() - startedTime;
                    console.log("[Done] " + time.toFixed(0) + "ms.");
                    return [3 /*break*/, 11];
                case 10:
                    for (_e = 0, files_1 = files; _e < files_1.length; _e++) {
                        f = files_1[_e];
                        File.close(f);
                    }
                    return [7 /*endfinally*/];
                case 11: return [2 /*return*/];
            }
        });
    });
}
