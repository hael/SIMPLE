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
var fs = require("fs");
var path = require("path");
var DataFormat = require("./data-format");
exports.IsNativeEndianLittle = new Uint16Array(new Uint8Array([0x12, 0x34]).buffer)[0] === 0x3412;
function openRead(filename) {
    return __awaiter(this, void 0, void 0, function () {
        var _this = this;
        return __generator(this, function (_a) {
            return [2 /*return*/, new Promise(function (res, rej) {
                    fs.open(filename, 'r', function (err, file) { return __awaiter(_this, void 0, void 0, function () {
                        return __generator(this, function (_a) {
                            if (err) {
                                rej(err);
                                return [2 /*return*/];
                            }
                            try {
                                res(file);
                            }
                            catch (e) {
                                fs.close(file);
                            }
                            return [2 /*return*/];
                        });
                    }); });
                })];
        });
    });
}
exports.openRead = openRead;
function readBuffer(file, position, sizeOrBuffer, size, byteOffset) {
    return new Promise(function (res, rej) {
        if (typeof sizeOrBuffer === 'number') {
            var buff = new Buffer(new ArrayBuffer(sizeOrBuffer));
            fs.read(file, buff, 0, sizeOrBuffer, position, function (err, bytesRead, buffer) {
                if (err) {
                    rej(err);
                    return;
                }
                res({ bytesRead: bytesRead, buffer: buffer });
            });
        }
        else {
            if (size === void 0) {
                rej('readBuffer: Specify size.');
                return;
            }
            fs.read(file, sizeOrBuffer, byteOffset ? +byteOffset : 0, size, position, function (err, bytesRead, buffer) {
                if (err) {
                    rej(err);
                    return;
                }
                res({ bytesRead: bytesRead, buffer: buffer });
            });
        }
    });
}
exports.readBuffer = readBuffer;
function writeBuffer(file, position, buffer, size) {
    return new Promise(function (res, rej) {
        fs.write(file, buffer, 0, size !== void 0 ? size : buffer.length, position, function (err, written) {
            if (err)
                rej(err);
            else
                res(written);
        });
    });
}
exports.writeBuffer = writeBuffer;
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
function exists(filename) {
    return fs.existsSync(filename);
}
exports.exists = exists;
function createFile(filename) {
    return new Promise(function (res, rej) {
        if (fs.existsSync(filename))
            fs.unlinkSync(filename);
        makeDir(path.dirname(filename));
        fs.open(filename, 'w', function (err, file) {
            if (err)
                rej(err);
            else
                res(file);
        });
    });
}
exports.createFile = createFile;
var __emptyFunc = function () { };
function close(file) {
    try {
        if (file !== void 0)
            fs.close(file, __emptyFunc);
    }
    catch (e) {
    }
}
exports.close = close;
var smallBuffer = new Buffer(8);
function writeInt(file, value, position) {
    return __awaiter(this, void 0, void 0, function () {
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    smallBuffer.writeInt32LE(value, 0);
                    return [4 /*yield*/, writeBuffer(file, position, smallBuffer, 4)];
                case 1:
                    _a.sent();
                    return [2 /*return*/];
            }
        });
    });
}
exports.writeInt = writeInt;
function getElementByteSize(type) {
    if (type === DataFormat.ValueType.Float32)
        return 4;
    if (type === DataFormat.ValueType.Int16)
        return 2;
    return 1;
}
function makeTypedArray(type, buffer) {
    if (type === DataFormat.ValueType.Float32)
        return new Float32Array(buffer);
    if (type === DataFormat.ValueType.Int16)
        return new Int16Array(buffer);
    return new Int8Array(buffer);
}
function createTypedArrayBufferContext(size, type) {
    var elementByteSize = getElementByteSize(type);
    var arrayBuffer = new ArrayBuffer(elementByteSize * size);
    var readBuffer = new Buffer(arrayBuffer);
    var valuesBuffer = exports.IsNativeEndianLittle ? arrayBuffer : new ArrayBuffer(elementByteSize * size);
    return {
        type: type,
        elementByteSize: elementByteSize,
        readBuffer: readBuffer,
        valuesBuffer: new Uint8Array(valuesBuffer),
        values: makeTypedArray(type, valuesBuffer)
    };
}
exports.createTypedArrayBufferContext = createTypedArrayBufferContext;
function flipByteOrder(source, target, byteCount, elementByteSize, offset) {
    for (var i = 0, n = byteCount; i < n; i += elementByteSize) {
        for (var j = 0; j < elementByteSize; j++) {
            target[offset + i + elementByteSize - j - 1] = source[offset + i + j];
        }
    }
}
function readTypedArray(ctx, file, position, count, valueOffset, littleEndian) {
    return __awaiter(this, void 0, void 0, function () {
        var byteCount, byteOffset;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0:
                    byteCount = ctx.elementByteSize * count;
                    byteOffset = ctx.elementByteSize * valueOffset;
                    return [4 /*yield*/, readBuffer(file, position, ctx.readBuffer, byteCount, byteOffset)];
                case 1:
                    _a.sent();
                    if (ctx.elementByteSize > 1 && ((littleEndian !== void 0 && littleEndian !== exports.IsNativeEndianLittle) || !exports.IsNativeEndianLittle)) {
                        // fix the endian 
                        flipByteOrder(ctx.readBuffer, ctx.valuesBuffer, byteCount, ctx.elementByteSize, byteOffset);
                    }
                    return [2 /*return*/, ctx.values];
            }
        });
    });
}
exports.readTypedArray = readTypedArray;
function ensureLittleEndian(source, target, byteCount, elementByteSize, offset) {
    if (exports.IsNativeEndianLittle)
        return;
    if (!byteCount || elementByteSize <= 1)
        return;
    flipByteOrder(source, target, byteCount, elementByteSize, offset);
}
exports.ensureLittleEndian = ensureLittleEndian;
