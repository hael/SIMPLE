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
var File = require("./file");
var Schema = require("./binary-schema");
var ValueType;
(function (ValueType) {
    ValueType.Float32 = 'float32';
    ValueType.Int8 = 'int8';
    ValueType.Int16 = 'int16';
})(ValueType = exports.ValueType || (exports.ValueType = {}));
var _schema;
(function (_schema) {
    var array = Schema.array, obj = Schema.obj, int = Schema.int, bool = Schema.bool, float = Schema.float, str = Schema.str;
    _schema.schema = obj([
        ['formatVersion', str],
        ['axisOrder', array(int)],
        ['origin', array(float)],
        ['dimensions', array(float)],
        ['spacegroup', obj([
                ['number', int],
                ['size', array(float)],
                ['angles', array(float)],
                ['isPeriodic', bool],
            ])],
        ['channels', array(str)],
        ['valueType', str],
        ['blockSize', int],
        ['sampling', array(obj([
                ['byteOffset', float],
                ['rate', int],
                ['valuesInfo', array(obj([
                        ['mean', float],
                        ['sigma', float],
                        ['min', float],
                        ['max', float]
                    ]))],
                ['sampleCount', array(int)]
            ]))]
    ]);
})(_schema || (_schema = {}));
var headerSchema = _schema.schema;
function getValueByteSize(type) {
    if (type === ValueType.Float32)
        return 4;
    if (type === ValueType.Int16)
        return 2;
    return 1;
}
exports.getValueByteSize = getValueByteSize;
function createValueArray(type, size) {
    switch (type) {
        case ValueType.Float32: return new Float32Array(new ArrayBuffer(4 * size));
        case ValueType.Int8: return new Int8Array(new ArrayBuffer(1 * size));
        case ValueType.Int16: return new Int16Array(new ArrayBuffer(2 * size));
    }
    throw Error(type + " is not a supported value format.");
}
exports.createValueArray = createValueArray;
function encodeHeader(header) {
    return Schema.encode(headerSchema, header);
}
exports.encodeHeader = encodeHeader;
function readHeader(file) {
    return __awaiter(this, void 0, void 0, function () {
        var buffer, headerSize, header;
        return __generator(this, function (_a) {
            switch (_a.label) {
                case 0: return [4 /*yield*/, File.readBuffer(file, 0, 4 * 4096)];
                case 1:
                    buffer = (_a.sent()).buffer;
                    headerSize = buffer.readInt32LE(0);
                    if (!(headerSize > buffer.byteLength - 4)) return [3 /*break*/, 3];
                    return [4 /*yield*/, File.readBuffer(file, 0, headerSize + 4)];
                case 2:
                    buffer = (_a.sent()).buffer;
                    _a.label = 3;
                case 3:
                    header = Schema.decode(headerSchema, buffer, 4);
                    return [2 /*return*/, { header: header, dataOffset: headerSize + 4 }];
            }
        });
    });
}
exports.readHeader = readHeader;
