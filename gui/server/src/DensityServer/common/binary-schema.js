"use strict";
/*
 * Copyright (c) 2017 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
Object.defineProperty(exports, "__esModule", { value: true });
var CIF = require("../lib/cif-tools");
var MsgPack = CIF.Binary.MessagePack;
exports.bool = { kind: 'bool' };
exports.int = { kind: 'int' };
exports.float = { kind: 'float' };
exports.str = { kind: 'string' };
// tslint:disable-next-line:array-type
function array(element) { return { kind: 'array', element: element }; }
exports.array = array;
function obj(schema) {
    return {
        kind: 'object',
        props: schema.map(function (s) { return ({
            element: s[1],
            prop: s[0]
        }); })
    };
}
exports.obj = obj;
function byteCount(e, src) {
    var size = 0;
    switch (e.kind) {
        case 'bool':
            size += 1;
            break;
        case 'int':
            size += 4;
            break;
        case 'float':
            size += 8;
            break;
        case 'string':
            size += 4 + MsgPack.utf8ByteCount(src);
            break;
        case 'array': {
            size += 4; // array length
            for (var _i = 0, src_1 = src; _i < src_1.length; _i++) {
                var x = src_1[_i];
                size += byteCount(e.element, x);
            }
            break;
        }
        case 'object': {
            for (var _a = 0, _b = e.props; _a < _b.length; _a++) {
                var p = _b[_a];
                size += byteCount(p.element, src[p.prop]);
            }
            break;
        }
    }
    return size;
}
function writeElement(e, buffer, src, offset) {
    switch (e.kind) {
        case 'bool':
            buffer.writeInt8(src ? 1 : 0, offset);
            offset += 1;
            break;
        case 'int':
            buffer.writeInt32LE(src | 0, offset);
            offset += 4;
            break;
        case 'float':
            buffer.writeDoubleLE(+src, offset);
            offset += 8;
            break;
        case 'string': {
            var val = '' + src;
            var size = MsgPack.utf8ByteCount(val);
            buffer.writeInt32LE(size, offset);
            offset += 4; // str len
            var str_1 = new Uint8Array(size);
            MsgPack.utf8Write(str_1, 0, val);
            for (var _i = 0, _a = str_1; _i < _a.length; _i++) {
                var b = _a[_i];
                buffer.writeUInt8(b, offset);
                offset++;
            }
            break;
        }
        case 'array': {
            buffer.writeInt32LE(src.length, offset);
            offset += 4; // array length
            for (var _b = 0, src_2 = src; _b < src_2.length; _b++) {
                var x = src_2[_b];
                offset = writeElement(e.element, buffer, x, offset);
            }
            break;
        }
        case 'object': {
            for (var _c = 0, _d = e.props; _c < _d.length; _c++) {
                var p = _d[_c];
                offset = writeElement(p.element, buffer, src[p.prop], offset);
            }
            break;
        }
    }
    return offset;
}
function write(element, src) {
    var size = byteCount(element, src);
    var buffer = new Buffer(size);
    writeElement(element, buffer, src, 0);
    return buffer;
}
function encode(element, src) {
    return write(element, src);
}
exports.encode = encode;
function decodeElement(e, buffer, offset, target) {
    switch (e.kind) {
        case 'bool':
            target.value = !!buffer.readInt8(offset);
            offset += 1;
            break;
        case 'int':
            target.value = buffer.readInt32LE(offset);
            offset += 4;
            break;
        case 'float':
            target.value = buffer.readDoubleLE(offset);
            offset += 8;
            break;
        case 'string': {
            var size = buffer.readInt32LE(offset);
            offset += 4; // str len
            var str_2 = new Uint8Array(size);
            for (var i = 0; i < size; i++) {
                str_2[i] = buffer.readUInt8(offset);
                offset++;
            }
            target.value = MsgPack.utf8Read(str_2, 0, size);
            break;
        }
        case 'array': {
            var array_1 = [];
            var count = buffer.readInt32LE(offset);
            var element = { value: void 0 };
            offset += 4;
            for (var i = 0; i < count; i++) {
                offset = decodeElement(e.element, buffer, offset, element);
                array_1.push(element.value);
            }
            target.value = array_1;
            break;
        }
        case 'object': {
            var t = Object.create(null);
            var element = { value: void 0 };
            for (var _i = 0, _a = e.props; _i < _a.length; _i++) {
                var p = _a[_i];
                offset = decodeElement(p.element, buffer, offset, element);
                t[p.prop] = element.value;
            }
            target.value = t;
            break;
        }
    }
    return offset;
}
function decode(element, buffer, offset) {
    var target = { value: void 0 };
    decodeElement(element, buffer, offset | 0, target);
    return target.value;
}
exports.decode = decode;
