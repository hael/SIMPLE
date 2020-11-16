"use strict";
/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
Object.defineProperty(exports, "__esModule", { value: true });
function formatTime(t) {
    if (isNaN(t))
        return 'n/a';
    var h = Math.floor(t / (60 * 60 * 1000)), m = Math.floor(t / (60 * 1000) % 60), s = Math.floor(t / 1000 % 60), ms = Math.floor(t % 1000).toString();
    while (ms.length < 3)
        ms = '0' + ms;
    if (h > 0)
        return h + "h" + m + "m" + s + "." + ms + "s";
    if (m > 0)
        return m + "m" + s + "." + ms + "s";
    if (s > 0)
        return s + "." + ms + "s";
    return t.toFixed(0) + "ms";
}
exports.formatTime = formatTime;
function logPlain(tag, msg) {
    console.log("[" + tag + "] " + msg);
}
exports.logPlain = logPlain;
function log(guid, tag, msg) {
    console.log("[" + guid + "][" + tag + "] " + msg);
}
exports.log = log;
function errorPlain(ctx, e) {
    console.error("[Error] (" + ctx + ") " + e);
    if (e.stack)
        console.error(e.stack);
}
exports.errorPlain = errorPlain;
function error(guid, e) {
    console.error("[" + guid + "][Error] " + e);
    if (e.stack)
        console.error(e.stack);
}
exports.error = error;
