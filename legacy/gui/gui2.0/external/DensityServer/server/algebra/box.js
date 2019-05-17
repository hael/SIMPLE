"use strict";
/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
Object.defineProperty(exports, "__esModule", { value: true });
var Coords = require("./coordinate");
///////////////////////////////////////////
// CONVERSIONS
///////////////////////////////////////////
function cartesianToFractional(box, spacegroup) {
    var l = box.a, r = box.b;
    var corners = [
        [l[0], l[1], l[2]],
        [r[0], l[1], l[2]],
        [l[0], r[1], l[2]],
        [l[0], l[1], r[2]],
        [r[0], r[1], l[2]],
        [r[0], l[1], r[2]],
        [l[0], r[1], r[2]],
        [r[0], r[1], r[2]],
    ].map(function (c) { return Coords.cartesianToFractional(Coords.cartesian(c[0], c[1], c[2]), spacegroup); });
    return bounding(corners);
}
exports.cartesianToFractional = cartesianToFractional;
function fractionalToGrid(box, domain) {
    return { a: Coords.fractionalToGrid(box.a, domain, 'bottom'), b: Coords.fractionalToGrid(box.b, domain, 'top') };
}
exports.fractionalToGrid = fractionalToGrid;
function gridToFractional(box) {
    return { a: Coords.gridToFractional(box.a), b: Coords.gridToFractional(box.b) };
}
exports.gridToFractional = gridToFractional;
function fractionalBoxReorderAxes(box, axisOrder) {
    var a = box.a, b = box.b;
    return {
        a: Coords.withCoord(a, a[axisOrder[0]], a[axisOrder[1]], a[axisOrder[2]]),
        b: Coords.withCoord(b, b[axisOrder[0]], b[axisOrder[1]], b[axisOrder[2]])
    };
}
exports.fractionalBoxReorderAxes = fractionalBoxReorderAxes;
function expandGridBox(box, by) {
    var a = box.a, b = box.b;
    return {
        a: Coords.withCoord(a, a[0] - by, a[1] - by, a[2] - by),
        b: Coords.withCoord(b, b[0] + by, b[1] + by, b[2] + by)
    };
}
exports.expandGridBox = expandGridBox;
///////////////////////////////////////////
// MISC
///////////////////////////////////////////
function shift(box, offset) {
    return { a: Coords.add(box.a, offset), b: Coords.add(box.b, offset) };
}
exports.shift = shift;
function clampGridToSamples(box) {
    return { a: Coords.clampGridToSamples(box.a), b: Coords.clampGridToSamples(box.b) };
}
exports.clampGridToSamples = clampGridToSamples;
function fractionalToDomain(box, kind, delta) {
    var ds = Coords.fractional(box.b[0] - box.a[0], box.b[1] - box.a[1], box.b[2] - box.a[2]);
    return Coords.domain(kind, {
        delta: delta,
        origin: box.a,
        dimensions: ds,
        sampleCount: Coords.sampleCounts(ds, delta)
    });
}
exports.fractionalToDomain = fractionalToDomain;
function fractionalFromBlock(block) {
    var domain = block.domain;
    var a = Coords.gridToFractional(block);
    var b = Coords.add(a, domain.delta);
    for (var i = 0; i < 3; i++) {
        b[i] = Math.min(b[i], domain.origin[i] + domain.dimensions[i]);
    }
    return { a: a, b: b };
}
exports.fractionalFromBlock = fractionalFromBlock;
function bounding(xs) {
    var a = Coords.clone(xs[0]);
    var b = Coords.clone(xs[0]);
    for (var _i = 0, xs_1 = xs; _i < xs_1.length; _i++) {
        var x = xs_1[_i];
        for (var i = 0; i < 3; i++) {
            a[i] = Math.min(a[i], x[i]);
            b[i] = Math.max(b[i], x[i]);
        }
    }
    return { a: a, b: b };
}
exports.bounding = bounding;
function areIntersecting(box1, box2) {
    for (var i = 0; i < 3; i++) {
        var x = box1.a[i], y = box1.b[i];
        var u = box2.a[i], v = box2.b[i];
        if (x > v || y < u)
            return false;
    }
    return true;
}
exports.areIntersecting = areIntersecting;
function intersect(box1, box2) {
    var a = Coords.clone(box1.a);
    var b = Coords.clone(box1.a);
    for (var i = 0; i < 3; i++) {
        var x = box1.a[i], y = box1.b[i];
        var u = box2.a[i], v = box2.b[i];
        if (x > v || y < u)
            return void 0;
        a[i] = Math.max(x, u);
        b[i] = Math.min(y, v);
    }
    return { a: a, b: b };
}
exports.intersect = intersect;
function dimensions(box) {
    return [box.b[0] - box.a[0], box.b[1] - box.a[1], box.b[2] - box.a[2]];
}
exports.dimensions = dimensions;
function volume(box) {
    return (box.b[0] - box.a[0]) * (box.b[1] - box.a[1]) * (box.b[2] - box.a[2]);
}
exports.volume = volume;
