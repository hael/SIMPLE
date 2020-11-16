"use strict";
/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
Object.defineProperty(exports, "__esModule", { value: true });
var LA = require("./linear");
/** Constructs spacegroup skew matrix from supplied info */
function spacegroup(info) {
    var cellAngles = info.angles, cellSize = info.size;
    var alpha = (Math.PI / 180.0) * cellAngles[0];
    var beta = (Math.PI / 180.0) * cellAngles[1];
    var gamma = (Math.PI / 180.0) * cellAngles[2];
    var xScale = cellSize[0], yScale = cellSize[1], zScale = cellSize[2];
    var z1 = Math.cos(beta);
    var z2 = (Math.cos(alpha) - Math.cos(beta) * Math.cos(gamma)) / Math.sin(gamma);
    var z3 = Math.sqrt(1.0 - z1 * z1 - z2 * z2);
    var x = [xScale, 0.0, 0.0];
    var y = [Math.cos(gamma) * yScale, Math.sin(gamma) * yScale, 0.0];
    var z = [z1 * zScale, z2 * zScale, z3 * zScale];
    var fromFrac = LA.Matrix4.ofRows([
        [x[0], y[0], z[0], 0],
        [0, y[1], z[1], 0],
        [0, 0, z[2], 0],
        [0, 0, 0, 1.0]
    ]);
    var toFrac = LA.Matrix4.invert(LA.Matrix4.empty(), fromFrac);
    return { angles: info.angles, size: info.size, number: info.number, toFrac: toFrac, fromFrac: fromFrac };
}
exports.spacegroup = spacegroup;
///////////////////////////////////////////
// CONSTRUCTORS
///////////////////////////////////////////
function domain(kind, info) {
    var sc = info.sampleCount;
    return {
        kind: kind,
        delta: info.delta,
        dimensions: info.dimensions,
        origin: info.origin,
        sampleCount: info.sampleCount,
        sampleVolume: sc[0] * sc[1] * sc[2]
    };
}
exports.domain = domain;
function cartesian(x, y, z) {
    return { 0: x, 1: y, 2: z, kind: 0 /* Cartesian */ };
}
exports.cartesian = cartesian;
function fractional(x, y, z) {
    return { 0: x, 1: y, 2: z, kind: 1 /* Fractional */ };
}
exports.fractional = fractional;
function grid(domain, x, y, z) {
    return { 0: x, 1: y, 2: z, kind: 2 /* Grid */, domain: domain };
}
exports.grid = grid;
function withCoord(a, x, y, z) {
    switch (a.kind) {
        case 0 /* Cartesian */: return cartesian(x, y, z);
        case 1 /* Fractional */: return fractional(x, y, z);
        case 2 /* Grid */: return grid(a.domain, x, y, z);
    }
}
exports.withCoord = withCoord;
function clone(a) {
    return withCoord(a, a[0], a[1], a[2]);
}
exports.clone = clone;
///////////////////////////////////////////
// CONVERSIONS
///////////////////////////////////////////
function cartesianToFractional(a, spacegroup) {
    var coord = Helpers.transform(a, spacegroup.toFrac);
    return fractional(coord[0], coord[1], coord[2]);
}
exports.cartesianToFractional = cartesianToFractional;
function fractionalToGrid(a, domain, snap) {
    var origin = domain.origin, delta = domain.delta;
    var coord = grid(domain, 0.1, 0.1, 0.1);
    for (var i = 0; i < 3; i++) {
        coord[i] = Helpers.snap((a[i] - origin[i]) / delta[i], snap);
    }
    return coord;
}
exports.fractionalToGrid = fractionalToGrid;
function gridToFractional(a) {
    var _a = a.domain, origin = _a.origin, delta = _a.delta;
    var coord = fractional(0.1, 0.1, 0.1);
    for (var i = 0; i < 3; i++) {
        coord[i] = a[i] * delta[i] + origin[i];
    }
    return coord;
}
exports.gridToFractional = gridToFractional;
///////////////////////////////////////////
// MISC
///////////////////////////////////////////
function clampGridToSamples(a) {
    var sampleCount = a.domain.sampleCount;
    var coord = withCoord(a, 0, 0, 0);
    for (var i = 0; i < 3; i++) {
        if (a[i] < 0)
            coord[i] = 0;
        else if (a[i] > sampleCount[i])
            coord[i] = sampleCount[i];
        else
            coord[i] = a[i];
    }
    return coord;
}
exports.clampGridToSamples = clampGridToSamples;
function add(a, b) {
    return withCoord(a, a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
exports.add = add;
function sub(a, b) {
    return withCoord(a, a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
exports.sub = sub;
function invert(a) {
    return withCoord(a, -a[0], -a[1], -a[2]);
}
exports.invert = invert;
/** Maps each grid point to a unique integer */
function linearGridIndex(a) {
    var samples = a.domain.sampleCount;
    return a[0] + samples[0] * (a[1] + a[2] * samples[1]);
}
exports.linearGridIndex = linearGridIndex;
function gridMetrics(dimensions) {
    return {
        sizeX: dimensions[0],
        sizeXY: dimensions[0] * dimensions[1],
        sizeXYZ: dimensions[0] * dimensions[1] * dimensions[2]
    };
}
exports.gridMetrics = gridMetrics;
function sampleCounts(dimensions, delta) {
    return [
        Helpers.snap(dimensions[0] / delta[0], 'top'),
        Helpers.snap(dimensions[1] / delta[1], 'top'),
        Helpers.snap(dimensions[2] / delta[2], 'top')
    ];
}
exports.sampleCounts = sampleCounts;
// to prevent floating point rounding errors
function round(v) {
    return Math.round(10000000 * v) / 10000000;
}
exports.round = round;
var Helpers;
(function (Helpers) {
    var applyTransform = LA.Matrix4.transformVector3;
    function transform(x, matrix) {
        return applyTransform([0.1, 0.1, 0.1], x, matrix);
    }
    Helpers.transform = transform;
    function snap(v, to) {
        switch (to) {
            case 'bottom': return Math.floor(round(v)) | 0;
            case 'top': return Math.ceil(round(v)) | 0;
        }
    }
    Helpers.snap = snap;
})(Helpers || (Helpers = {}));
