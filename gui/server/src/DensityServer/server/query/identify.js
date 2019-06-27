"use strict";
/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
Object.defineProperty(exports, "__esModule", { value: true });
var Coords = require("../algebra/coordinate");
var Box = require("../algebra/box");
var collections_1 = require("../utils/collections");
/** Find a list of unique blocks+offsets that overlap with the query region. */
function findUniqueBlocks(data, sampling, queryBox) {
    var translations = data.header.spacegroup.isPeriodic
        // find all query box translations that overlap with the unit cell.
        ? findDataOverlapTranslationList(queryBox, sampling.dataDomain)
        // no translations
        : [Coords.fractional(0, 0, 0)];
    var blocks = collections_1.FastMap.create();
    for (var _i = 0, translations_1 = translations; _i < translations_1.length; _i++) {
        var t = translations_1[_i];
        findUniqueBlocksOffset(data, sampling, queryBox, t, blocks);
    }
    var blockList = blocks.forEach(function (b, _, ctx) { ctx.push(b); }, []);
    // sort the data so that the first coodinate changes the fastest 
    // this is because that's how the data is laid out in the underlaying 
    // data format and reading the data 'in order' makes it faster.
    blockList.sort(function (a, b) {
        var x = a.coord, y = b.coord;
        for (var i = 2; i >= 0; i--) {
            if (x[i] !== y[i])
                return x[i] - y[i];
        }
        return 0;
    });
    return blockList;
}
exports.default = findUniqueBlocks;
/**
 * Find the integer interval [x, y] so that for all k \in [x, y]
 * [a + k, b + k] intersects with (u, v)
 */
function overlapMultiplierRange(a, b, u, v) {
    var x = Math.ceil(u - b) | 0, y = Math.floor(v - a) | 0;
    if (Coords.round(b + x) <= Coords.round(u))
        x++;
    if (Coords.round(a + y) >= Coords.round(v))
        y--;
    if (x > y)
        return void 0;
    return [x, y];
}
/**
 * Finds that list of "unit" offsets (in fractional space) so that
 * shift(box, offset) has non-empty interaction with the region
 * described in the give domain.
 */
function findDataOverlapTranslationList(box, domain) {
    var ranges = [];
    var translations = [];
    var origin = domain.origin, dimensions = domain.dimensions;
    for (var i = 0; i < 3; i++) {
        var range = overlapMultiplierRange(box.a[i], box.b[i], origin[i], origin[i] + dimensions[i]);
        if (!range)
            return translations;
        ranges[i] = range;
    }
    var u = ranges[0], v = ranges[1], w = ranges[2];
    for (var k = w[0]; k <= w[1]; k++) {
        for (var j = v[0]; j <= v[1]; j++) {
            for (var i = u[0]; i <= u[1]; i++) {
                translations.push(Coords.fractional(i, j, k));
            }
        }
    }
    return translations;
}
function addUniqueBlock(blocks, coord, offset) {
    var hash = Coords.linearGridIndex(coord);
    if (blocks.has(hash)) {
        var entry = blocks.get(hash);
        entry.offsets.push(offset);
    }
    else {
        blocks.set(hash, { coord: coord, offsets: [offset] });
    }
}
function findUniqueBlocksOffset(data, sampling, queryBox, offset, blocks) {
    var shifted = Box.shift(queryBox, offset);
    var intersection = Box.intersect(shifted, data.dataBox);
    // Intersection can be empty in the case of "aperiodic spacegroups"
    if (!intersection)
        return;
    var blockDomain = sampling.blockDomain;
    // this gets the "3d range" of block indices that contain data that overlaps 
    // with the query region.
    //
    // Clamping the data makes sure we avoid silly rounding errors (hopefully :))
    var _a = Box.clampGridToSamples(Box.fractionalToGrid(intersection, blockDomain)), min = _a.a, max = _a.b;
    for (var i = min[0]; i < max[0]; i++) {
        for (var j = min[1]; j < max[1]; j++) {
            for (var k = min[2]; k < max[2]; k++) {
                addUniqueBlock(blocks, Coords.grid(blockDomain, i, j, k), offset);
            }
        }
    }
}
