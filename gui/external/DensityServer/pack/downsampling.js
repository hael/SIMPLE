"use strict";
/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
Object.defineProperty(exports, "__esModule", { value: true });
/**
 * Downsamples each slice of input data and checks if there is enough data to perform
 * higher rate downsampling.
 */
function downsampleLayer(ctx) {
    for (var i = 0, _ii = ctx.sampling.length - 1; i < _ii; i++) {
        var s = ctx.sampling[i];
        downsampleSlice(ctx, s);
        if (canDownsampleBuffer(s, false)) {
            downsampleBuffer(ctx.kernel, s, ctx.sampling[i + 1], ctx.blockSize);
        }
        else {
            break;
        }
    }
}
exports.downsampleLayer = downsampleLayer;
/**
 * When the "native" (rate = 1) sampling is finished, there might still
 * be some data left to be processed for the higher rate samplings.
 */
function finalize(ctx) {
    for (var i = 0, _ii = ctx.sampling.length - 1; i < _ii; i++) {
        var s = ctx.sampling[i];
        // skip downsampling the 1st slice because that is guaranteed to be done in "downsampleLayer"
        if (i > 0)
            downsampleSlice(ctx, s);
        // this is different from downsample layer in that it does not need 2 extra slices but just 1 is enough.
        if (canDownsampleBuffer(s, true)) {
            downsampleBuffer(ctx.kernel, s, ctx.sampling[i + 1], ctx.blockSize);
        }
        else {
            break;
        }
    }
}
exports.finalize = finalize;
/**
 * The functions downsampleH and downsampleHK both essentially do the
 * same thing: downsample along H (1st axis in axis order) and K (2nd axis in axis order) axes respectively.
 *
 * The reason there are two copies of almost the same code is performance:
 * Both functions use a different memory layout to improve cache coherency
 *  - downsampleU uses the H axis as the fastest moving one
 *  - downsampleUV uses the K axis as the fastest moving one
 */
function conv(w, c, src, b, i0, i1, i2, i3, i4) {
    return w * (c[0] * src[b + i0] + c[1] * src[b + i1] + c[2] * src[b + i2] + c[3] * src[b + i3] + c[4] * src[b + i4]);
}
/**
 * Map from L-th slice in src to an array of dimensions (srcDims[1], (srcDims[0] / 2), 1),
 * flipping the 1st and 2nd axis in the process to optimize cache coherency for downsampleUV call
 * (i.e. use (K, H, L) axis order).
 */
function downsampleH(kernel, srcDims, src, srcLOffset, buffer) {
    var target = buffer.downsampleH;
    var sizeH = srcDims[0], sizeK = srcDims[1], srcBaseOffset = srcLOffset * sizeH * sizeK;
    var targetH = Math.floor((sizeH + 1) / 2);
    var isEven = sizeH % 2 === 0;
    var w = 1.0 / kernel.coefficientSum;
    var c = kernel.coefficients;
    for (var k = 0; k < sizeK; k++) {
        var srcOffset = srcBaseOffset + k * sizeH;
        var targetOffset = k;
        target[targetOffset] = conv(w, c, src, srcOffset, 0, 0, 0, 1, 2);
        for (var h = 1; h < targetH - 1; h++) {
            srcOffset += 2;
            targetOffset += sizeK;
            target[targetOffset] = conv(w, c, src, srcOffset, -2, -1, 0, 1, 2);
        }
        srcOffset += 2;
        targetOffset += sizeK;
        if (isEven)
            target[targetOffset] = conv(w, c, src, srcOffset, -2, -1, 0, 1, 1);
        else
            target[targetOffset] = conv(w, c, src, srcOffset, -2, -1, 0, 0, 0);
    }
}
/**
 * Downsample first axis in the slice present in buffer.downsampleH
 * The result is written into the "cyclical" downsampleHk buffer
 * in the (L, H, K) axis order.
 */
function downsampleHK(kernel, dimsX, buffer) {
    var src = buffer.downsampleH, target = buffer.downsampleHK, slicesWritten = buffer.slicesWritten;
    var kernelSize = kernel.size;
    var sizeH = dimsX[0], sizeK = dimsX[1];
    var targetH = Math.floor((sizeH + 1) / 2);
    var isEven = sizeH % 2 === 0;
    var targetSliceSize = kernelSize * sizeK;
    var targetBaseOffset = slicesWritten % kernelSize;
    var w = 1.0 / kernel.coefficientSum;
    var c = kernel.coefficients;
    for (var k = 0; k < sizeK; k++) {
        var sourceOffset = k * sizeH;
        var targetOffset = targetBaseOffset + k * kernelSize;
        target[targetOffset] = conv(w, c, src, sourceOffset, 0, 0, 0, 1, 2);
        for (var h = 1; h < targetH - 1; h++) {
            sourceOffset += 2;
            targetOffset += targetSliceSize;
            target[targetOffset] = conv(w, c, src, sourceOffset, -2, -1, 0, 1, 2);
        }
        sourceOffset += 2;
        targetOffset += targetSliceSize;
        if (isEven)
            target[targetOffset] = conv(w, c, src, sourceOffset, -2, -1, 0, 1, 1);
        else
            target[targetOffset] = conv(w, c, src, sourceOffset, -2, -1, 0, 0, 0);
    }
    buffer.slicesWritten++;
}
/** Calls downsampleH and downsampleHk for each input channel separately. */
function downsampleSlice(ctx, sampling) {
    var dimsU = [sampling.sampleCount[1], Math.floor((sampling.sampleCount[0] + 1) / 2)];
    for (var i = 0, _ii = sampling.blocks.values.length; i < _ii; i++) {
        downsampleH(ctx.kernel, sampling.sampleCount, sampling.blocks.values[i], sampling.blocks.slicesWritten - 1, sampling.downsampling[i]);
        downsampleHK(ctx.kernel, dimsU, sampling.downsampling[i]);
    }
}
/** Determine if a buffer has enough data to be downsampled */
function canDownsampleBuffer(source, finishing) {
    var buffer = source.downsampling[0];
    var delta = buffer.slicesWritten - buffer.startSliceIndex;
    return (finishing && delta > 0) || (delta > 2 && (delta - 3) % 2 === 0);
}
/** Downsample data in the buffer */
function downsampleBuffer(kernel, source, target, blockSize) {
    var downsampling = source.downsampling;
    var _a = downsampling[0], slicesWritten = _a.slicesWritten, startSliceIndex = _a.startSliceIndex;
    var sizeH = target.sampleCount[0], sizeK = target.sampleCount[1], sizeHK = sizeH * sizeK;
    var kernelSize = kernel.size;
    var w = 1.0 / kernel.coefficientSum;
    var c = kernel.coefficients;
    // Indices to the 1st dimeninsion in the cyclic buffer.
    var i0 = Math.max(0, startSliceIndex - 2) % kernelSize;
    var i1 = Math.max(0, startSliceIndex - 1) % kernelSize;
    var i2 = startSliceIndex % kernelSize;
    var i3 = Math.min(slicesWritten, startSliceIndex + 1) % kernelSize;
    var i4 = Math.min(slicesWritten, startSliceIndex + 2) % kernelSize;
    var channelCount = downsampling.length;
    var valuesBaseOffset = target.blocks.slicesWritten * sizeHK;
    for (var channelIndex = 0; channelIndex < channelCount; channelIndex++) {
        var src = downsampling[channelIndex].downsampleHK;
        var values = target.blocks.values[channelIndex];
        for (var k = 0; k < sizeK; k++) {
            var valuesOffset = valuesBaseOffset + k * sizeH;
            for (var h = 0; h < sizeH; h++) {
                var sO = kernelSize * h + kernelSize * k * sizeH;
                var s = conv(w, c, src, sO, i0, i1, i2, i3, i4);
                values[valuesOffset + h] = s;
            }
        }
        // we have "consume" two layers of the buffer.
        downsampling[channelIndex].startSliceIndex += 2;
    }
    target.blocks.slicesWritten++;
}
