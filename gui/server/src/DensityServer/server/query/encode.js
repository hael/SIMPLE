"use strict";
/*
 * Copyright (c) 2016 - now, David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
Object.defineProperty(exports, "__esModule", { value: true });
var CIF = require("../../lib/cif-tools");
var version_1 = require("../version");
var DataFormat = require("../../common/data-format");
function encode(query, output) {
    var w = query.params.asBinary
        ? new CIF.Binary.Writer("DensityServer " + version_1.default)
        : new CIF.Text.Writer();
    write(w, query);
    w.encode();
    w.flush(output);
}
exports.default = encode;
var E = CIF.Binary.Encoder;
function string(name, str, isSpecified) {
    if (isSpecified) {
        return { name: name, string: str, presence: function (data, i) { return isSpecified(data, i) ? 0 /* Present */ : 1 /* NotSpecified */; } };
    }
    return { name: name, string: str };
}
function int32(name, num) {
    return { name: name, string: function (data, i) { return '' + num(data, i); }, number: num, typedArray: Int32Array, encoder: E.by(E.byteArray) };
}
function float64(name, num, precision) {
    if (precision === void 0) { precision = 1000000; }
    return { name: name, string: function (data, i) { return '' + Math.round(precision * num(data, i)) / precision; }, number: num, typedArray: Float64Array, encoder: E.by(E.byteArray) };
}
var _volume_data_3d_info_fields = [
    string('name', function (ctx) { return ctx.header.channels[ctx.channelIndex]; }),
    int32('axis_order[0]', function (ctx) { return ctx.header.axisOrder[0]; }),
    int32('axis_order[1]', function (ctx) { return ctx.header.axisOrder[1]; }),
    int32('axis_order[2]', function (ctx) { return ctx.header.axisOrder[2]; }),
    float64('origin[0]', function (ctx) { return ctx.grid.origin[0]; }),
    float64('origin[1]', function (ctx) { return ctx.grid.origin[1]; }),
    float64('origin[2]', function (ctx) { return ctx.grid.origin[2]; }),
    float64('dimensions[0]', function (ctx) { return ctx.grid.dimensions[0]; }),
    float64('dimensions[1]', function (ctx) { return ctx.grid.dimensions[1]; }),
    float64('dimensions[2]', function (ctx) { return ctx.grid.dimensions[2]; }),
    int32('sample_rate', function (ctx) { return ctx.sampleRate; }),
    int32('sample_count[0]', function (ctx) { return ctx.grid.sampleCount[0]; }),
    int32('sample_count[1]', function (ctx) { return ctx.grid.sampleCount[1]; }),
    int32('sample_count[2]', function (ctx) { return ctx.grid.sampleCount[2]; }),
    int32('spacegroup_number', function (ctx) { return ctx.header.spacegroup.number; }),
    float64('spacegroup_cell_size[0]', function (ctx) { return ctx.header.spacegroup.size[0]; }, 1000),
    float64('spacegroup_cell_size[1]', function (ctx) { return ctx.header.spacegroup.size[1]; }, 1000),
    float64('spacegroup_cell_size[2]', function (ctx) { return ctx.header.spacegroup.size[2]; }, 1000),
    float64('spacegroup_cell_angles[0]', function (ctx) { return ctx.header.spacegroup.angles[0]; }, 1000),
    float64('spacegroup_cell_angles[1]', function (ctx) { return ctx.header.spacegroup.angles[1]; }, 1000),
    float64('spacegroup_cell_angles[2]', function (ctx) { return ctx.header.spacegroup.angles[2]; }, 1000),
    float64('mean_source', function (ctx) { return ctx.globalValuesInfo.mean; }),
    float64('mean_sampled', function (ctx) { return ctx.sampledValuesInfo.mean; }),
    float64('sigma_source', function (ctx) { return ctx.globalValuesInfo.sigma; }),
    float64('sigma_sampled', function (ctx) { return ctx.sampledValuesInfo.sigma; }),
    float64('min_source', function (ctx) { return ctx.globalValuesInfo.min; }),
    float64('min_sampled', function (ctx) { return ctx.sampledValuesInfo.min; }),
    float64('max_source', function (ctx) { return ctx.globalValuesInfo.max; }),
    float64('max_sampled', function (ctx) { return ctx.sampledValuesInfo.max; })
];
function _volume_data_3d_info(result) {
    var ctx = {
        header: result.query.data.header,
        channelIndex: result.channelIndex,
        grid: result.query.samplingInfo.gridDomain,
        sampleRate: result.query.samplingInfo.sampling.rate,
        globalValuesInfo: result.query.data.header.sampling[0].valuesInfo[result.channelIndex],
        sampledValuesInfo: result.query.data.header.sampling[result.query.samplingInfo.sampling.index].valuesInfo[result.channelIndex]
    };
    return {
        data: ctx,
        count: 1,
        desc: {
            name: '_volume_data_3d_info',
            fields: _volume_data_3d_info_fields
        }
    };
}
function _volume_data_3d_str(ctx, i) {
    return '' + Math.round(1000000 * ctx[i]) / 1000000;
}
function _volume_data_3d_number(ctx, i) {
    return ctx[i];
}
function _volume_data_3d(ctx) {
    var data = ctx.query.values[ctx.channelIndex];
    var encoder;
    var typedArray;
    if (ctx.query.data.header.valueType === DataFormat.ValueType.Float32 || ctx.query.data.header.valueType === DataFormat.ValueType.Int16) {
        var min = void 0, max = void 0;
        min = data[0], max = data[0];
        for (var i = 0, n = data.length; i < n; i++) {
            var v = data[i];
            if (v < min)
                min = v;
            else if (v > max)
                max = v;
        }
        typedArray = Float32Array;
        // encode into 255 steps and store each value in 1 byte.
        encoder = E.by(E.intervalQuantizaiton(min, max, 255, Uint8Array)).and(E.byteArray);
    }
    else {
        typedArray = Int8Array;
        // just encode the bytes
        encoder = E.by(E.byteArray);
    }
    var fields = [{
            name: 'values',
            string: _volume_data_3d_str,
            number: _volume_data_3d_number,
            typedArray: typedArray,
            encoder: encoder
        }];
    return {
        data: data,
        count: data.length,
        desc: {
            name: '_volume_data_3d',
            fields: fields
        }
    };
}
function pickQueryBoxDimension(ctx, e, d) {
    var box = ctx.params.box;
    switch (box.kind) {
        case 'Cartesian':
        case 'Fractional':
            return "" + Math.round(1000000 * box[e][d]) / 1000000;
        default: return '';
    }
}
function queryBoxDimension(e, d) {
    return string("query_box_" + e + "[" + d + "]", function (ctx) { return pickQueryBoxDimension(ctx, e, d); }, function (ctx) { return ctx.params.box.kind !== 'Cell'; });
}
var _density_server_result_fields = [
    string('server_version', function (ctx) { return version_1.default; }),
    string('datetime_utc', function (ctx) { return new Date().toISOString().replace(/T/, ' ').replace(/\..+/, ''); }),
    string('guid', function (ctx) { return ctx.guid; }),
    string('is_empty', function (ctx) { return ctx.kind === 'Empty' || ctx.kind === 'Error' ? 'yes' : 'no'; }),
    string('has_error', function (ctx) { return ctx.kind === 'Error' ? 'yes' : 'no'; }),
    string('error', function (ctx) { return ctx.kind === 'Error' ? ctx.message : ''; }, function (ctx) { return ctx.kind === 'Error'; }),
    string('query_source_id', function (ctx) { return ctx.params.sourceId; }),
    string('query_type', function (ctx) { return 'box'; }),
    string('query_box_type', function (ctx) { return ctx.params.box.kind.toLowerCase(); }),
    queryBoxDimension('a', 0),
    queryBoxDimension('a', 1),
    queryBoxDimension('a', 2),
    queryBoxDimension('b', 0),
    queryBoxDimension('b', 1),
    queryBoxDimension('b', 2)
];
function _density_server_result(ctx) {
    return {
        data: ctx,
        count: 1,
        desc: {
            name: '_density_server_result',
            fields: _density_server_result_fields
        }
    };
}
function write(writer, query) {
    writer.startDataBlock('SERVER');
    writer.writeCategory(_density_server_result, [query]);
    switch (query.kind) {
        case 'Data':
    }
    if (query.kind === 'Data') {
        var header = query.data.header;
        for (var i = 0; i < header.channels.length; i++) {
            writer.startDataBlock(header.channels[i]);
            var ctx = [{ query: query, channelIndex: i }];
            writer.writeCategory(_volume_data_3d_info, ctx);
            writer.writeCategory(_volume_data_3d, ctx);
        }
    }
}
