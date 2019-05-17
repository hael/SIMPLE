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
// import * as Data from './server/query/data-model'
// import * as Api from './server/api'
// import * as Coordinate from './server/algebra/coordinate'
var LocalApi = require("./server/local-api");
// import * as fs from 'fs'
// function wrapResponse(fn: string) {
//     const w = {
//         open(this: any) {
//             if (this.opened) return;
//             this.file = fs.openSync(fn, 'w');
//             this.opened = true;
//         },
//         writeBinary(this: any, data: Uint8Array) {
//             this.open();
//             fs.writeSync(this.file, new Buffer(data));
//             return true;
//         },
//         writeString(this: any, data: string) {
//             this.open();
//             fs.writeSync(this.file, data);
//             return true;
//         },
//         end(this: any) {
//             if (!this.opened || this.ended) return;
//             fs.close(this.file, function() { });
//             this.ended = true;
//         },
//         file: 0,
//         ended: false,
//         opened: false
//     };
//     return w;
// }
// async function query(src: string, id: string, asBinary: boolean, box: Data.QueryParamsBox, detail: number = 0, forcedSamplingLevel?: number) {
//     // const job: LocalApi.JobEntry = {
//     //     sourceFilename: `g:/test/mdb/${src}-${id}.mdb`,
//     //     sourceName: src,
//     //     sourceId: id,
//     // }
//     const params: Data.QueryParams = {
//         sourceFilename: `g:/test/mdb/${src}-${id}.mdb`,
//         sourceId: `${src}/${id}`,
//         asBinary,
//         box,
//         detail,
//         forcedSamplingLevel
//     };
//     const res = () => wrapResponse(`g:/test/${Api.getOutputFilename(src, id, params)}`);
//     await Api.queryBox(params, res)
// }
function run() {
    return __awaiter(this, void 0, void 0, function () {
        var job1;
        return __generator(this, function (_a) {
            job1 = {
                source: {
                    filename: "g:/test/mdb/emd-5725.mdb",
                    name: 'emd',
                    id: '5725',
                },
                query: {
                    kind: 'cell'
                },
                params: {
                    // detail: 4,
                    forcedSamplingLevel: 1,
                    asBinary: true
                },
                outputFolder: 'g:/test/local-test'
            };
            LocalApi.run([job1]);
            return [2 /*return*/];
        });
    });
}
run();
