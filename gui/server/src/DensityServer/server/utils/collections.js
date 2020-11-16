"use strict";
/*
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */
Object.defineProperty(exports, "__esModule", { value: true });
function createMapObject() {
    var map = Object.create(null);
    // to cause deoptimization as we don't want to create hidden classes
    map['__'] = void 0;
    delete map['__'];
    return map;
}
var FastMap;
(function (FastMap) {
    function forEach(data, f, ctx) {
        var hasOwn = Object.prototype.hasOwnProperty;
        for (var _i = 0, _a = Object.keys(data); _i < _a.length; _i++) {
            var key = _a[_i];
            if (!hasOwn.call(data, key))
                continue;
            var v = data[key];
            if (v === void 0)
                continue;
            f(v, key, ctx);
        }
        return ctx;
    }
    var fastMap = {
        set: function (key, v) {
            if (this.data[key] === void 0 && v !== void 0) {
                this.size++;
            }
            this.data[key] = v;
        },
        get: function (key) {
            return this.data[key];
        },
        delete: function (key) {
            if (this.data[key] === void 0)
                return false;
            delete this.data[key];
            this.size--;
            return true;
        },
        has: function (key) {
            return this.data[key] !== void 0;
        },
        clear: function () {
            this.data = createMapObject();
            this.size = 0;
        },
        forEach: function (f, ctx) {
            return forEach(this.data, f, ctx !== void 0 ? ctx : void 0);
        }
    };
    /**
     * Creates an empty map.
     */
    function create() {
        var ret = Object.create(fastMap);
        ret.data = createMapObject();
        ret.size = 0;
        return ret;
    }
    FastMap.create = create;
    /**
     * Create a map from an array of the form [[key, value], ...]
     */
    function ofArray(data) {
        var ret = create();
        for (var _i = 0, data_1 = data; _i < data_1.length; _i++) {
            var xs = data_1[_i];
            ret.set(xs[0], xs[1]);
        }
        return ret;
    }
    FastMap.ofArray = ofArray;
    /**
     * Create a map from an object of the form { key: value, ... }
     */
    function ofObject(data) {
        var ret = create();
        var hasOwn = Object.prototype.hasOwnProperty;
        for (var _i = 0, _a = Object.keys(data); _i < _a.length; _i++) {
            var key = _a[_i];
            if (!hasOwn.call(data, key))
                continue;
            var v = data[key];
            ret.set(key, v);
        }
        return ret;
    }
    FastMap.ofObject = ofObject;
})(FastMap = exports.FastMap || (exports.FastMap = {}));
var FastSet;
(function (FastSet) {
    function forEach(data, f, ctx) {
        var hasOwn = Object.prototype.hasOwnProperty;
        for (var _i = 0, _a = Object.keys(data); _i < _a.length; _i++) {
            var p = _a[_i];
            if (!hasOwn.call(data, p) || data[p] !== null)
                continue;
            f(p, ctx);
        }
        return ctx;
    }
    /**
     * Uses null for present values.
     */
    var fastSet = {
        add: function (key) {
            if (this.data[key] === null)
                return false;
            this.data[key] = null;
            this.size++;
            return true;
        },
        delete: function (key) {
            if (this.data[key] !== null)
                return false;
            delete this.data[key];
            this.size--;
            return true;
        },
        has: function (key) {
            return this.data[key] === null;
        },
        clear: function () {
            this.data = createMapObject();
            this.size = 0;
        },
        forEach: function (f, ctx) {
            return forEach(this.data, f, ctx !== void 0 ? ctx : void 0);
        }
    };
    /**
     * Create an empty set.
     */
    function create() {
        var ret = Object.create(fastSet);
        ret.data = createMapObject();
        ret.size = 0;
        return ret;
    }
    FastSet.create = create;
    /**
     * Create a set of an "array like" sequence.
     */
    function ofArray(xs) {
        var ret = create();
        for (var i = 0, l = xs.length; i < l; i++) {
            ret.add(xs[i]);
        }
        return ret;
    }
    FastSet.ofArray = ofArray;
})(FastSet = exports.FastSet || (exports.FastSet = {}));
