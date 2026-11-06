#!/usr/bin/env python3
"""probe_entity_users.py -- find the users of Fortran entities with the compiler, not with grep.

apply renames each entity in its module family (the module and its non-tester submodules) and in
the files given with --accept, after backing every edited file up. A build (BUILD_TESTS=ON) then
fails at each reference outside those files; classify sorts the errors of the build log into
production and test-only users; restore puts the files back. A build that succeeds proves that
no other file references the entities; a failing build lists some users, not necessarily all
(make stops scheduling after the first failure): accept or fix them and build again.

usage:
  scripts/probe_entity_users.py apply ROOT ENTITY... [--accept FILE]... [--backup DIR]
  scripts/probe_entity_users.py classify LOG
  scripts/probe_entity_users.py restore ROOT [--backup DIR]
ENTITY: <module>:<name> (a procedure, constant, variable, type or interface of the module) or
        <type>%<binding> (a type-bound procedure or generic binding)
Renamed names end in _probe. Warnings name reach the compiler cannot see (generics, overrides of a
parent's binding, bind(c)).
"""

import os
import re
import shutil
import sys

SUFFIX = "_probe"
DEFAULT_BACKUP = os.path.join(os.environ.get("TMPDIR", "/tmp"), "simple_probe_backup")
# test-only code, as src/CMakeLists.txt excludes it without BUILD_TESTS
TEST_ONLY = [re.compile(p) for p in (r"_tester\.[Ff]90$", r"/main/commanders/test/",
                                      r"/main/exec/simple_test_exec[^/]*\.[Ff]90$",
                                      r"/main/apis/simple_test_exec_api\.[Ff]90$")]


def is_test_only(path):
    return any(p.search(path) for p in TEST_ONLY)


def sources(root):
    return sorted(os.path.relpath(os.path.join(d, n), root) for top in ("src", "production")
                  for d, _, names in os.walk(os.path.join(root, top)) for n in names if n.endswith((".f90", ".F90")))


def read(root, rel):
    with open(os.path.join(root, rel), errors="replace") as f:
        return f.read()


def module_family(root, files, module):
    """the module's file and its submodules (any depth), testers excluded"""
    mod_re = re.compile(r"(?im)^\s*module\s+%s\s*$" % re.escape(module))
    sub_re = re.compile(r"(?im)^\s*submodule\s*\(\s*(\w+)\s*(?::\s*(\w+))?\s*\)\s*(\w+)")
    home = [f for f in files if mod_re.search(read(root, f))]
    if len(home) != 1:
        sys.exit("module %s: %d defining files" % (module, len(home)))
    parent_of, file_of = {}, {}
    for f in files:
        for m in sub_re.finditer(read(root, f)):
            parent_of[m.group(3).lower()] = (m.group(2) or m.group(1)).lower()
            file_of[m.group(3).lower()] = f
    family = list(home)
    for sub, f in file_of.items():
        anc, seen = parent_of[sub], set()
        while anc in parent_of and anc not in seen:
            seen.add(anc)
            anc = parent_of[anc]
        if anc == module.lower() and not f.endswith(("_tester.f90", "_tester.F90")):
            family.append(f)
    return family


def type_home(root, files, tname):
    t_re = re.compile(r"(?im)^\s*type\s*(?:,[^:\n]*)?::\s*%s\b" % re.escape(tname))
    hits = [f for f in files if not is_test_only(f) and t_re.search(read(root, f))]
    if len(hits) != 1:
        sys.exit("type %s: %d defining files %s" % (tname, len(hits), hits))
    module = re.search(r"(?im)^\s*module\s+(\w+)\s*$", read(root, hits[0])).group(1)
    return module


def type_block(text, tname):
    m = re.search(r"(?ims)^\s*type\s*(,[^:\n]*)?::\s*%s\b.*?^\s*end\s*type" % re.escape(tname), text)
    return m


def warn_binding(root, files, text, tname, b):
    warnings = []
    blk = type_block(text, tname)
    body = blk.group(0)
    for g in re.finditer(r"(?im)^\s*generic[^:\n]*::\s*([\w()=.+\-*/ ]+?)\s*=>\s*(.*)$", body):
        if re.search(r"(?i)\b%s\b" % re.escape(b), g.group(2)):
            warnings.append("%%%s is a specific of generic %s: calls through the generic reach it; "
                            "probe %s%%%s too" % (b, g.group(1).strip(), tname, g.group(1).strip()))
    ext = re.search(r"(?i)extends\s*\(\s*(\w+)", blk.group(1) or "")
    if ext:
        parent = ext.group(1)
        for f in files:
            pt = type_block(read(root, f), parent)
            if pt and re.search(r"(?im)^\s*(procedure|generic)[^:\n]*::.*\b%s\b" % re.escape(b), pt.group(0)):
                warnings.append("%%%s overrides %s%%%s (%s): calls through class(%s) dispatch here and "
                                "the probe cannot see them" % (b, parent, b, f, parent))
    return warnings


def apply(root, entities, accepted, backup):
    if os.path.isdir(backup) and os.listdir(backup):
        sys.exit("%s holds a probe in progress: restore it first" % backup)
    os.makedirs(backup, exist_ok=True)
    files = sources(root)
    edits, warnings = {}, []

    def text_of(f):
        if f not in edits:
            edits[f] = read(root, f)
        return edits[f]

    for ent in entities:
        if "%" in ent:
            tname, b = ent.split("%", 1)
            module = type_home(root, files, tname)
            family = module_family(root, files, module)
            home = family[0]
            text = text_of(home)
            warnings += warn_binding(root, files, text, tname, b)
            blk = type_block(text, tname)
            body = blk.group(0)
            line_re = re.compile(r"(?im)^(\s*(?:procedure|generic)[^:\n]*::\s*)(%s\b)(\s*=>\s*\w+)?" % re.escape(b))
            m = line_re.search(body)
            if not m:
                sys.exit("%s: no binding %s in type %s" % (home, b, tname))
            target = m.group(3) or " => " + b
            new_body = body[:m.start()] + m.group(1) + b + SUFFIX + target + body[m.end():]
            # a generic of the type that lists the binding follows the new name
            new_body = re.sub(r"(?im)^(\s*generic[^:\n]*::[^=\n]*=>)(.*)$",
                              lambda g: g.group(1) + re.sub(r"(?i)\b%s\b" % re.escape(b), b + SUFFIX, g.group(2)),
                              new_body)
            edits[home] = text[:blk.start()] + new_body + text[blk.end():]
            for f in family + accepted:
                edits[f] = re.sub(r"(?i)%%(\s*)%s\b" % re.escape(b), r"%\g<1>" + b + SUFFIX, text_of(f))
        else:
            module, name = ent.split(":", 1)
            family = module_family(root, files, module)
            for f in family:
                t = text_of(f)
                if re.search(r"(?im)^.*\bbind\s*\(\s*c\b.*\b%s\b" % re.escape(name), t):
                    warnings.append("%s:%s is bind(c): C callers are invisible to the probe" % (module, name))
                for blk in re.finditer(r"(?ims)^\s*interface\s+(\w+|operator\s*\([^)]*\)|assignment\s*\(=\))\s*$"
                                       r"(.*?)^\s*end\s*interface", t):
                    if re.search(r"(?i)\b%s\b" % re.escape(name), blk.group(2)):
                        warnings.append("%s:%s is a specific of generic %s: calls through the generic "
                                        "reach it; probe the generic too" % (module, name, blk.group(1)))
            for f in family + accepted:
                edits[f] = re.sub(r"(?i)\b%s\b" % re.escape(name), name + SUFFIX, text_of(f))
    for f, t in edits.items():
        dst = os.path.join(backup, f.replace("/", "__"))
        shutil.copy2(os.path.join(root, f), dst)
        with open(os.path.join(root, f), "w") as fh:
            fh.write(t)
    print("probed %d entities in %d files (backups in %s)" % (len(entities), len(edits), backup))
    for w in sorted(set(warnings)):
        print("WARNING:", w)


def classify(log):
    loc = re.compile(r"^(\S+\.[Ff]90):(\d+):\d+:\s*$")
    seen, pending = set(), None
    with open(log, errors="replace") as fh:
        for line in fh:
            m = loc.match(line.strip())
            if m:
                pending = (m.group(1), m.group(2))
                continue
            if pending and line.startswith("Error:"):
                key = pending
                if key not in seen:
                    seen.add(key)
                    kind = "TEST      " if is_test_only(key[0]) else "PRODUCTION"
                    print("%s %s:%s  %s" % (kind, key[0], key[1], line.strip()[7:]))
                pending = None
    if not seen:
        print("no compile errors in", log)


def restore(root, backup):
    if not os.path.isdir(backup):
        sys.exit("no backups in " + backup)
    n = 0
    for name in os.listdir(backup):
        dst = os.path.join(root, name.replace("__", "/"))
        shutil.copyfile(os.path.join(backup, name), dst)
        # a new time stamp: a tree that compiled the probed file must compile it again
        os.utime(dst, None)
        n += 1
    shutil.rmtree(backup)
    print("restored %d files" % n)


def main():
    args = sys.argv[1:]
    if not args or args[0] not in ("apply", "classify", "restore"):
        sys.exit(__doc__)
    backup = DEFAULT_BACKUP
    if "--backup" in args:
        i = args.index("--backup")
        backup = os.path.abspath(args[i + 1])
        del args[i:i + 2]
    accepted = []
    while "--accept" in args:
        i = args.index("--accept")
        accepted.append(args[i + 1])
        del args[i:i + 2]
    if args[0] == "apply":
        apply(os.path.abspath(args[1]), args[2:], accepted, backup)
    elif args[0] == "classify":
        classify(args[1])
    else:
        restore(os.path.abspath(args[1]), backup)


if __name__ == "__main__":
    main()
