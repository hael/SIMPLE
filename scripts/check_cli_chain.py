#!/usr/bin/env python3
"""Static check of the command-line chain: UI inputs, accepted keys, bindings, private programs.

Run by scripts/run_fast_gate.sh before CTest:

    python3 scripts/check_cli_chain.py [REPO_ROOT]

A command-line key passes through four places: its declaration in type(parameters)
(src/main/params/simple_parameters.f90), the accepted-key list generated from the
declarations (scripts/simple_args_generator.pl), its binding onto the field
(src/main/params/simple_parameters_parse.f90) and the program descriptors that offer
it (src/main/ui). A key missing from one of them fails without a word or at a user's
first try, so this check fails the gate on each of these:

- a UI input, or a key of a private program's table, that is no accepted key
  (giving it stops the program);
- an accepted key that parameters never binds (it is accepted and ignored), unless
  parameters reads it from the command line itself or it is one of KEYS_READ_OUTSIDE,
  which a program reads from its command line directly;
- a key bound twice;
- a yes/no key (declared with the choices "(yes|no)") bound without add_yesno, which
  refuses any other value, add_yesno on another key, or a yes/no (binary) UI input whose
  declaration names no (yes|no) choices;
- a key of the stream launcher written by hand (print_stream_ui_json in
  src/main/ui/simple_ui.f90) that is no accepted key;
- a program of simple_private_exec that has neither a public UI descriptor nor an
  entry in src/main/ui/simple_private_prgs.f90 (its command line is never checked),
  a table entry no dispatch case runs, or a name in the table twice.

A field that parameters derives or sets itself is no key: its declaration's comment
says "not a command-line key" and the generator leaves it out. The accepted keys are
the generator's own output (run here into a temporary directory), so the check and
the build agree by construction.

Python 3.6, standard library only (the gate runs python3).
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

# Accepted keys that a program reads from its command line itself, not through
# parameters (the field of that name is derived): each must be read from a cmdline
# somewhere outside src/main/params, which the check verifies.
KEYS_READ_OUTSIDE = {
    'nrots': 'make_pickrefs reads the number of in-plane rotations of its references '
             '(simple_pick_strategy.f90); the field is the search grid the builder derives',
}

PARAMS = os.path.join('src', 'main', 'params', 'simple_parameters.f90')
PARSE = os.path.join('src', 'main', 'params', 'simple_parameters_parse.f90')
PARAMS_DIR = os.path.join('src', 'main', 'params')
GENERATOR = os.path.join('scripts', 'simple_args_generator.pl')
UI_DIR = os.path.join('src', 'main', 'ui')
PRIVATE_PRGS = os.path.join('src', 'main', 'ui', 'simple_private_prgs.f90')
PRIVATE_DRIVER = os.path.join('production', 'simple_private_exec_driver.f90')
CMDLINE = os.path.join('src', 'main', 'simple_cmdline.f90')
SELF = os.path.join('scripts', 'check_cli_chain.py')
UI_JSON = os.path.join('src', 'main', 'ui', 'simple_ui.f90')


def read(root, rel):
    with open(os.path.join(root, rel), encoding='utf-8', errors='replace') as handle:
        return handle.read()


def strip_comments(text):
    """Fortran source without comments (a '!' outside a string starts one)."""
    out = []
    for line in text.split('\n'):
        quote = None
        for i, char in enumerate(line):
            if quote:
                if char == quote:
                    quote = None
            elif char in '\'"':
                quote = char
            elif char == '!':
                line = line[:i]
                break
        out.append(line)
    return '\n'.join(out)


def join_continuations(text):
    """Free-form continuation lines joined into one."""
    return re.sub(r'&[ \t]*\n[ \t]*&?', ' ', text)


def fortran_files(root, rel_dir):
    for dirpath, _, files in os.walk(os.path.join(root, rel_dir)):
        for name in sorted(files):
            if name.endswith('.f90'):
                yield os.path.relpath(os.path.join(dirpath, name), root)


def accepted_keys(root):
    """The generator's accepted keys: (all keys, keys from declarations)."""
    tmp = tempfile.mkdtemp(prefix='check_cli_chain_')
    try:
        env = dict(os.environ)
        env['SIMPLE_PATH'] = tmp
        proc = subprocess.run(['perl', os.path.join(root, GENERATOR), tmp],
                              cwd=os.path.join(root, PARAMS_DIR), env=env,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True)
        if proc.returncode:
            sys.stderr.write(proc.stdout + proc.stderr)
            raise SystemExit('check_cli_chain: the accepted-key generator failed')
        with open(os.path.join(tmp, 'simple_varlist.txt')) as handle:
            declared = [line.strip() for line in handle if line.strip()]
        with open(os.path.join(tmp, 'simple_args.f90')) as handle:
            accepted = set(re.findall(r"self%args\(\d+\) = '([^']+)'", handle.read()))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return accepted, set(declared) - {'vol%d' % i for i in range(1, 21)}


def declared_choices(root):
    """The choices "(a|b|...)" a character component of type(parameters) declares in its comment."""
    text = read(root, PARAMS)
    start = re.search(r'^\s*type\s*(::)?\s*parameters\b', text, re.M | re.I)
    end = re.search(r'^\s*contains\b', text[start.end():], re.M)
    choices = {}
    for line in text[start.end():start.end() + end.start()].split('\n'):
        match = re.match(r"\s*character\(len=\w+\)\s*::\s*(\w+)\s*=?[^!]*!(.*)$", line)
        if not match:
            continue
        found = re.search(r'\(([^()]*\|[^()]*)\)', match.group(2))
        if found:
            choices[match.group(1)] = set(found.group(1).split('|'))
    return choices


def bindings(root):
    """(kind, key, line) of every registry binding in the parse submodule."""
    out = []
    for lineno, line in enumerate(read(root, PARSE).split('\n'), 1):
        code = strip_comments(line)
        for kind, key in re.findall(r"reg%add_(\w+)\(\s*'([^']+)'", code):
            out.append((kind, key, lineno))
    return out


CLINE_READ = re.compile(r"%(?:defined|get_carg|get_iarg|get_rarg)\(\s*'([^']+)'")


def read_by_parameters(root):
    """Keys parameters reads from the command line itself (outside the registry)."""
    keys = set()
    for rel in fortran_files(root, PARAMS_DIR):
        if rel.endswith('_tester.f90'):
            continue
        keys.update(CLINE_READ.findall(strip_comments(read(root, rel))))
    return keys


def read_outside_parameters(root, key):
    pattern = re.compile(r"%(?:defined|get_carg|get_iarg|get_rarg)\(\s*'" + re.escape(key) + "'")
    for rel_dir in ('src', 'production'):
        for rel in fortran_files(root, rel_dir):
            if rel.startswith(PARAMS_DIR) or rel.endswith('_tester.f90'):
                continue
            if pattern.search(strip_comments(read(root, rel))):
                return True
    return False


def ui_inputs(root):
    """(key, file) of every UI input and every key of the private programs' table."""
    texts = {rel: join_continuations(strip_comments(read(root, rel))) for rel in fortran_files(root, UI_DIR)
             if not rel.endswith('_tester.f90')}
    param_keys = {}
    for text in texts.values():
        for var, key in re.findall(r"call\s+(\w+)%set_param\(\s*'([^']+)'", text):
            param_keys[var.lower()] = key
    inputs = []
    unresolved = []
    binary = set()
    for text in texts.values():
        binary.update(re.findall(r"add_input\(\s*\w+\s*,\s*'([^']+)'\s*,\s*'binary'", text))
        binary.update(re.findall(r"set_param\(\s*'([^']+)'\s*,\s*'binary'", text))
    for rel, text in texts.items():
        for arg in re.findall(r"%add_input\(\s*\w+\s*,\s*('[^']+'|\w+)", text):
            if arg.startswith("'"):
                inputs.append((arg.strip("'"), rel))
            elif arg.lower() in param_keys:
                inputs.append((param_keys[arg.lower()], rel))
            else:
                unresolved.append((arg, rel))
        for key in re.findall(r"(?:%push_(?:req|opt)_key|\bpush_(?:req|opt))\(\s*'([^']+)'", text):
            inputs.append((key, rel))
    return inputs, unresolved, binary


def stream_launcher_keys(root):
    """(key, line) of the keys the hand-written stream launcher names."""
    text = read(root, UI_JSON)
    start = text.index('subroutine print_stream_ui_json')
    end = text.index('end subroutine print_stream_ui_json')
    first = text[:start].count('\n') + 1
    keys = []
    for offset, line in enumerate(strip_comments(text[start:end]).split('\n')):
        found = re.findall(r"json%add\(\s*input\s*,\s*'key'\s*,\s*'([^']+)'", line)
        found += re.findall(r"json%add\(\s*process_inputs\s*,\s*''\s*,\s*'([^'=]+)'\s*\)", line)
        found += re.findall(r"json%add\(\s*process_inputs\s*,\s*''\s*,\s*'([^'=]+)=", line)
        keys.extend((key.strip(), first + offset) for key in found)
    return keys


def public_programs(root):
    names = set()
    for rel in fortran_files(root, UI_DIR):
        names.update(re.findall(r"add_ui_program\(\s*'([^']+)'", strip_comments(read(root, rel))))
    return names


def private_dispatch(root):
    """Program names of dispatch_private_prg's select case."""
    text = strip_comments(read(root, PRIVATE_DRIVER))
    body = re.search(r'subroutine dispatch_private_prg(.*?)end subroutine dispatch_private_prg', text, re.S)
    names = []
    for case in re.findall(r"^\s*case\s*\(([^)]*)\)", body.group(1), re.M):
        names.extend(re.findall(r"'([^']+)'", case))
    return names


def private_table(root):
    return re.findall(r"(?:call\s+add_private_prg|%set_name)\(\s*'([^']+)'", strip_comments(read(root, PRIVATE_PRGS)))


def keyless_utilities(root):
    """The utility programs cmdline parses without keys (no table entry needed)."""
    text = strip_comments(read(root, CMDLINE))
    body = re.search(r'subroutine parse_private\b(.*?)end subroutine parse_private\b', text, re.S)
    case = re.search(r"case\s*\(([^)]*)\)\s*\n\s*return", body.group(1))
    return set(re.findall(r"'([^']+)'", case.group(1))) if case else set()


def main():
    root = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..'))
    errors = []

    accepted, declared = accepted_keys(root)
    bound = bindings(root)
    bound_keys = {key for _, key, _ in bound}

    # a key bound twice
    seen = {}
    for kind, key, lineno in bound:
        if key in seen:
            errors.append('%s:%d: %s is bound twice (first at line %d)' % (PARSE, lineno, key, seen[key]))
        else:
            seen[key] = lineno

    # every accepted key is bound
    by_parameters = read_by_parameters(root)
    for key in sorted(declared):
        if key in bound_keys or key in by_parameters or key in KEYS_READ_OUTSIDE:
            continue
        errors.append('%s: %s is an accepted key that parameters never binds: bind it in %s, or, if parameters '
                      'derives the field, say "not a command-line key" in its declaration' % (PARAMS, key, PARSE))
    for key, why in sorted(KEYS_READ_OUTSIDE.items()):
        if key not in accepted:
            errors.append('%s: %s in KEYS_READ_OUTSIDE is no accepted key' % (SELF, key))
        elif key in bound_keys:
            errors.append('%s: %s in KEYS_READ_OUTSIDE is bound: remove it from the list' % (SELF, key))
        elif not read_outside_parameters(root, key):
            errors.append('%s: %s in KEYS_READ_OUTSIDE is read from no command line (%s)' % (SELF, key, why))
    for kind, key, lineno in bound:
        if key not in accepted:
            errors.append('%s:%d: %s is bound but no accepted key' % (PARSE, lineno, key))

    # yes/no keys are bound with add_yesno
    choices = declared_choices(root)
    yes_no = {key for key, options in choices.items() if options == {'yes', 'no'}}
    for kind, key, lineno in bound:
        if key in yes_no and kind != 'yesno':
            errors.append('%s:%d: %s is a yes/no key bound with add_%s: bind it with add_yesno' % (PARSE, lineno, key, kind))
        if kind == 'yesno' and key not in yes_no:
            errors.append('%s:%d: %s is bound with add_yesno but not declared with the choices (yes|no)' % (PARSE, lineno, key))

    # every UI input is an accepted key
    inputs, unresolved, binary = ui_inputs(root)
    for key in sorted(binary & declared):
        if not {'yes', 'no'} <= choices.get(key, set()):
            errors.append('%s: %s is a yes/no (binary) UI input, but its declaration names no choices (yes|no): '
                          'add them, and bind it with add_yesno' % (PARAMS, key))
    for arg, rel in unresolved:
        errors.append('%s: add_input names %s, which no set_param in %s defines' % (rel, arg, UI_DIR))
    reported = set()
    for key, rel in inputs:
        if key not in accepted and (key, rel) not in reported:
            reported.add((key, rel))
            errors.append('%s: input %s is no accepted key: giving it stops the program' % (rel, key))

    # the hand-written stream launcher names accepted keys only
    launcher = stream_launcher_keys(root)
    for key, lineno in launcher:
        if key not in accepted:
            errors.append('%s:%d: the stream launcher names %s, which is no accepted key' % (UI_JSON, lineno, key))

    # every private program has its keys checked
    public = public_programs(root)
    table = private_table(root)
    dispatched = private_dispatch(root)
    keyless = keyless_utilities(root)
    for name in dispatched:
        if name not in public and name not in table and name not in keyless:
            errors.append('%s: private program %s has neither a public UI descriptor nor an entry in %s: '
                          'its command line is never checked' % (PRIVATE_DRIVER, name, PRIVATE_PRGS))
    for name in sorted(set(table)):
        if table.count(name) > 1:
            errors.append('%s: %s has %d entries' % (PRIVATE_PRGS, name, table.count(name)))
        if name not in dispatched:
            errors.append('%s: %s has an entry but no case in dispatch_private_prg' % (PRIVATE_PRGS, name))
        if name in public:
            errors.append('%s: %s has an entry but its public descriptor already checks its keys' % (PRIVATE_PRGS, name))

    if errors:
        sys.stderr.write('check_cli_chain: %d problem(s) in the command-line chain:\n' % len(errors))
        for error in errors:
            sys.stderr.write('  ' + error + '\n')
        return 1
    print('Command-line chain consistent: %d accepted keys, %d bindings, %d UI inputs, %d stream launcher keys, '
          '%d private programs.' % (len(accepted), len(bound), len(inputs), len(launcher), len(dispatched)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
