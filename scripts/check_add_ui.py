#!/usr/bin/env python3
import re
from pathlib import Path
root=Path('src/main/ui')
pat=re.compile(r"add_ui_program\(\s*'([^']+)'\s*,\s*([A-Za-z0-9_]+)\s*,")
miss=[]
for f in root.rglob('*.f90'):
    text=f.read_text()
    for m in pat.finditer(text):
        s=m.group(1)
        inst=m.group(2)
        if s!=inst:
            # record line number
            lineno=text[:m.start()].count('\n')+1
            miss.append((str(f),lineno,s,inst))
if not miss:
    print('No mismatches found')
else:
    for f,ln,s,inst in miss:
        print(f+":"+str(ln)+": '{}' vs {}".format(s,inst))
    print('\nTotal mismatches:', len(miss))
