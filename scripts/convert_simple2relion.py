import os
import sys

nargs = len(sys.argv)
print nargs
if nargs != 3 :
    print 'Missing arguments. First argument is the simple document, second is the stack_name'
    sys.exit(-1)

# ARGUMENTS
simple_doc = sys.argv[1]   # simple doc
stack_name = sys.argv[2]   # name of stack

# hack to manually add cs,kv,frac (also see below)
# ctf_str  = ' cs={:3.1f}'.format(2.7)
# ctf_str += ' kv={:5.1f}'.format(300.)
# ctf_str += ' fraca={:5.3f}'.format(.07)
# end hack

# PARSE DOC
if not os.path.isfile( simple_doc ):
    print 'Could not find file:', simple_doc
    sys.exit(-1)
# trim and optionally hack
lines = open(simple_doc,'r').readlines()
for i in range(len(lines)):
    line = lines[i].strip()
    # hack to manually add cs,kv,frac
    # line += ctf_str
    # end hack
    lines[i] = line

# keys equivalence
key_dict = { 
             'e1'    :'AngleRot',    # to comment to parse ctf parms only
             'e2'    :'AngleTilt',   # to comment to parse ctf parms only
             'e3'    :'AnglePsi',    # to comment to parse ctf parms only
             'x'     :'OriginX',     # to comment to parse ctf parms only
             'y'     :'OriginY',     # to comment to parse ctf parms only
             'dfx'   :'DefocusU',
             'dfy'   :'DefocusV',
             'angast':'DefocusAngle',
             'kv'    :'Voltage',
             'cs'    :'SphericalAberration',
             'fraca' :'AmplitudeContrast',
            }

def write_header():
    print 'data_images\n\nloop_'
    line_here = lines[0]
    fields = line_here.strip().split()
    cnt    = 1
    for field in fields:
        key = field.split('=')[0]
        if key in key_dict.keys():
            print '_rln'+key_dict[ key ]+' #'+str(cnt)
            cnt += 1
    print '_rlnImageName #'+str(cnt)   # not sure this is necessary 
    return None

def write_body():
    cnt = 0
    for line in lines:
        cnt += 1
        fields = line.split()
        strout = " "
        for field in fields:
            vals = field.split('=')
            key  = vals[0]
            if key in key_dict.keys():
                arg = vals[1]
                val = float(arg)
                if key=='dfx' or key=='dfy':
                    val = val * 10000.
                strout += ("%.6f " % val)
        strout += ("%d@" % cnt) +stack_name
        print strout
    return None

# MAIN
write_header()
write_body()
