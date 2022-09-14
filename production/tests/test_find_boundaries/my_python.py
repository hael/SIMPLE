# Python program to test calling python from Fortran

print("Entering Python")

import sys
import numpy as np
import ruptures as rpt


if len(sys.argv) != 5:
    raise Exception("Execute as python my_python.py filename outfile min_size pen")
fn = sys.argv[1]
if not fn.endswith(".csv"):
    raise Exception("Input correlation time series must be a .csv file")
fn_out = sys.argv[2]
if not fn_out.endswith(".txt"):
    raise Exception("Output boundary file must be a .txt file")
min_size = int(sys.argv[3])
if min_size <= 0:
    raise Exception("Argument min_size must be a positive integer")
pen = float(sys.argv[4])
if pen <= 0:
    raise Exception("Penalty must be positive")
f = open(fn_out, 'w')
f.write("file: " + fn + "\n")
f.write("min_size: " + str(min_size) + "\n")
f.write("pen: " + str(pen) + "\n")
f.close()

ptcls, corrs = np.transpose(np.genfromtxt(fn, delimiter=',', skip_header = 1))
rbf = rpt.KernelCPD(kernel='rbf', min_size=min_size).fit(corrs)
bounds = rbf.predict(pen=pen)[:-1]
print("Writing boundaries to " + fn_out)
f = open(fn_out, 'a')
for bound in bounds:
    f.write(str(ptcls[bound - 1]) + "\n")
f.close()

print("Exiting Python")