#! /usr/bin/env python

import sys, glob

samples = {}
with open(sys.argv[1]) as f:
  for line in f:
    line = line.strip()
    fields = line.split()
    samples[fields[0]] = fields[1]

print("sample,pop,Hobs")
for i in glob.glob("*.sfs"):
  sname = i.split('.')[0]
  with open(i) as f:
    for line in f:
      line = line.strip()
      fields = line.split()
      hom, het = float(fields[0]), float(fields[1])
      print(f"{sname},{samples[sname]},{het/(hom+het)}")
