#! /usr/bin/env python

import sys

pop = sys.argv[2]

print("chrm,winMid,tP,nSites,pi,pop")
with open(sys.argv[1]) as f:
  next(f)
  for line in f:
    line = line.strip()
    fields = line.split()
    tP, sites = float(fields[3]), int(fields[-1])
    chrm, mid = fields[1], fields[2]
    if sites != 0:
      print(f"{chrm},{mid},{tP},{sites},{tP/sites},{pop}")
    
    
