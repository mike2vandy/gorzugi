#! /usr/bin/env python

import sys, gzip

sites = {}
with open(sys.argv[1]) as f:
  for line in f:
    line = line.strip()
    sites[line] = ''

with gzip.open(sys.argv[2],'rt') as f:
  for line in f:
    line = line.strip()
    if line.startswith('marker'):
      print(line)
    else:
      key = line.split()[0].replace('_',':')
      if key in sites:
        print(line)
