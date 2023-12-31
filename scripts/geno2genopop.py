#! /usr/bin/env python

import sys, gzip

variants = {
  "A" : '01',
  "C" : '02',
  "G" : '03',
  "T" : '04',
  "N" : '00'
}

main = []
header = ['locus']

with open(sys.argv[1]) as f:
  for line in f:
    line = line.strip()
    header.append(line) 

RG = []
Pecos = []
with open(sys.argv[3]) as f:
  for line in f:
    line = line.strip()
    fields = line.split()
    if fields[1] == 'RioGrande':
      RG.append(fields[0])
    elif fields[1] == 'Pecos':
      Pecos.append(fields[0])


main.append(header)

with gzip.open(sys.argv[2], 'rt') as f:
  for line in f:
    line = line.strip()
    fields = line.split()
    locus = fields[0] + '_' + fields[1]
    genos = fields[2:]
    newGeno = [locus]
    for i in genos:
      tmp = variants[i[0]] + variants[i[1]]
      newGeno.append(tmp)
    if newGeno.count('0000') / len(newGeno) <= 0.1: 
      main.append(newGeno) 

transposed = list(map(list, zip(*main)))
print("#customGenePopfromANGSD")
print(','.join(transposed[0][1:]))
print("pop")
for i in transposed[1:]:
  if i[0] in RG:
    print(f"{i[0]},\t{'\t'.join(i[1:])}")

print("pop")
for i in transposed[1:]:
  if i[0] in Pecos:
    print(f"{i[0]},\t{'\t'.join(i[1:])}") 


