#!/usr/bin/env python

import fileinput

lastChr = '';
lastStart = 0;
lastEnd = 0;
lastCoverage = 0;

secondLastChr = '';
secondLastStart = 0;
secondLastEnd = 0;
secondLastCoverage = 0;

for line in fileinput.input():
  fields = line.split("\t")

  if (lastCoverage > int(fields[3]) and lastCoverage > secondLastCoverage and lastChr == fields[0] and secondLastChr == fields[0]) or (secondLastEnd != lastStart and lastEnd != int(fields[1])):
    print(lastChr + "\t" + str(lastStart) + "\t" + str(lastEnd))

  secondLastChr = lastChr
  secondLastStart = lastStart
  secondLastEnd = lastEnd
  secondLastCoverage = lastCoverage

  lastChr = fields[0]
  lastStart = int(fields[1])
  lastEnd = int(fields[2])
  lastCoverage = int(fields[3])