#!/usr/bin/env python

import fileinput
import re

# Chromosome lengths from Ensembl GRCh37
chromosomeLengths = {
1: 249250621,
2: 243199373,
3: 198022430,
4: 191154276,
5: 180915260,
6: 171115067,
7: 159138663,
8: 146364022,
9: 141213431,
10: 135534747,
11: 135006516,
12: 133851895,
13: 115169878,
14: 107349540,
15: 102531392,
16: 90354753,
17: 81195210,
18: 78077248,
19: 59128983,
20: 63025520,
21: 48129895,
22: 51304566,
'X': 155270560,
'Y': 59373566
}

for line in fileinput.input():
  fields = line.split("\t")

  # only select sample CNVs not CNVRs
  if fields[8][:7] == "ID=essv" or fields[8][:7] == "ID=nssv" :

      # choose those with clinical interpretations
      if "clinical_int=" in fields[8] :
        clinicalInt = re.findall("(?:clinical_int=)(.+?)(?=;{1})", line)

        # convert to chromosome name
        def convertChr(x):
          return {
            'NC_000001.10': 1,
            'NC_000002.11': 2,
            'NC_000003.11': 3,
            'NC_000004.11': 4,
            'NC_000005.9': 5,
            'NC_000006.11': 6,
            'NC_000007.13': 7,
            'NC_000008.10': 8,
            'NC_000009.11': 9,
            'NC_000010.10': 10,
            'NC_000011.9': 11,
            'NC_000012.11': 12,
            'NC_000013.10': 13,
            'NC_000014.8': 14,
            'NC_000015.9': 15,
            'NC_000016.9': 16,
            'NC_000017.10': 17,
            'NC_000018.9': 18,
            'NC_000019.9': 19,
            'NC_000020.10': 20,
            'NC_000021.8': 21,
            'NC_000022.10': 22,
            'NC_000023.10': 'X',
            'NC_000024.9': 'Y'
          }.get(fields[0], None)

        # filter CNVs with lengths less than a tenth of a chromosome into into different files based on type (gain/loss) and clinical interpretation (benign/pathogenic)

        if clinicalInt and int(fields[3]) != int(fields[4]) and int(fields[4]) - int(fields[3]) < int((chromosomeLengths[convertChr(fields[0])]) / 10) and fields[0] != 'NC_000023.10' and fields[0] != 'NC_000024.9':
          if clinicalInt[0] == "Benign" or clinicalInt[0] == "no known pathogenicity":
            if fields[2] == "copy_number_loss" or fields[2] == "deletion":
              with open("files/dbVarCNLsBenign.txt", "a") as benignCNLFile:
                benignCNLFile.write("chr{}\t{}\t{}\n".format(convertChr(fields[0]), fields[3], fields[4]))
            if fields[2] == "copy_number_gain" or fields[2] == "duplication":
              with open("files/dbVarCNGsBenign.txt", "a") as benignCNGFile:
                benignCNGFile.write("chr{}\t{}\t{}\n".format(convertChr(fields[0]), fields[3], fields[4]))
          if clinicalInt[0] == "Pathogenic" or clinicalInt[0] == "pathogenic":
            if fields[2] == "copy_number_loss" or fields[2] == "deletion":
              with open("files/dbVarCNLsPathogenic.txt", "a") as pathogenicCNLFile:
                pathogenicCNLFile.write("chr{}\t{}\t{}\n".format(convertChr(fields[0]), fields[3], fields[4]))
            if fields[2] == "copy_number_gain" or fields[2] == "duplication":
              with open("files/dbVarCNGsPathogenic.txt", "a") as pathogenicCNGFile:
                pathogenicCNGFile.write("chr{}\t{}\t{}\n".format(convertChr(fields[0]), fields[3], fields[4]))