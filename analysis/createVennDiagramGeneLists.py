#!/usr/bin/env python

import pydash

benignGain = []
benignLoss = []
pathoGain = []
pathoLoss = []

lists = {
  'BGBLPG': [],
  'BGBLPL': [],
  'BGPGPL': [],
  'BLPGPL': [],
  'BGBL': [],
  'BGPG': [],
  'BGPL': [],
  'BLPG': [],
  'BLPL': [],
  'PGPL': [],
  'BG': [],
  'BL': [],
  'PG': [],
  'PL': [],
  'BGBLPGPL': []
}

with open('files/gene lists/dbVarCNGsBenignPCGenes.txt', 'r') as dbVarCNGsBenignPCGenes:
    benignGain = dbVarCNGsBenignPCGenes.read().splitlines()

with open('files/gene lists/dbVarCNLsBenignPCGenes.txt', 'r') as dbVarCNLsBenignPCGenes:
    benignLoss = dbVarCNLsBenignPCGenes.read().splitlines()

with open('files/gene lists/dbVarCNGPeaksPathogenicPCGenes.txt', 'r') as dbVarCNGPeaksPathogenicPCGenes:
    pathoGain = dbVarCNGPeaksPathogenicPCGenes.read().splitlines()

with open('files/gene lists/dbVarCNLPeaksPathogenicPCGenes.txt', 'r') as dbVarCNLPeaksPathogenicPCGenes:
    pathoLoss = dbVarCNLPeaksPathogenicPCGenes.read().splitlines()

lists['BGBLPGPL'] = pydash.intersection(benignGain, benignLoss, pathoGain, pathoLoss)

_BGBLPG = pydash.intersection(benignGain, benignLoss, pathoGain)
lists['BGBLPG'] = pydash.difference(_BGBLPG, pathoLoss)

_BGBLPL = pydash.intersection(benignGain, benignLoss, pathoLoss)
lists['BGBLPL'] = pydash.difference(_BGBLPL, pathoGain)

_BGPGPL = pydash.intersection(benignGain, pathoGain, pathoLoss)
lists['BGPGPL'] = pydash.difference(_BGPGPL, benignLoss)

_BLPGPL = pydash.intersection(benignLoss, pathoGain, pathoLoss)
lists['BLPGPL'] = pydash.difference(_BLPGPL, benignGain)

_BGBL = pydash.intersection(benignGain, benignLoss)
lists['BGBL'] = pydash.difference(_BGBL, pathoGain, pathoLoss)

_BGPG = pydash.intersection(benignGain, pathoGain)
lists['BGPG'] = pydash.difference(_BGPG, benignLoss, pathoLoss)

_BGPL = pydash.intersection(benignGain, pathoLoss)
lists['BGPL'] = pydash.difference(_BGPL, benignLoss, pathoGain)

_BLPG = pydash.intersection(benignLoss, pathoGain)
lists['BLPG'] = pydash.difference(_BLPG, benignGain, pathoLoss)

_BLPL = pydash.intersection(benignLoss, pathoLoss)
lists['BLPL'] = pydash.difference(_BLPL, benignGain, pathoGain)

_PGPL = pydash.intersection(pathoGain, pathoLoss)
lists['PGPL'] = pydash.difference(_PGPL, benignGain, benignLoss)

lists['BG'] = pydash.difference(benignGain, benignLoss, pathoGain, pathoLoss)
lists['BL'] = pydash.difference(benignLoss, benignGain, pathoGain, pathoLoss)
lists['PG'] = pydash.difference(pathoGain, benignGain, benignLoss, pathoLoss)
lists['PL'] = pydash.difference(pathoLoss, benignGain, benignLoss, pathoGain)

for key in lists:
	with open('files/gene lists/vennDiagram/'+key+'.txt', 'a') as output:
		output.write("\n".join(lists[key])+"\n")
