#!/usr/bin/env python

import os
import re
import json
import pydash

# Species in Ensembl to include in analysis
speciesSubset = ['bos_taurus', 'callithrix_jacchus', 'canis_familiaris',
                 'equus_caballus', 'felis_catus', 'gorilla_gorilla',
                 'macaca_mulatta', 'mus_musculus', 'oryctolagus_cuniculus',
                 'ovis_aries', 'pan_troglodytes', 'rattus_norvegicus',
                 'sus_scrofa']

# The taxonomy nodes from the mammalian root to human
taxonomyLevels = ['Homo sapiens', 'Homininae', 'Hominidae', 'Hominoidea',
                  'Catarrhini', 'Simiiformes', 'Haplorrhini', 'Primates',
                  'Euarchontoglires']

# Initialise variables
gainLossTrees = []
proteinsToGenes = {}


# Function to convert Ensembl protein IDs to Ensembl gene IDs
def convertProteinToGene(protein):
    gene = proteinsToGenes[protein]
    return gene


# Function to extract Ensembl protein IDs from gene trees in Newick format
def getProteinList(geneTreeFile):
    with open('datasets/geneTreesNewickGrch37/' + geneTreeFile + '.nh',
              'r') as geneTree:
        proteins = re.findall("(ENSP[0-9]+)", geneTree.read())
        return proteins


# Returns the intersection of two lists
intersection = lambda l, r: list(set(l).intersection(r))

# Returns merged values, removes duplicates
union = lambda l, r: list(set(l).union(r))

# Main function for analysis
def copyNumberAnalysis(name, geneList, gainLossTree, geneTree):

    # Initialise variables

    independentGenes = []
    affectedGenes = []
    affectedOrthologs = {}

    hsapProteins = re.findall(r"(ENSP\d+)", geneTree)
    ptroProteins = re.findall(r"(ENSPTRP\d+)", geneTree)
    ggorProteins = re.findall(r"(ENSGGOP\d+)", geneTree)
    mmulProteins = re.findall(r"(ENSMMUP\d+)", geneTree)
    cjacProteins = re.findall(r"(ENSCJAP\d+)", geneTree)
    rnorProteins = re.findall(r"(ENSRNOP\d+)", geneTree)
    mmusProteins = re.findall(r"(ENSMUSP\d+)", geneTree)
    ocunProteins = re.findall(r"(ENSOCUP\d+)", geneTree)
    cfamProteins = re.findall(r"(ENSCAFP\d+)", geneTree)
    fcatProteins = re.findall(r"(ENSFCAP\d+)", geneTree)
    ecabProteins = re.findall(r"(ENSECAP\d+)", geneTree)
    oariProteins = re.findall(r"(ENSOARP\d+)", geneTree)
    btauProteins = re.findall(r"(ENSBTAP\d+)", geneTree)
    sscrProteins = re.findall(r"(ENSSSCP\d+)", geneTree)

    hsapCopyNumber = len(hsapProteins)
    ptroCopyNumber = len(ptroProteins)
    ggorCopyNumber = len(ggorProteins)
    mmulCopyNumber = len(mmulProteins)
    cjacCopyNumber = len(cjacProteins)
    rnorCopyNumber = len(rnorProteins)
    mmusCopyNumber = len(mmusProteins)
    ocunCopyNumber = len(ocunProteins)
    cfamCopyNumber = len(cfamProteins)
    fcatCopyNumber = len(fcatProteins)
    ecabCopyNumber = len(ecabProteins)
    oariCopyNumber = len(oariProteins)
    btauCopyNumber = len(btauProteins)
    sscrCopyNumber = len(sscrProteins)

    ptroOrthologsFamily = []
    ggorOrthologsFamily = []
    mmulOrthologsFamily = []
    cjacOrthologsFamily = []
    rnorOrthologsFamily = []
    mmusOrthologsFamily = []
    ocunOrthologsFamily = []
    cfamOrthologsFamily = []
    fcatOrthologsFamily = []
    ecabOrthologsFamily = []
    oariOrthologsFamily = []
    btauOrthologsFamily = []
    sscrOrthologsFamily = []

    contractionPValueRegexMatches = re.findall("(?:40093268|40093269|40093275|40093291|40093295|40093297|40093299|40093301|40093302|40093303|40093305|Hsap) => \d+ \(([01]\.\d+)\) \[contraction\]", gainLossTree)

    # if contractionPValueRegexMatches:
    #     print("contractionPValueRegexMatches " + str(len(contractionPValueRegexMatches)))

    for gene in geneList:
        with open('datasets/ensemblHomologiesGrch37/' + gene + '.json',
                  'r') as geneFile:
            homologyData = json.load(geneFile)
            homologies = homologyData['data'][0]['homologies']

            numOneToOne = 0
            oneToOne = []
            oneToMany = []
            manyToMany = []
            affected = False
            ptroOrthologs = []
            ggorOrthologs = []
            mmulOrthologs = []
            cjacOrthologs = []
            rnorOrthologs = []
            mmusOrthologs = []
            ocunOrthologs = []
            cfamOrthologs = []
            fcatOrthologs = []
            ecabOrthologs = []
            oariOrthologs = []
            btauOrthologs = []
            sscrOrthologs = []

            for homolog in homologies:
                if (homolog['type'] == 'within_species_paralog' or
                        homolog['type'] == 'other_paralog' or
                        homolog['type'] == 'gene_split' or
                        homolog['type'] == 'between_species_paralog') and any(
                            homolog['taxonomy_level'] in s
                            for s in taxonomyLevels):
                    affected = True
                else:
                    if homolog['target']['species'] == 'bos_taurus':
                        btauOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'callithrix_jacchus':
                        cjacOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'canis_familiaris':
                        cfamOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'equus_caballus':
                        ecabOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'felis_catus':
                        fcatOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'gorilla_gorilla':
                        ggorOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'macaca_mulatta':
                        mmulOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'mus_musculus':
                        mmusOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'oryctolagus_cuniculus':
                        ocunOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'ovis_aries':
                        oariOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'pan_troglodytes':
                        ptroOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'rattus_norvegicus':
                        rnorOrthologs.append(homolog['target']['protein_id'])
                    elif homolog['target']['species'] == 'sus_scrofa':
                        sscrOrthologs.append(homolog['target']['protein_id'])

                    if (homolog['type'] == 'ortholog_one2one' or homolog['type'] == 'apparent_ortholog_one2one') and any(homolog['target']['species'] in s for s in speciesSubset):
                        numOneToOne += 1
                        oneToOne.append(homolog['target']['species'])
                    elif (homolog['type'] == 'ortholog_one2many' or homolog['type'] == 'apparent_ortholog_one2many') and any(homolog['target']['species'] in s for s in speciesSubset):
                        oneToMany.append(homolog['target']['species'])
                    elif (homolog['type'] == 'ortholog_many2many' or homolog['type'] == 'apparent_ortholog_many2many') and any(homolog['target']['species'] in s for s in speciesSubset):
                        manyToMany.append(homolog['target']['species'])

            if affected:
                affectedGenes.append(gene)
                orthologs = ptroOrthologs + ggorOrthologs + mmulOrthologs + cjacOrthologs + rnorOrthologs + mmusOrthologs + ocunOrthologs + cfamOrthologs + fcatOrthologs + ecabOrthologs + oariOrthologs + sscrOrthologs + btauOrthologs
                orthologs.sort()
                affectedOrthologs[gene] = orthologs
            else:
                ptroOrthologsFamily.append(ptroOrthologs)
                ggorOrthologsFamily.append(ggorOrthologs)
                mmulOrthologsFamily.append(mmulOrthologs)
                cjacOrthologsFamily.append(cjacOrthologs)
                rnorOrthologsFamily.append(rnorOrthologs)
                mmusOrthologsFamily.append(mmusOrthologs)
                ocunOrthologsFamily.append(ocunOrthologs)
                cfamOrthologsFamily.append(cfamOrthologs)
                fcatOrthologsFamily.append(fcatOrthologs)
                ecabOrthologsFamily.append(ecabOrthologs)
                oariOrthologsFamily.append(oariOrthologs)
                sscrOrthologsFamily.append(sscrOrthologs)
                btauOrthologsFamily.append(btauOrthologs)

                independentGenes.append(gene)

                oneToMany = list(set(oneToMany))
                manyToMany = list(set(manyToMany))
                numOneToMany = len(oneToMany)
                numManyToMany = len(manyToMany)
                numDup = numOneToMany + numManyToMany
                numNoOrtho = 13 - (numOneToMany + numOneToOne + numManyToMany)
                print(gene + '\t' + str(numOneToOne) + '\t' + str(numDup) + '\t' + str(numNoOrtho) + '\tindependent\t' + name)

    if affectedGenes:
        independentAndAffectedGroup = len(independentGenes) + 1
        significantContraction = False
        for pvalue in contractionPValueRegexMatches:
            if pvalue < 0.01:
                significantContraction = True
        ptroCopyNumber -= len(pydash.uniq(ptroOrthologsFamily))
        ggorCopyNumber -= len(pydash.uniq(ggorOrthologsFamily))
        mmulCopyNumber -= len(pydash.uniq(mmulOrthologsFamily))
        cjacCopyNumber -= len(pydash.uniq(cjacOrthologsFamily))
        rnorCopyNumber -= len(pydash.uniq(rnorOrthologsFamily))
        mmusCopyNumber -= len(pydash.uniq(mmusOrthologsFamily))
        ocunCopyNumber -= len(pydash.uniq(ocunOrthologsFamily))
        cfamCopyNumber -= len(pydash.uniq(cfamOrthologsFamily))
        fcatCopyNumber -= len(pydash.uniq(fcatOrthologsFamily))
        ecabCopyNumber -= len(pydash.uniq(ecabOrthologsFamily))
        oariCopyNumber -= len(pydash.uniq(oariOrthologsFamily))
        sscrCopyNumber -= len(pydash.uniq(sscrOrthologsFamily))
        btauCopyNumber -= len(pydash.uniq(btauOrthologsFamily))

        if significantContraction:
            #print('significant contraction ' + name)
            pass
        else:

            groups = []

            for key1, value1 in affectedOrthologs.iteritems():
                added = False
                for groupIndex, groupValue in enumerate(groups):
                    if len(intersection(value1, groupValue['orthologs'])):
                        groupGenes = union(groupValue['genes'],[key1])
                        groupOrthologs = union(value1, groupValue['orthologs'])

                        del groups[groupIndex]

                        for groupIndex2, groupValue2, in enumerate(groups):

                            if len(intersection(groupOrthologs, groupValue2['orthologs'])):
                                groupGenes = union(groupGenes,groupValue2['genes'])
                                groupOrthologs = union(groupOrthologs, groupValue2['orthologs'])

                                del groups[groupIndex2]

                                groups.append({'genes': groupGenes, 'orthologs': groupOrthologs})
                                added = True
                        if not added:
                            groups.append({'genes': groupGenes, 'orthologs': groupOrthologs})
                            added = True
                if not added:
                    groups.append({'genes':[key1], 'orthologs': value1})

            for group in groups:
                groupGenes = group['genes']
                groupOrthologs = group['orthologs']
                ptroCopyNumber -= len(pydash.uniq(ptroOrthologsFamily))
                ggorCopyNumber -= len(pydash.uniq(ggorOrthologsFamily))
                mmulCopyNumber -= len(pydash.uniq(mmulOrthologsFamily))
                cjacCopyNumber -= len(pydash.uniq(cjacOrthologsFamily))
                rnorCopyNumber -= len(pydash.uniq(rnorOrthologsFamily))
                mmusCopyNumber -= len(pydash.uniq(mmusOrthologsFamily))
                ocunCopyNumber -= len(pydash.uniq(ocunOrthologsFamily))
                cfamCopyNumber -= len(pydash.uniq(cfamOrthologsFamily))
                fcatCopyNumber -= len(pydash.uniq(fcatOrthologsFamily))
                ecabCopyNumber -= len(pydash.uniq(ecabOrthologsFamily))
                oariCopyNumber -= len(pydash.uniq(oariOrthologsFamily))
                sscrCopyNumber -= len(pydash.uniq(sscrOrthologsFamily))
                btauCopyNumber -= len(pydash.uniq(btauOrthologsFamily))

                numBtau = 0
                numCjac = 0
                numCfam = 0
                numEcab = 0
                numFcat = 0
                numGgor = 0
                numMmul = 0
                numMmus = 0
                numOcun = 0
                numOari = 0
                numPtro = 0
                numRnor = 0
                numSscr = 0

                for ortholog in groupOrthologs:
                    orthologMatch = re.findall(r"(\D+)", ortholog)
                    if (orthologMatch[0] == 'ENSBTAP'):
                        numBtau += 1
                    elif (orthologMatch[0] == 'ENSCJAP'):
                        numCjac += 1
                    elif (orthologMatch[0] == 'ENSCAFP'):
                        numCfam += 1
                    elif (orthologMatch[0] == 'ENSECAP'):
                        numEcab += 1
                    elif (orthologMatch[0] == 'ENSFCAP'):
                        numFcat += 1
                    elif (orthologMatch[0] == 'ENSGGOP'):
                        numGgor += 1
                    elif (orthologMatch[0] == 'ENSMMUP'):
                        numMmul += 1
                    elif (orthologMatch[0] == 'ENSMUSP'):
                        numMmus += 1
                    elif (orthologMatch[0] == 'ENSOCUP'):
                        numOcun += 1
                    elif (orthologMatch[0] == 'ENSOARP'):
                        numOari += 1
                    elif (orthologMatch[0] == 'ENSPTRP'):
                        numPtro += 1
                    elif (orthologMatch[0] == 'ENSRNOP'):
                        numRnor += 1
                    elif (orthologMatch[0] == 'ENSSSCP'):
                        numSscr += 1
                numOneToOne = 0
                numDup = 0
                numNoOrtho = 0

                for count in [numBtau, numCjac, numCfam, numEcab, numFcat, numGgor, numMmul, numMmus, numOcun, numOari, numPtro, numRnor, numSscr]:
                    if count == 1:
                        numOneToOne += 1
                    elif count > 1:
                        numDup += 1
                    else:
                        numNoOrtho += 1

                for groupGene in groupGenes:
                    print(groupGene + '\t' + str(numOneToOne) + '\t' + str(numDup) + '\t' + str(numNoOrtho) + '\taffected\t' + name)


with open('datasets/ensGrch37ProteinGeneList.txt', 'r') as ensGrch37ProteinGeneList:
    lines = ensGrch37ProteinGeneList.read().splitlines()
    for line in lines:
        fields = line.split('\t')
        proteinsToGenes[fields[0]] = fields[1]

for fileName in os.listdir("datasets/gainLossTreesGrch37/"):
    if fileName.startswith("ENSGT") and fileName.endswith(".txt"):
        gainLossTrees.append(fileName)

for gainLossTreeFile in gainLossTrees:
    name = gainLossTreeFile[:-4]

    with open('datasets/geneTreesNewickGrch37/' + name + '.nh', 'r') as geneTree, open('datasets/gainLossTreesGrch37/' + gainLossTreeFile, 'r') as gainLossTree:

        gainLossTree = gainLossTree.read()
        geneTree = geneTree.read()

        # Check that there is non-zero copy number at mammalian root as only want to include these gene families that existed at mammalian root
        copyNumberAtMammalianRoot = re.findall("(40093267) => [1-9]\d*", gainLossTree)
        if not copyNumberAtMammalianRoot:
            # Copy number at mammalian root is zero so exclude these genes from the analysis and take note of them
            with open('geneFamiliesWithZeroAtMammalianRoot.txt', 'a') as geneFamiliesWithZeroAtMammalianRootFile:
                geneFamiliesWithZeroAtMammalianRootFile.write(name + '\n')
            proteins = getProteinList(name)
            genes = map(lambda x: convertProteinToGene(x), proteins)
            with open('genesFromFamiliesWithZeroAtMammalianRoot.txt', 'a') as genesFromFamiliesWithZeroAtMammalianRootPython:
                for gene in genes:
                    genesFromFamiliesWithZeroAtMammalianRootPython.write(gene + '\n')
        else:
            # Copy number at mammalian root is not zero so proceed with analysis

            proteins = getProteinList(name)
            genes = map(lambda x: convertProteinToGene(x), proteins)
            copyNumberAnalysis(name, genes, gainLossTree, geneTree)
