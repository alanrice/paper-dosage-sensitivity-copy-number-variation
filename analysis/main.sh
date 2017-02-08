#!/bin/sh

# ---Some housekeeping before running main analysis---

set -euo pipefail
IFS=$'\n\t'

version=1.0.0

# Check dependencies are installed
hash bedtools 2>/dev/null || { echo >&2 "Bedtools required but it's not installed. Install and try again."; exit 1; }
hash bedops 2>/dev/null || { echo >&2 "Bedop required but it's not installed. Install and try again."; exit 1; }
hash python 2>/dev/null || { echo >&2 "Python required but it's not installed. Install and try again."; exit 1; }

#---End housekeeping---

echo "Start of analysis"
echo "Version $version"
echo "-------------"

mkdir -p files
mkdir -p files/gene\ lists
mkdir -p files/gene\ lists/vennDiagram
mkdir -p files/speciesCopyNumberAnalysis
mkdir -p files/haploinsufficiencyScores
mkdir -p files/maxExpression
mkdir -p files/medianExpression

echo "Unzipping datasets..."
file1md5=`md5 -q datasets/GRCh37.2013_10_31.remap.all.germline.gvf.gz`
if [ $file1md5 != "764edc5dbf0d0a67e116aa4d60db7dd3" ]
then
  echo >&2 "GRCh37.2013_10_31.remap.all.germline.gvf.gz seems to be corrupted. Try download again."; exit 1;
fi

file2md5=`md5 -q datasets/GRCh37.2013_10_31.submitted.all.germline.gvf.gz`
if [ $file2md5 != "4ccd5634f7fa9aca564a303596aa6886" ]
then
  echo >&2 "GRCh37.2013_10_31.submitted.all.germline.gvf.gz seems to be corrupted. Try download again."; exit 1;
fi

file3md5=`md5 -q datasets/ensemblHomologiesGrch37.zip`
if [ $file3md5 != "3ed39862115238bf4a66cc4992805c7c" ]
then
  echo >&2 "ensemblHomologiesGrch37.zip seems to be corrupted. Try download again."; exit 1;
fi

file4md5=`md5 -q datasets/geneTreesNewickGrch37.zip`
if [ $file4md5 != "1052e44cd6b835417420316c258b7e39" ]
then
  echo >&2 "geneTreesNewickGrch37.zip seems to be corrupted. Try download again."; exit 1;
fi

file5md5=`md5 -q datasets/gainLossTreesGrch37.zip`
if [ $file5md5 != "a71bf15cc2831722ba0340da22032f4d" ]
then
  echo >&2 "gainLossTreesGrch37.zip seems to be corrupted. Try download again."; exit 1;
fi

echo "Decompressing dbVar CNV data..."
gunzip --keep datasets/GRCh37.2013_10_31.remap.all.germline.gvf.gz
gunzip --keep datasets/GRCh37.2013_10_31.submitted.all.germline.gvf.gz
unzip -q datasets/ensemblHomologiesGrch37.zip -d datasets/
unzip -q datasets/geneTreesNewickGrch37.zip -d datasets/
unzip -q datasets/gainLossTreesGrch37.zip -d datasets/

echo "Filtering dbVar CNV file..."
cat datasets/GRCh37.2013_10_31.submitted.all.germline.gvf datasets/GRCh37.2013_10_31.remap.all.germline.gvf | grep "\t" | python dbVarFilter.py

cat files/dbVarCNGsBenign.txt files/dbVarCNLsBenign.txt > files/dbVarCNVsBenign.txt
cat files/dbVarCNGsPathogenic.txt files/dbVarCNLsPathogenic.txt > files/dbVarCNVsPathogenic.txt

echo "Sorting CNV files..."
cat files/dbVarCNGsBenign.txt | sort-bed - > files/dbVarCNGsBenign.bed
cat files/dbVarCNLsBenign.txt | sort-bed - > files/dbVarCNLsBenign.bed
cat files/dbVarCNVsBenign.txt | sort-bed - > files/dbVarCNVsBenign.bed
cat files/dbVarCNGsPathogenic.txt | sort-bed - > files/dbVarCNGsPathogenic.bed
cat files/dbVarCNLsPathogenic.txt | sort-bed - > files/dbVarCNLsPathogenic.bed
cat files/dbVarCNVsPathogenic.txt | sort-bed - > files/dbVarCNVsPathogenic.bed

wc -l files/dbVarCNGsBenign.bed
wc -l files/dbVarCNLsBenign.bed
wc -l files/dbVarCNVsBenign.bed
wc -l files/dbVarCNGsPathogenic.bed
wc -l files/dbVarCNLsPathogenic.bed
wc -l files/dbVarCNVsPathogenic.bed

echo "Creating unique CNV files..."
cat files/dbVarCNGsBenign.txt | sort-bed - | uniq > files/dbVarCNGsBenign.unique.bed
cat files/dbVarCNLsBenign.txt | sort-bed - | uniq > files/dbVarCNLsBenign.unique.bed
cat files/dbVarCNVsBenign.txt | sort-bed - | uniq > files/dbVarCNVsBenign.unique.bed
cat files/dbVarCNGsPathogenic.txt | sort-bed - | uniq > files/dbVarCNGsPathogenic.unique.bed
cat files/dbVarCNLsPathogenic.txt | sort-bed - | uniq > files/dbVarCNLsPathogenic.unique.bed
cat files/dbVarCNVsPathogenic.txt | sort-bed - | uniq > files/dbVarCNVsPathogenic.unique.bed

wc -l files/dbVarCNGsBenign.unique.bed
wc -l files/dbVarCNLsBenign.unique.bed
wc -l files/dbVarCNVsBenign.unique.bed
wc -l files/dbVarCNGsPathogenic.unique.bed
wc -l files/dbVarCNLsPathogenic.unique.bed
wc -l files/dbVarCNVsPathogenic.unique.bed

echo "Creating CNV regions (CNVRs)..."
bedtools merge -i files/dbVarCNGsBenign.unique.bed > files/dbVarCNGsBenign.regions.bed
bedtools merge -i files/dbVarCNLsBenign.unique.bed > files/dbVarCNLsBenign.regions.bed
bedtools merge -i files/dbVarCNVsBenign.unique.bed > files/dbVarCNVsBenign.regions.bed
bedtools merge -i files/dbVarCNGsPathogenic.unique.bed > files/dbVarCNGsPathogenic.regions.bed
bedtools merge -i files/dbVarCNLsPathogenic.unique.bed > files/dbVarCNLsPathogenic.regions.bed
bedtools merge -i files/dbVarCNVsPathogenic.unique.bed > files/dbVarCNVsPathogenic.regions.bed

wc -l files/dbVarCNGsBenign.regions.bed
wc -l files/dbVarCNLsBenign.regions.bed
wc -l files/dbVarCNVsBenign.regions.bed
wc -l files/dbVarCNGsPathogenic.regions.bed
wc -l files/dbVarCNLsPathogenic.regions.bed
wc -l files/dbVarCNVsPathogenic.regions.bed

echo "Calculating CNV coverage..."
bedops --partition files/dbVarCNGsBenign.bed > files/dbVarCNGsBenign.partitions.bed
bedops --partition files/dbVarCNLsBenign.bed > files/dbVarCNLsBenign.partitions.bed
bedops --partition files/dbVarCNVsBenign.bed > files/dbVarCNVsBenign.partitions.bed
bedops --partition files/dbVarCNGsPathogenic.bed > files/dbVarCNGsPathogenic.partitions.bed
bedops --partition files/dbVarCNLsPathogenic.bed > files/dbVarCNLsPathogenic.partitions.bed
bedops --partition files/dbVarCNVsPathogenic.bed > files/dbVarCNVsPathogenic.partitions.bed

bedtools coverage -a files/dbVarCNGsBenign.bed -b files/dbVarCNGsBenign.partitions.bed > files/dbVarCNGsBenign.coverage.bed
bedtools coverage -a files/dbVarCNLsBenign.bed -b files/dbVarCNLsBenign.partitions.bed > files/dbVarCNLsBenign.coverage.bed
bedtools coverage -a files/dbVarCNVsBenign.bed -b files/dbVarCNVsBenign.partitions.bed > files/dbVarCNVsBenign.coverage.bed
bedtools coverage -a files/dbVarCNGsPathogenic.bed -b files/dbVarCNGsPathogenic.partitions.bed > files/dbVarCNGsPathogenic.coverage.bed
bedtools coverage -a files/dbVarCNLsPathogenic.bed -b files/dbVarCNLsPathogenic.partitions.bed > files/dbVarCNLsPathogenic.coverage.bed
bedtools coverage -a files/dbVarCNVsPathogenic.bed -b files/dbVarCNVsPathogenic.partitions.bed > files/dbVarCNVsPathogenic.coverage.bed

echo "Sort CNV coverage regions..."
cat files/dbVarCNGsBenign.coverage.bed | sort-bed - > files/dbVarCNGsBenign.coverageSorted.bed
cat files/dbVarCNLsBenign.coverage.bed | sort-bed - > files/dbVarCNLsBenign.coverageSorted.bed
cat files/dbVarCNVsBenign.coverage.bed | sort-bed - > files/dbVarCNVsBenign.coverageSorted.bed
cat files/dbVarCNGsPathogenic.coverage.bed | sort-bed - > files/dbVarCNGsPathogenic.coverageSorted.bed
cat files/dbVarCNLsPathogenic.coverage.bed | sort-bed - > files/dbVarCNLsPathogenic.coverageSorted.bed
cat files/dbVarCNVsPathogenic.coverage.bed | sort-bed - > files/dbVarCNVsPathogenic.coverageSorted.bed

echo "Finding peak CNV coverage regions..."

cat files/dbVarCNGsBenign.coverageSorted.bed | python findPeaksInCNVCoverage.py > files/dbVarCNGsBenign.peakRegions.bed
cat files/dbVarCNLsBenign.coverageSorted.bed | python findPeaksInCNVCoverage.py > files/dbVarCNLsBenign.peakRegions.bed
cat files/dbVarCNVsBenign.coverageSorted.bed | python findPeaksInCNVCoverage.py > files/dbVarCNVsBenign.peakRegions.bed
cat files/dbVarCNGsPathogenic.coverageSorted.bed | python findPeaksInCNVCoverage.py > files/dbVarCNGsPathogenic.peakRegions.bed
cat files/dbVarCNLsPathogenic.coverageSorted.bed | python findPeaksInCNVCoverage.py > files/dbVarCNLsPathogenic.peakRegions.bed
cat files/dbVarCNVsPathogenic.coverageSorted.bed | python findPeaksInCNVCoverage.py > files/dbVarCNVsPathogenic.peakRegions.bed

wc -l files/dbVarCNGsBenign.peakRegions.bed
wc -l files/dbVarCNLsBenign.peakRegions.bed
wc -l files/dbVarCNVsBenign.peakRegions.bed
wc -l files/dbVarCNGsPathogenic.peakRegions.bed
wc -l files/dbVarCNLsPathogenic.peakRegions.bed
wc -l files/dbVarCNVsPathogenic.peakRegions.bed

echo "Printing length of CNV regions..."
awk '{ print $4 = $3 - $2 "\tBenign CNG"}' files/dbVarCNGsBenign.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tBenign CNL"}' files/dbVarCNLsBenign.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tBenign CNV"}' files/dbVarCNVsBenign.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tPathogenic CNG"}' files/dbVarCNGsPathogenic.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tPathogenic CNL"}' files/dbVarCNLsPathogenic.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tPathogenic CNV"}' files/dbVarCNVsPathogenic.bed | awk '{ sum += $1 } END { print sum }'

awk '{ print $4 = $3 - $2 "\tBenign CNG"}' files/dbVarCNGsBenign.regions.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tBenign CNL"}' files/dbVarCNLsBenign.regions.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tBenign CNV"}' files/dbVarCNVsBenign.regions.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tPathogenic CNG"}' files/dbVarCNGsPathogenic.regions.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tPathogenic CNL"}' files/dbVarCNLsPathogenic.regions.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tPathogenic CNV"}' files/dbVarCNVsPathogenic.regions.bed | awk '{ sum += $1 } END { print sum }'

awk '{ print $4 = $3 - $2 "\tBenign CNG"}' files/dbVarCNGsBenign.peakRegions.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tBenign CNL"}' files/dbVarCNLsBenign.peakRegions.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tBenign CNV"}' files/dbVarCNVsBenign.peakRegions.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tPathogenic CNG"}' files/dbVarCNGsPathogenic.peakRegions.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tPathogenic CNL"}' files/dbVarCNLsPathogenic.peakRegions.bed | awk '{ sum += $1 } END { print sum }'
awk '{ print $4 = $3 - $2 "\tPathogenic CNV"}' files/dbVarCNVsPathogenic.peakRegions.bed | awk '{ sum += $1 } END { print sum }'


echo "Creating PC and dev gene list files..."
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNGsBenign.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNGsBenignPCGenes.txt
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNLsBenign.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNLsBenignPCGenes.txt
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNVsBenign.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNVsBenignPCGenes.txt
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNGsPathogenic.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNGsPathogenicPCGenes.txt
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNLsPathogenic.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNLsPathogenicPCGenes.txt
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNVsPathogenic.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNVsPathogenicPCGenes.txt

bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNGsBenign.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNGsBenignPCDevGenes.txt
bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNLsBenign.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNLsBenignPCDevGenes.txt
bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNVsBenign.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNVsBenignPCDevGenes.txt
bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNGsPathogenic.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNGsPathogenicPCDevGenes.txt
bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNLsPathogenic.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNLsPathogenicPCDevGenes.txt
bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNVsPathogenic.regions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNVsPathogenicPCDevGenes.txt

bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNGsBenign.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNGPeaksBenignPCGenes.txt
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNLsBenign.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNLPeaksBenignPCGenes.txt
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNVsBenign.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNVPeaksBenignPCGenes.txt
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNGsPathogenic.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNGPeaksPathogenicPCGenes.txt
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNLsPathogenic.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNLPeaksPathogenicPCGenes.txt
bedtools intersect -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNVsPathogenic.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNVPeaksPathogenicPCGenes.txt

bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNGsBenign.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNGPeaksBenignPCDevGenes.txt
bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNLsBenign.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNLPeaksBenignPCDevGenes.txt
bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNVsBenign.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNVPeaksBenignPCDevGenes.txt
bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNGsPathogenic.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNGPeaksPathogenicPCDevGenes.txt
bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNLsPathogenic.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNLPeaksPathogenicPCDevGenes.txt
bedtools intersect -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed -b files/dbVarCNVsPathogenic.peakRegions.bed -u | awk '{print $4}' > files/gene\ lists/dbVarCNVPeaksPathogenicPCDevGenes.txt

echo "Number of PC and dev genes in CNV..."

wc -l files/gene\ lists/dbVarCNGsBenignPCGenes.txt
wc -l files/gene\ lists/dbVarCNLsBenignPCGenes.txt
wc -l files/gene\ lists/dbVarCNVsBenignPCGenes.txt
wc -l files/gene\ lists/dbVarCNGsPathogenicPCGenes.txt
wc -l files/gene\ lists/dbVarCNLsPathogenicPCGenes.txt
wc -l files/gene\ lists/dbVarCNVsPathogenicPCGenes.txt

wc -l files/gene\ lists/dbVarCNGsBenignPCDevGenes.txt
wc -l files/gene\ lists/dbVarCNLsBenignPCDevGenes.txt
wc -l files/gene\ lists/dbVarCNVsBenignPCDevGenes.txt
wc -l files/gene\ lists/dbVarCNGsPathogenicPCDevGenes.txt
wc -l files/gene\ lists/dbVarCNLsPathogenicPCDevGenes.txt
wc -l files/gene\ lists/dbVarCNVsPathogenicPCDevGenes.txt

wc -l files/gene\ lists/dbVarCNGPeaksBenignPCGenes.txt
wc -l files/gene\ lists/dbVarCNLPeaksBenignPCGenes.txt
wc -l files/gene\ lists/dbVarCNVPeaksBenignPCGenes.txt
wc -l files/gene\ lists/dbVarCNGPeaksPathogenicPCGenes.txt
wc -l files/gene\ lists/dbVarCNLPeaksPathogenicPCGenes.txt
wc -l files/gene\ lists/dbVarCNVPeaksPathogenicPCGenes.txt

wc -l files/gene\ lists/dbVarCNGPeaksBenignPCDevGenes.txt
wc -l files/gene\ lists/dbVarCNLPeaksBenignPCDevGenes.txt
wc -l files/gene\ lists/dbVarCNVPeaksBenignPCDevGenes.txt
wc -l files/gene\ lists/dbVarCNGPeaksPathogenicPCDevGenes.txt
wc -l files/gene\ lists/dbVarCNLPeaksPathogenicPCDevGenes.txt
wc -l files/gene\ lists/dbVarCNVPeaksPathogenicPCDevGenes.txt

echo "Containing 1+ dev genes..."

bedtools intersect -a files/dbVarCNGsBenign.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNLsBenign.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNVsBenign.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNGsPathogenic.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNLsPathogenic.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNVsPathogenic.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l

bedtools intersect -a files/dbVarCNGsBenign.regions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNLsBenign.regions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNVsBenign.regions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNGsPathogenic.regions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNLsPathogenic.regions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNVsPathogenic.regions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l

bedtools intersect -a files/dbVarCNGsBenign.peakRegions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNLsBenign.peakRegions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNVsBenign.peakRegions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNGsPathogenic.peakRegions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNLsPathogenic.peakRegions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l
bedtools intersect -a files/dbVarCNVsPathogenic.peakRegions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -u | wc -l

# Create files to be read by R for Mann-Whitney U test of dev gene fraction overlap (test of benign & pathogenic CNVs independent of CNV length)
bedtools coverage -b files/dbVarCNVsBenign.bed -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed | cut -f 4,6,7 > files/dbVarCNVsBenign.devGeneOverlap.txt
bedtools subtract -a files/dbVarCNVsPathogenic.bed -b files/dbVarCNVsBenign.regions.bed -A > files/dbVarCNVsPathogenic.withoutBenign.bed
# the number of features in A that overlapped (by at least one base pair) the B interval, the length of the entry in B, the fraction of bases in B that had non-zero coverage from features in A.
bedtools coverage -b files/dbVarCNVsPathogenic.withoutBenign.bed -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed | cut -f 4,6,7 > files/dbVarCNVsPathogenic.withoutBenign.devGeneOverlap.txt

# alternative where normalised by gene count instead of proportion of length
bedtools intersect -a files/dbVarCNVsPathogenic.withoutBenign.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -c > files/dbVarCNVsPathogenic.withoutBenign.devGeneCount.txt
bedtools intersect -a files/dbVarCNVsPathogenic.withoutBenign.bed -b datasets/ensGRCh37PCGenes1-22.sorted.bed -c > files/dbVarCNVsPathogenic.withoutBenign.geneCount.txt
cut -f 4 files/dbVarCNVsPathogenic.withoutBenign.devGeneCount.txt | paste files/dbVarCNVsPathogenic.withoutBenign.geneCount.txt - > files/dbVarCNVsPathogenic.withoutBenign.geneCountWithDevCount.txt

bedtools intersect -a files/dbVarCNVsBenign.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -c > files/dbVarCNVsBenign.devGeneCount.txt
bedtools intersect -a files/dbVarCNVsBenign.bed -b datasets/ensGRCh37PCGenes1-22.sorted.bed -c > files/dbVarCNVsBenign.geneCount.txt
cut -f 4 files/dbVarCNVsBenign.devGeneCount.txt | paste files/dbVarCNVsBenign.geneCount.txt - > files/dbVarCNVsBenign.geneCountWithDevCount.txt

# Create files to be read by R for Mann-Whitney U test of dev gene fraction overlap for peak regions (test of benign & pathogenic CNVs independent of CNV length)
bedtools subtract -a files/dbVarCNVsPathogenic.peakRegions.bed -b files/dbVarCNVsBenign.regions.bed -A > files/dbVarCNVsPathogenic.peakRegions.withoutBenign.bed
# the number of features in A that overlapped (by at least one base pair) the B interval, the length of the entry in B, the fraction of bases in B that had non-zero coverage from features in A.
bedtools coverage -b files/dbVarCNVsPathogenic.peakRegions.withoutBenign.bed -a datasets/ensGRCh37PCDevGenes1-22.sorted.bed | cut -f 4,6,7 > files/dbVarCNVsPathogenic.peakRegions.withoutBenign.devGeneOverlap.txt

# alternative where normalised by gene count instead of proportion of length
bedtools intersect -a files/dbVarCNVsPathogenic.peakRegions.withoutBenign.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -c > files/dbVarCNVsPathogenic.peakRegions.withoutBenign.devGeneCount.txt
bedtools intersect -a files/dbVarCNVsPathogenic.peakRegions.withoutBenign.bed -b datasets/ensGRCh37PCGenes1-22.sorted.bed -c > files/dbVarCNVsPathogenic.peakRegions.withoutBenign.geneCount.txt
cut -f 4 files/dbVarCNVsPathogenic.peakRegions.withoutBenign.devGeneCount.txt | paste files/dbVarCNVsPathogenic.peakRegions.withoutBenign.geneCount.txt - > files/dbVarCNVsPathogenic.peakRegions.withoutBenign.geneCountWithDevCount.txt

bedtools intersect -a files/dbVarCNVsBenign.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -c > files/dbVarCNVsBenign.devGeneCount.txt
bedtools intersect -a files/dbVarCNVsBenign.bed -b datasets/ensGRCh37PCGenes1-22.sorted.bed -c > files/dbVarCNVsBenign.geneCount.txt
cut -f 4 files/dbVarCNVsBenign.devGeneCount.txt | paste files/dbVarCNVsBenign.geneCount.txt - > files/dbVarCNVsBenign.geneCountWithDevCount.txt

bedtools intersect -a files/dbVarCNVsBenign.regions.bed -b datasets/ensGRCh37PCDevGenes1-22.sorted.bed -c > files/dbVarCNVsBenign.regions.devGeneCount.txt
bedtools intersect -a files/dbVarCNVsBenign.regions.bed -b datasets/ensGRCh37PCGenes1-22.sorted.bed -c > files/dbVarCNVsBenign.regions.geneCount.txt
cut -f 4 files/dbVarCNVsBenign.regions.devGeneCount.txt | paste files/dbVarCNVsBenign.regions.geneCount.txt - > files/dbVarCNVsBenign.regions.geneCountWithDevCount.txt


echo "(CNV/Peak) Regions containing just one non-passenger gene"
bedtools subtract -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNGsBenign.regions.bed -f 1 > files/dbVarCNGs.nonPassengerGenes.bed
bedtools subtract -a datasets/ensGRCh37PCGenes1-22.sorted.bed -b files/dbVarCNLsBenign.regions.bed -f 1 > files/dbVarCNLs.nonPassengerGenes.bed

# Region containing just one non-passenger gene
bedtools coverage -counts -a files/dbVarCNGs.nonPassengerGenes.bed -b files/dbVarCNGsPathogenic.regions.bed | grep -E "\t1$" | wc -l
bedtools coverage -counts -a files/dbVarCNGs.nonPassengerGenes.bed -b files/dbVarCNGsPathogenic.peakRegions.bed | grep -E "\t1$" | wc -l
bedtools coverage -counts -a files/dbVarCNLs.nonPassengerGenes.bed -b files/dbVarCNLsPathogenic.regions.bed | grep -E "\t1$" | wc -l
bedtools coverage -counts -a files/dbVarCNLs.nonPassengerGenes.bed -b files/dbVarCNLsPathogenic.peakRegions.bed | grep -E "\t1$" | wc -l

# Gene lists
bedtools coverage -counts -a files/dbVarCNGs.nonPassengerGenes.bed -b files/dbVarCNGsPathogenic.peakRegions.bed | grep -E "\t1$" | bedtools intersect -a files/dbVarCNGs.nonPassengerGenes.bed -b stdin -u >> files/solitaryNonPassengerGenes.bed
bedtools coverage -counts -a files/dbVarCNLs.nonPassengerGenes.bed -b files/dbVarCNLsPathogenic.peakRegions.bed | grep -E "\t1$" | bedtools intersect -a files/dbVarCNLs.nonPassengerGenes.bed -b stdin -u >> files/solitaryNonPassengerGenes.bed
cat files/solitaryNonPassengerGenes.bed | sort-bed - | uniq > files/solitaryNonPassengerGenes.unique.bed

echo "Segmenting genes between types of CNVs..."

python createVennDiagramGeneLists.py
wc -l files/gene\ lists/vennDiagram/BG.txt
wc -l files/gene\ lists/vennDiagram/BL.txt
wc -l files/gene\ lists/vennDiagram/PG.txt
wc -l files/gene\ lists/vennDiagram/PL.txt
wc -l files/gene\ lists/vennDiagram/BGBL.txt
wc -l files/gene\ lists/vennDiagram/BGPG.txt
wc -l files/gene\ lists/vennDiagram/BGPL.txt
wc -l files/gene\ lists/vennDiagram/BLPG.txt
wc -l files/gene\ lists/vennDiagram/BLPL.txt
wc -l files/gene\ lists/vennDiagram/PGPL.txt
wc -l files/gene\ lists/vennDiagram/BGBLPG.txt
wc -l files/gene\ lists/vennDiagram/BGBLPL.txt
wc -l files/gene\ lists/vennDiagram/BGPGPL.txt
wc -l files/gene\ lists/vennDiagram/BLPGPL.txt
wc -l files/gene\ lists/vennDiagram/BGBLPGPL.txt

echo "Create gene groups from merged Venn diagram segments..."
cat files/gene\ lists/vennDiagram/BG.txt files/gene\ lists/vennDiagram/BL.txt files/gene\ lists/vennDiagram/BGBL.txt > files/gene\ lists/vennDiagram/benignOnly.txt
cat files/gene\ lists/vennDiagram/PG.txt files/gene\ lists/vennDiagram/PL.txt files/gene\ lists/vennDiagram/PGPL.txt > files/gene\ lists/vennDiagram/pathogenicOnly.txt
cat files/gene\ lists/vennDiagram/BGPG.txt files/gene\ lists/vennDiagram/BLPL.txt files/gene\ lists/vennDiagram/BGBLPG.txt files/gene\ lists/vennDiagram/BGBLPL.txt files/gene\ lists/vennDiagram/BGPGPL.txt files/gene\ lists/vennDiagram/BLPGPL.txt files/gene\ lists/vennDiagram/BGBLPGPL.txt > files/gene\ lists/vennDiagram/passengers.txt
cat files/gene\ lists/vennDiagram/BLPG.txt > files/gene\ lists/vennDiagram/hypothesisedAggProne.txt
cat files/gene\ lists/vennDiagram/BGPL.txt > files/gene\ lists/vennDiagram/hypothesisedHaploinsufficient.txt

echo "Testing segments for enrichment of gene categories..."

echo " - Dev genes"
sort datasets/ensGRCh37PCDevGenes1-22.txt files/gene\ lists/vennDiagram/benignOnly.txt | uniq -d | wc -l
sort datasets/ensGRCh37PCDevGenes1-22.txt files/gene\ lists/vennDiagram/pathogenicOnly.txt | uniq -d | wc -l
sort datasets/ensGRCh37PCDevGenes1-22.txt files/gene\ lists/vennDiagram/passengers.txt | uniq -d | wc -l
sort datasets/ensGRCh37PCDevGenes1-22.txt files/gene\ lists/vennDiagram/hypothesisedAggProne.txt | uniq -d | wc -l
sort datasets/ensGRCh37PCDevGenes1-22.txt files/gene\ lists/vennDiagram/hypothesisedHaploinsufficient.txt | uniq -d | wc -l

echo " - Ohnologs"
sort datasets/ohnologs.txt files/gene\ lists/vennDiagram/benignOnly.txt | uniq -d | wc -l
sort datasets/ohnologs.txt files/gene\ lists/vennDiagram/pathogenicOnly.txt | uniq -d | wc -l
sort datasets/ohnologs.txt files/gene\ lists/vennDiagram/passengers.txt | uniq -d | wc -l
sort datasets/ohnologs.txt files/gene\ lists/vennDiagram/hypothesisedAggProne.txt | uniq -d | wc -l
sort datasets/ohnologs.txt files/gene\ lists/vennDiagram/hypothesisedHaploinsufficient.txt | uniq -d | wc -l

echo " - Haploinsufficiency scores from ExAC"
perl -e 'open(A,"files/gene\ lists/vennDiagram/benignOnly.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithExACpLI.txt > files/haploinsufficiencyScores/benignOnly.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/pathogenicOnly.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithExACpLI.txt > files/haploinsufficiencyScores/pathogenicOnly.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/passengers.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithExACpLI.txt > files/haploinsufficiencyScores/passengers.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/hypothesisedAggProne.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithExACpLI.txt > files/haploinsufficiencyScores/hypothesisedAggProne.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/hypothesisedHaploinsufficient.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithExACpLI.txt > files/haploinsufficiencyScores/hypothesisedHaploinsufficient.txt

echo " - Gene expression from GTEx"
perl -e 'open(A,"files/gene\ lists/vennDiagram/benignOnly.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithGTExMaxExp.txt > files/maxExpression/benignOnly.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/pathogenicOnly.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithGTExMaxExp.txt > files/maxExpression/pathogenicOnly.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/passengers.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithGTExMaxExp.txt > files/maxExpression/passengers.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/hypothesisedAggProne.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithGTExMaxExp.txt > files/maxExpression/hypothesisedAggProne.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/hypothesisedHaploinsufficient.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithGTExMaxExp.txt > files/maxExpression/hypothesisedHaploinsufficient.txt

perl -e 'open(A,"files/gene\ lists/vennDiagram/benignOnly.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithGTExMedianExp.txt > files/medianExpression/benignOnly.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/pathogenicOnly.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithGTExMedianExp.txt > files/medianExpression/pathogenicOnly.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/passengers.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithGTExMedianExp.txt > files/medianExpression/passengers.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/hypothesisedAggProne.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithGTExMedianExp.txt > files/medianExpression/hypothesisedAggProne.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/hypothesisedHaploinsufficient.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' datasets/ensGRCh37WithGTExMedianExp.txt > files/medianExpression/hypothesisedHaploinsufficient.txt

echo " - Protein complex members"
cut -f 8,9 datasets/uniprot-reviewed%3Ayes+AND+organism%3A-Human+%5B9606%5D-.tab | grep "\S\t\S" | grep -i -E "(?:\S+mer(?:s|ic|iz\S+|\W)|complex|self-associates)" | cut -f 2 | tr \; '\n' | sort | uniq > datasets/ensTranscriptsInvolvedInUniprotProteinComplexes.txt
sort datasets/ensGenesInvolvedInUniportProteinComplexes.txt files/gene\ lists/vennDiagram/benignOnly.txt | uniq -d | wc -l
sort datasets/ensGenesInvolvedInUniportProteinComplexes.txt files/gene\ lists/vennDiagram/pathogenicOnly.txt | uniq -d | wc -l
sort datasets/ensGenesInvolvedInUniportProteinComplexes.txt files/gene\ lists/vennDiagram/passengers.txt | uniq -d | wc -l
sort datasets/ensGenesInvolvedInUniportProteinComplexes.txt files/gene\ lists/vennDiagram/hypothesisedAggProne.txt | uniq -d | wc -l
sort datasets/ensGenesInvolvedInUniportProteinComplexes.txt files/gene\ lists/vennDiagram/hypothesisedHaploinsufficient.txt | uniq -d | wc -l

echo "Copy number analysis of gene segments..."
python copyNumberAnalysis.py | cut -f 1,2,3,4 > files/speciesCopyNumberAnalysis.txt
perl -e 'open(A,"datasets/ensGRCh37PCGenes1-Y.ids.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.txt > files/speciesCopyNumberAnalysis.PCGenes1-Y.txt

perl -e 'open(A,"files/gene\ lists/vennDiagram/BG.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BG.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BL.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BL.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/PG.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/PG.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/PL.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/PL.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BGBL.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BGBL.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BGPG.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BGPG.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BGPL.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BGPL.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BLPG.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BLPG.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BLPL.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BLPL.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/PGPL.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/PGPL.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BGBLPG.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BGBLPG.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BGBLPL.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BGBLPL.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BGPGPL.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BGPGPL.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BLPGPL.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BLPGPL.txt
perl -e 'open(A,"files/gene\ lists/vennDiagram/BGBLPGPL.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/BGBLPGPL.txt
perl -e 'open(A,"files/gene\ lists/solitaryNonPassengers.txt"); while(<A>){chomp; $k{$_}++} while(<>){@a=split(/\t/); print if defined $k{$a[0]}}' files/speciesCopyNumberAnalysis.PCGenes1-Y.txt > files/speciesCopyNumberAnalysis/solitaryNonPassengers.txt


echo "Species copy number analysis across genome..."
cat files/speciesCopyNumberAnalysis.PCGenes1-Y.txt | sort > files/speciesCopyNumberAnalysis.sorted.txt
join datasets/ensGRCh37PCGenes1-Y.sorted.txt files/speciesCopyNumberAnalysis.sorted.txt > files/speciesCopyNumberAnalysis.withGenomePosition.txt
awk '{ print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $1 }' files/speciesCopyNumberAnalysis.withGenomePosition.txt > files/speciesCopyNumberAnalysis.withGenomePosition.bed
sort-bed files/speciesCopyNumberAnalysis.withGenomePosition.bed > files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed

echo "start\tdup\tnoortho" > files/1q21.1del.txt
echo "start\tdup\tnoortho" > files/3q29del.txt
echo "start\tdup\tnoortho" > files/7q11.23dup.txt
echo "start\tdup\tnoortho" > files/7q36.1del.txt
echo "start\tdup\tnoortho" > files/7q36.3dup.txt
echo "start\tdup\tnoortho" > files/10q11.22-23dup.txt
echo "start\tdup\tnoortho" > files/15q11.2del.txt
echo "start\tdup\tnoortho" > files/15q11-q13dup.txt
echo "start\tdup\tnoortho" > files/15q13.3del.txt
echo "start\tdup\tnoortho" > files/16p11.2dup.txt
echo "start\tdup\tnoortho" > files/16p12.1del.txt
echo "start\tdup\tnoortho" > files/16p13.11dupdel.txt
echo "start\tdup\tnoortho" > files/17q12del.txt
echo "start\tdup\tnoortho" > files/22q11del.txt

sed -n '/ENSG00000131788/,/ENSG00000203852/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/1q21.1del.txt
sed -n '/ENSG00000145014/,/ENSG00000185621/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/3q29del.txt
sed -n '/ENSG00000241258/,/ENSG00000106178/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/7q11.23dup.txt
sed -n '/ENSG00000221938/,/ENSG00000204946/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/7q36.1del.txt
sed -n '/ENSG00000164778/,/ENSG00000106018/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/7q36.3dup.txt
sed -n '/ENSG00000165511/,/ENSG00000122873/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/10q11.22-23dup.txt
sed -n '/ENSG00000215405/,/ENSG00000214265/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/15q11.2del.txt
sed -n '/ENSG00000215405/,/ENSG00000034053/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/15q11-q13dup.txt
sed -n '/ENSG00000114062/,/ENSG00000140199/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/15q13.3del.txt
sed -n '/ENSG00000177548/,/ENSG00000270466/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/16p11.2dup.txt
sed -n '/ENSG00000188215/,/ENSG00000083093/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/16p12.1del.txt
sed -n '/ENSG00000103342/,/ENSG00000103534/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/16p13.11dupdel.txt
sed -n '/ENSG00000172660/,/ENSG00000108292/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/17q12del.txt
sed -n '/ENSG00000093072/,/ENSG00000169575/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $2 "\t" $5 "\t" $6 }' >> files/22q11del.txt

sed -n '/ENSG00000131788/,/ENSG00000203852/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000145014/,/ENSG00000185621/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000241258/,/ENSG00000106178/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000221938/,/ENSG00000204946/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000164778/,/ENSG00000106018/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000165511/,/ENSG00000122873/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000215405/,/ENSG00000214265/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000215405/,/ENSG00000034053/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000114062/,/ENSG00000140199/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000177548/,/ENSG00000270466/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000188215/,/ENSG00000083093/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000103342/,/ENSG00000103534/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000172660/,/ENSG00000108292/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
sed -n '/ENSG00000093072/,/ENSG00000169575/p' files/speciesCopyNumberAnalysis.withGenomePosition.sorted.bed | awk '{ print $7 "\t" $4 }' >> files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt
cat files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.txt | sort | uniq > files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.unique.txt

grep -Fwf files/genesFromCritRegionInExampleRegions.txt files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.unique.txt > files/genesFromCritRegionInExampleRegionsWithSpeciesUnchanged.txt
grep -vFwf files/genesFromCritRegionInExampleRegions.txt files/genesFromExampleRegionsWithFlankingGenesAndSpeciesUnchanged.unique.txt > files/genesFlankingCritExampleRegionsWithSpeciesUnchanged.txt


# Clustering of dev genes and effect on number of species unchanged (Supp Figure)
split -l 1 -a 3 datasets/ensGRCh37PCDevGenes1-22.bed files/devGenes/F
for file in files/devGenes/F*; do cat "$file" | grep -v -f /dev/stdin datasets/ensGRCh37PCDevGenes1-22.sorted.bed | closest-features --dist --closest --delim "\t" "$file" - >> files/closestDevGene.txt; done;
cat files/closestDevGene.txt | cut -f 4,9 > files/closestDevGeneSubset.txt
sort files/closestDevGeneSubset.txt > files/closestDevGeneSubset.sorted.txt
join files/closestDevGeneSubset.sorted.txt files/speciesCopyNumberAnalysis.sorted.txt > files/closestDevGeneWithUnchanged.txt

# Randomised location of pathogenic CNVs
for (( i=1; i <= 1000; i++ )); do bedtools shuffle -i files/dbVarCNVsPathogenic.bed -g dataset/bedtoolsChromosomes.hg19.bed > files/randomPathoCNVPositions/dbVarCNVsPathogenic.$1.bed; done;

