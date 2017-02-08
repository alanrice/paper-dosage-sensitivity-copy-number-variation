# Overlap between gene categories and Venn diagram groups

# Developmental genes
M <- as.table(rbind(c(154, 1606, 117, 19, 28), c(1075-154, 6367-1606, 524-117, 94-19, 110-28)))
dimnames(M) <- list(category = c("Dev", "NonDev"), group = c("Benign","Pathogenic", "Passenger", "AggPro", "Hap"))
(Xsq <- chisq.test(M))
Xsq$observed
Xsq$expected
Xsq$stdres

# Ohnologs
M <- as.table(rbind(c(289, 2395, 159, 28, 46), c(1075-289, 6367-2395, 524-159, 94-28, 110-46)))
dimnames(M) <- list(category = c("Ohno", "NonOhno"), group = c("Benign","Pathogenic", "Passenger", "AggPro", "Hap"))
(Xsq <- chisq.test(M))
Xsq$observed
Xsq$expected
Xsq$stdres

# Uniprot Protein complex membership
M <- as.table(rbind(c(251, 2156, 150, 31, 39), c(1075-251, 6367-2156, 524-150, 94-31, 110-39)))
dimnames(M) <- list(category = c("Complex Members", "Non Complex Members"), group = c("Benign","Pathogenic", "Passenger", "AggPro", "Hap"))
(Xsq <- chisq.test(M))
Xsq$observed
Xsq$expected
Xsq$stdres
Xsq$p.value

# Species copy number - proportion of unchanged across mammalian species
M <- as.table(rbind(c(1196, 994, 295, 154, 43, 14, 40, 26, 116), c(2960-1196, 2383-994, 765-295, 617-154, 186-43, 183-14, 104-40, 90-26, 485-116)))
dimnames(M) <- list(category = c("Totally conserved", "Not tot conserved"), group = c("PG","PL", "PGPL", "BG", "BL", "BGBL", "BGPL", "BLPG", "passengers"))
(Xsq <- chisq.test(M))

# Difference in variance between groups
total <-rbind(bgData, blData, pgData, plData, bgplData, blpgData, bgblData, pgplData,passengers)
totalSubset <- subset(total, variable == "onetoone")
totalSubset$segment <- as.factor(totalSubset$segment)
fl <- fligner.test(value ~ segment, data = totalSubset)
fl$p.value

# Example known pathogenic regions from Figure 3
# Comparing median and variance of species with unchanged copy number between critical regions (where applicable) and flanking regions
genesFlankingCritExampleRegions <- read.csv("files/genesFlankingCritExampleRegionsWithSpeciesUnchanged.txt", sep="\t", header=FALSE)
genesFlankingCritExampleRegions$status <- "flanking"
summary(genesFlankingCritExampleRegions)
genesFromCritRegionInExampleRegions <- read.csv("files/genesFromCritRegionInExampleRegionsWithSpeciesUnchanged.txt", sep="\t", header=FALSE)
genesFromCritRegionInExampleRegions$status <- "within"
summary(genesFromCritRegionInExampleRegions)
genesFromCritExampleRegionsCombined <- rbind(genesFlankingCritExampleRegions, genesFromCritRegionInExampleRegions)
summary(genesFromCritExampleRegionsCombined)                                
ggplot(genesFromCritExampleRegionsCombined, aes(factor(status), V2)) + geom_boxplot()
wilcox.test(V2 ~ factor(status), data=genesFromCritExampleRegionsCombined)
fligner.test(V2 ~ factor(status), data = genesFromCritExampleRegionsCombined)
var(genesFlankingCritExampleRegions$V2)
var(genesFromCritRegionInExampleRegions$V2)


# 7014 genes with conserved copy number across mammals

# Compare enrichment of OMIM among totally conserved genes versus rest of genome with copy number changes in other species
M <- as.table(rbind(c(1273, 1637), c(7014-1273, 12246-1637)))
dimnames(M) <- list(category = c("OMIM Disease", "Not OMIM"), group = c("Totally conserved","Rest of genome"))
(Xsq <- chisq.test(M))


# Overlap of totally conserved genes and benign CNVs (only considering protein coding genes on Chr1-22 [as that CNVs not considered on sex chromosomes] and genes included in species copy number analysis)
M <- as.table(rbind(c(393, 6416), c(1272, 10360)))
> dimnames(M) <- list(category = c("Totally conserved","Not totally conserved"), group = c("Overlapped by Benign CNV", "Not overlapped by Benign CNV"))
> (Xsq <- chisq.test(M))

#Overlap of totally conserved genes and inclusive CNV map from Zarrei et al. 2015
M <- as.table(rbind(c(2500, 4514), c(5104, 8196)))
dimnames(M) <- list(category = c("Totally conserved","Not totally conserved"), group = c("Overlapped by IncMap", "Not overlapped by IncMap"))
(Xsq <- chisq.test(M))

#Overlap of mouse orthologs of totally conserved genes and mouse CNVs
M <- as.table(rbind(c(1632, 6928-1632), c(5904-1632, 22559-6928-5904+1632)))
dimnames(M) <- list(category = c("Conserved ortho","Not conserved ortho"), group = c("Overlapped by CNV", "Not overlapped by CNV"))
(Xsq <- chisq.test(M))

#Genes excluded from species copy number analysis
M <- as.table(rbind(c(137, 852), c(1665, 16776)))
dimnames(M) <- list(category = c("Excluded","Not excluded"), group = c("Benign", "Not benign"))
(Xsq <- chisq.test(M))


# p-value for 1000 random simulations (Supplementary Fig. 3)
pnorm(0.954, mean=0.7487162, sd=0.005522641, lower.tail=FALSE) * 2
#[1] 1.984779e-302
pnorm(0.228967, mean=0.2264564, sd=0.0006062768, lower.tail=FALSE) * 2
#[1] 3.457756e-05
