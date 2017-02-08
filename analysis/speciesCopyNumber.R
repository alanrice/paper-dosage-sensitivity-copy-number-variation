library(ggplot2)
library(reshape2)

genomeAvg <- read.csv("files/speciesCopyNumberAnalysis.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(genomeAvg) <- c("id", "onetoone", "dup", "noortho")
genomeAvgData <- melt(genomeAvg, id=c("id"))
genomeAvgData$segment <- "genomeAvg"

bg <- read.csv("files/speciesCopyNumberAnalysis/BG.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(bg) <- c("id", "onetoone", "dup", "noortho")
bgData <- melt(bg, id=c("id"))
bgData$segment <- "bg"

bl <- read.csv("files/speciesCopyNumberAnalysis/BL.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(bl) <- c("id", "onetoone", "dup", "noortho")
blData <- melt(bl, id=c("id"))
blData$segment <- "bl"

pg <- read.csv("files/speciesCopyNumberAnalysis/PG.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(pg) <- c("id", "onetoone", "dup", "noortho")
pgData <- melt(pg, id=c("id"))
pgData$segment <- "pg"

pl <- read.csv("files/speciesCopyNumberAnalysis/PL.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(pl) <- c("id", "onetoone", "dup", "noortho")
plData <- melt(pl, id=c("id"))
plData$segment <- "pl"

bgbl <- read.csv("files/speciesCopyNumberAnalysis/BGBL.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(bgbl) <- c("id", "onetoone", "dup", "noortho")
bgblData <- melt(bgbl, id=c("id"))
bgblData$segment <- "bgbl"

bgpg <- read.csv("files/speciesCopyNumberAnalysis/BGPG.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(bgpg) <- c("id", "onetoone", "dup", "noortho")
bgpgData <- melt(bgpg, id=c("id"))
bgpgData$segment <- "bgpg"

bgpl <- read.csv("files/speciesCopyNumberAnalysis/BGPL.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(bgpl) <- c("id", "onetoone", "dup", "noortho")
bgplData <- melt(bgpl, id=c("id"))
bgplData$segment <- "bgpl"

blpg <- read.csv("files/speciesCopyNumberAnalysis/BLPG.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(blpg) <- c("id", "onetoone", "dup", "noortho")
blpgData <- melt(blpg, id=c("id"))
blpgData$segment <- "blpg"

blpl <- read.csv("files/speciesCopyNumberAnalysis/BLPL.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(blpl) <- c("id", "onetoone", "dup", "noortho")
blplData <- melt(blpl, id=c("id"))
blplData$segment <- "blpl"

pgpl <- read.csv("files/speciesCopyNumberAnalysis/PGPL.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(pgpl) <- c("id", "onetoone", "dup", "noortho")
pgplData <- melt(pgpl, id=c("id"))
pgplData$segment <- "pgpl"

bgblpg <- read.csv("files/speciesCopyNumberAnalysis/BGBLPG.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(bgblpg) <- c("id", "onetoone", "dup", "noortho")
bgblpgData <- melt(bgblpg, id=c("id"))
bgblpgData$segment <- "bgblpg"

bgblpl <- read.csv("files/speciesCopyNumberAnalysis/BGBLPL.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(bgblpl) <- c("id", "onetoone", "dup", "noortho")
bgblplData <- melt(bgblpl, id=c("id"))
bgblplData$segment <- "bgblpl"

bgpgpl <- read.csv("files/speciesCopyNumberAnalysis/BGPGPL.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(bgpgpl) <- c("id", "onetoone", "dup", "noortho")
bgpgplData <- melt(bgpgpl, id=c("id"))
bgpgplData$segment <- "bgpgpl"

blpgpl <- read.csv("files/speciesCopyNumberAnalysis/BLPGPL.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(blpgpl) <- c("id", "onetoone", "dup", "noortho")
blpgplData <- melt(blpgpl, id=c("id"))
blpgplData$segment <- "blpgpl"

bgblpgpl <- read.csv("files/speciesCopyNumberAnalysis/BGBLPGPL.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(bgblpgpl) <- c("id", "onetoone", "dup", "noortho")
bgblpgplData <- melt(bgblpgpl, id=c("id"))
bgblpgplData$segment <- "bgblpgpl"

passengers <- rbind(bgpgData, blplData, bgblpgData, bgblplData, bgpgplData, blpgplData, bgblpgplData)
passengers$segment <- "passengers"

solitaries <- read.csv("files/speciesCopyNumberAnalysis/solitaryNonPassengers.txt", sep="\t",header=FALSE, row.names=NULL)
colnames(solitaries) <- c("id", "onetoone", "dup", "noortho")
solitariesData <- melt(solitaries, id=c("id"))
solitariesData$segment <- "solitaries"

# colours from B. Wong, Nature Methods 8, 441 (2011)

total <- rbind(genomeAvgData,bgData, blData, pgData, plData, bgplData, blpgData, bgblData, pgplData,passengers, solitariesData)
svg(filename="figures/speciesCopyNumber.svg", width=9.5, height=4, pointsize=18)
ggplot(total, aes(x = factor(segment), y = value,fill=variable,color=variable)) +
  geom_point(position=position_jitterdodge(dodge.width=0.8,jitter.height=0.25),size = 0.5) +
  geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.8), notch = TRUE, notchwidth = 0.5) +
  coord_cartesian(ylim=c(-0.31,13.31)) +
  scale_y_continuous(breaks=seq(0, 13, 2)) +
  scale_color_manual(values=c("#000000","#009c73", "#e69f00"),name="Status",limits = c("onetoone","dup","noortho"),labels=c("Unchanged","Duplicated", "No ortholog")) + scale_fill_manual(values=c("#000000","#009c73", "#e69f00"),name="Status",limits = c("onetoone","dup","noortho"),labels=c("Unchanged","Duplicated", "No ortholog")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=14), axis.title.x = element_blank(), legend.position="bottom") +
  labs(y = "Number of species") +
  scale_x_discrete(limits = c("genomeAvg", "solitaries","pg","pl","pgpl","bg","bl","bgbl","bgpl", "blpg", "passengers"))
dev.off()

