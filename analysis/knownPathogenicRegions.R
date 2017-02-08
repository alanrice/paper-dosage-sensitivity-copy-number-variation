library(ggplot2)
library(reshape2)
library(scales)
library(grid)

# colours from B. Wong, Nature Methods 8, 441 (2011)
cbPalette <- c("#009E73","#e69f00")

region1 <- read.csv("files/1q21.1del.txt", sep="\t",header=TRUE, row.names=NULL)
region1.m <- melt(region1, id.vars = "start")
svg(filename="figures/1q21.1del.svg", width=6, height=4, pointsize=12)
ggplot(region1.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 146500000, xmax = 148700000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 146620000, xmax = 147580000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(146000000, 149000000, by = 3000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region2 <- read.csv("files/3q29del.txt", sep="\t",header=TRUE, row.names=NULL)
region2.m <- melt(region2, id.vars = "start")
svg(filename="figures/3q29del.svg", width=6, height=4, pointsize=12)
ggplot(region2.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 195800000, xmax = 197300200, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(194500000, 197500000, by = 3000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region3 <- read.csv("files/7q11.23dup.txt", sep="\t",header=TRUE, row.names=NULL)
region3.m <- melt(region3, id.vars = "start")
svg(filename="figures/7q11.23dup.svg", width=6, height=4, pointsize=12)
ggplot(region3.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 72400000, xmax = 74200000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(66000000, 75000000, by = 9000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region4 <- read.csv("files/7q36.3dup.txt", sep="\t",header=TRUE, row.names=NULL)
region4.m <- melt(region4, id.vars = "start")
svg(filename="figures/7q36.3dup.svg", width=6, height=4, pointsize=12)
ggplot(region4.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 157300000, xmax = 158875000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 158765000, xmax = 158875000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(155500000, 158500000, by = 3000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region5 <- read.csv("files/10q11.22-23dup.txt", sep="\t",header=TRUE, row.names=NULL)
region5.m <- melt(region5, id.vars = "start")
svg(filename="figures/10q11.22-23dup.svg", width=6, height=4, pointsize=12)
ggplot(region5.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 46800000, xmax = 52400000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 49360000, xmax = 51600000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(46000000, 59000000, by = 13000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region6 <- read.csv("files/15q11.2del.txt", sep="\t",header=TRUE, row.names=NULL)
region6.m <- melt(region6, id.vars = "start")
svg(filename="figures/15q11.2del.svg", width=6, height=4, pointsize=12)
ggplot(region6.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 22830000, xmax = 23250000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(21000000, 25000000, by = 4000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region7 <- read.csv("files/15q11-q13dup.txt", sep="\t",header=TRUE, row.names=NULL)
region7.m <- melt(region7, id.vars = "start")
svg(filename="figures/15q11-q13dup.svg", width=6, height=4, pointsize=12)
ggplot(region7.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 22400000, xmax = 29000000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 23800000, xmax = 28400000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(21000000, 29000000, by = 8000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region8 <- read.csv("files/15q13.3del.txt", sep="\t",header=TRUE, row.names=NULL)
region8.m <- melt(region8, id.vars = "start")
svg(filename="figures/15q13.3del.svg", width=6, height=4, pointsize=12)
ggplot(region8.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 29200000, xmax = 33000000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 32250000, xmax = 32390000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(26000000, 34000000, by = 8000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region9 <- read.csv("files/16p13.11dupdel.txt", sep="\t",header=TRUE, row.names=NULL)
region9.m <- melt(region9, id.vars = "start")
svg(filename="figures/16p13.11dupdel.svg", width=6, height=4, pointsize=12)
ggplot(region9.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 14900000, xmax = 18550000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 15480000, xmax = 16250000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(13000000, 19000000, by = 6000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region10 <- read.csv("files/16p11.2dup.txt", sep="\t",header=TRUE, row.names=NULL)
region10.m <- melt(region10, id.vars = "start")
svg(filename="figures/16p11.2dup.svg", width=6, height=4, pointsize=12)
ggplot(region10.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 29620000, xmax = 30200000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 29670000, xmax = 30117000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(29000000, 30300000, by = 1300000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region11 <- read.csv("files/16p12.1del.txt", sep="\t",header=TRUE, row.names=NULL)
region11.m <- melt(region11, id.vars = "start")
svg(filename="figures/16p12.1del.svg", width=6, height=4, pointsize=12)
ggplot(region11.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 21800000, xmax = 22400000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 21960000, xmax = 22360000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(21000000, 23500000, by = 2500000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region12 <- read.csv("files/17q12del.txt", sep="\t",header=TRUE, row.names=NULL)
region12.m <- melt(region12, id.vars = "start")
svg(filename="figures/17q12del.svg", width=6, height=4, pointsize=12)
ggplot(region12.m, aes(x=start, y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 34500000, xmax = 36050000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 34840000, xmax = 36050000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(34500000, 36500000, by = 2000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region13 <- read.csv("files/22q11del.txt", sep="\t",header=TRUE, row.names=NULL)
region13.m <- melt(region13, id.vars = "start")
svg(filename="figures/22q11del.svg", width=6, height=4, pointsize=12)
ggplot(region13.m,
  aes(x=start,y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 18890000, xmax = 21800000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 18890000, xmax = 20230000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(18000000, 22000000, by = 4000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

region14 <- read.csv("files/7q36.1del.txt", sep="\t",header=TRUE, row.names=NULL)
region14.m <- melt(region14, id.vars = "start")
svg(filename="figures/7q36.1del.svg", width=6, height=4, pointsize=12)
ggplot(region14.m,
  aes(x=start,y=value,colour=variable,shape=variable)) +
  annotate("rect", xmin = 145750000, xmax = 145880000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  annotate("rect", xmin = 145750000, xmax = 145880000, ymin=-0.5, ymax= 13.5,fill=alpha(c("black"), .1)) +
  geom_point(size = 6, alpha=0.85) +
  geom_line(size=2.5, alpha=0.7) + 
  coord_cartesian(ylim=c(0,13)) +
  scale_x_continuous(breaks = seq(144000000, 149000000, by = 5000000)) +
  scale_y_continuous(breaks=seq(0, 12, 6)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black",size = 3), axis.line.y = element_line(colour = "black",size = 3), legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.length=unit(0.5,"cm"), axis.ticks.x=element_line(colour = "black",size=3), axis.ticks.y=element_line(colour = "black",size=3), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank()) +
  scale_colour_manual(values=cbPalette)
dev.off()

