library(ggplot2)
suffix <- "__genome__C-40000-iced__chr_allchrom.scaled-1000000_dixonDI_w-2000000"

PANC1B_R1 <-read.table(paste('ENCODE3-PANC1B-HindIII-R1__hg19', suffix, sep=""), header=TRUE, sep="\t")
G401A_R1 <- read.table(paste('ENCODE3-G401A-HindIII-R1__hg19', suffix, sep=""), header=TRUE, sep="\t")

diff.df <- data.frame(PANC1B_R1$chrom, PANC1B_R1$start, PANC1B_R1$end, abs(PANC1B_R1$DI - G401A_R1$DI))
colnames(diff.df) = c("chrom","start","end", "diff")


sorted.diff.df <- diff.df[order(diff.df$diff, decreasing=TRUE),]

PANC1B_R1_ex <- PANC1B_R1[PANC1B_R1$chrom == "chr21",]

PANC1B_R1_ex <- PANC1B_R1_ex[PANC1B_R1_ex$start > 40000000 & PANC1B_R1_ex$start < 50000000,]

PANC1B_R1_ex[["sign"]] = ifelse(PANC1B_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("PANC1B_chr21_40Mb-50Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(PANC1B_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 21", y = "Directionality Index", size=12 )
       + theme_bw()
       + geom_vline(xintercept = 40480001, linetype=2, colour="blue", size = 1)
       + xlim(40000000, 50000000)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()


G401A_R1_ex <- G401A_R1[G401A_R1$chrom == "chr21",]

G401A_R1_ex <- G401A_R1_ex[G401A_R1_ex$start > 40000000 & G401A_R1_ex$start < 50000000,]

G401A_R1_ex[["sign"]] = ifelse(G401A_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("G401A_chr21_40Mb-50Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(G401A_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 21", y = "Directionality Index", size=12 )
       + theme_bw()
       + geom_vline(xintercept = 40480001, linetype=2, colour="blue", size = 1)
       + xlim(40000000, 50000000)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()


###########################################################################################

PANC1B_R1_ex <- PANC1B_R1[PANC1B_R1$chrom == "chr22",]

PANC1B_R1_ex <- PANC1B_R1_ex[PANC1B_R1_ex$start > 30000000 & PANC1B_R1_ex$start < 40000000,]

PANC1B_R1_ex[["sign"]] = ifelse(PANC1B_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("PANC1B_chr22_30Mb-40Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(PANC1B_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 22", y = "Directionality Index", size=12 )
       + theme_bw()
       + geom_vline(xintercept = 35520001, linetype=2, colour="blue", size = 1)
       + xlim(30000000, 40000000)
       + ylim(-700, 500)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()


G401A_R1_ex <- G401A_R1[G401A_R1$chrom == "chr22",]

G401A_R1_ex <- G401A_R1_ex[G401A_R1_ex$start > 30000000 & G401A_R1_ex$start < 40000000,]

G401A_R1_ex[["sign"]] = ifelse(G401A_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("G401A_chr22_30Mb-40Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(G401A_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 22", y = "Directionality Index", size=12 )
       + theme_bw()
       + geom_vline(xintercept = 35520001, linetype=2, colour="blue", size = 1)
       + xlim(30000000, 40000000)
       + ylim(-700, 500)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()
####################################################################################
PANC1B_G401 <- cbind(PANC1B_R1, G401A_R1$DI)
colnames(PANC1B_G401) <- c("chrom", "start", "end", "PANC1B_DI", "G401A_DI")
opp_signs <- PANC1B_G401[(PANC1B_G401$PANC1B_DI > 0 & PANC1B_G401$G401A_DI < 0) | (PANC1B_G401$PANC1B_DI < 0 & PANC1B_G401$G401A_DI > 0),]
opp_signs <- cbind(opp_signs,abs(opp_signs$PANC1B_DI - opp_signs$G401A_DI))
opp_signs <- opp_signs[abs(opp_signs$PANC1B_DI) > 10 & abs(opp_signs$G401A_DI) > 10,]
colnames(opp_signs) =  c("chrom", "start", "end", "PANC1B_DI", "G401A_DI", "diff")


sorted.opp_signs.df <- opp_signs[order(opp_signs$diff, decreasing=TRUE),]

PANC1B_R1_ex <- PANC1B_R1[PANC1B_R1$chrom == "chr14",]

PANC1B_R1_ex <- PANC1B_R1_ex[PANC1B_R1_ex$start > 70000000 & PANC1B_R1_ex$start < 80000000,]

PANC1B_R1_ex[["sign"]] = ifelse(PANC1B_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("PANC1B_chr14_70Mb-80Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(PANC1B_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 14", y = "Directionality Index", size=12 )
       + theme_bw()
       + geom_vline(xintercept = 72320001, linetype=2, colour="blue", size = 1)
       + xlim(70000000, 80000000)
       + ylim(-150, 175)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()


G401A_R1_ex <- G401A_R1[G401A_R1$chrom == "chr14",]

G401A_R1_ex <- G401A_R1_ex[G401A_R1_ex$start > 70000000 & G401A_R1_ex$start < 80000000,]

G401A_R1_ex[["sign"]] = ifelse(G401A_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("G401A_chr14_70Mb-80Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(G401A_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 14", y = "Directionality Index", size=12 )
       + theme_bw()
       + geom_vline(xintercept = 72320001, linetype=2, colour="blue", size = 1)
       + xlim(70000000, 80000000)
       + ylim(-150, 175)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()

PANC1B_R1_ex <- PANC1B_R1[PANC1B_R1$chrom == "chr19",]

PANC1B_R1_ex <- PANC1B_R1_ex[PANC1B_R1_ex$start > 38000000 & PANC1B_R1_ex$start < 50000000,]

PANC1B_R1_ex[["sign"]] = ifelse(PANC1B_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("PANC1B_chr19_38Mb-50Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(PANC1B_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 19", y = "Directionality Index", size=12 )
       + theme_bw()
       + geom_vline(xintercept = c(49040001, 38920001,40800001), linetype=rep(2,3), colour=rep("blue",3), size = rep(1,3))
       + xlim(38000000, 50000000)
       + ylim(-500, 500)
       #+ ylim(-150, 175)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()

G401A_R1_ex <- G401A_R1[G401A_R1$chrom == "chr19",]

G401A_R1_ex <- G401A_R1_ex[G401A_R1_ex$start > 38000000 & G401A_R1_ex$start < 50000000,]

G401A_R1_ex[["sign"]] = ifelse(G401A_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("G401A_chr19_38Mb-50Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(G401A_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 19", y = "Directionality Index", size=12 )
       + theme_bw()
       + geom_vline(xintercept = c(49040001, 38920001,40800001), linetype=rep(2,3), colour=rep("blue",3), size = rep(1,3))
       + xlim(38000000, 50000000)
       + ylim(-500, 500)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()

PANC1B_R1_ex <- PANC1B_R1[PANC1B_R1$chrom == "chr5",]

PANC1B_R1_ex <- PANC1B_R1_ex[PANC1B_R1_ex$start > 150000000 & PANC1B_R1_ex$start < 160000000,]

PANC1B_R1_ex[["sign"]] = ifelse(PANC1B_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("PANC1B_chr5_150Mb-160Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(PANC1B_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 5", y = "Directionality Index", size=12 )
       + theme_bw()
       + geom_vline(xintercept = 154360001, linetype=2, colour="blue", size = 1)
       + xlim(150000000, 160000000)
       + ylim(-60, 75)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()






G401A_R1_ex <- G401A_R1[G401A_R1$chrom == "chr5",]

G401A_R1_ex <- G401A_R1_ex[G401A_R1_ex$start > 150000000 & G401A_R1_ex$start < 160000000,]

G401A_R1_ex[["sign"]] = ifelse(G401A_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("G401A_chr5_150Mb-160Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(G401A_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 5", y = "Directionality Index", size=12 )
       + geom_vline(xintercept = 154360001, linetype=2, colour="blue", size = 1)
       + theme_bw()
       + xlim(150000000, 160000000)
       + ylim(-60, 75)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()

PANC1B_R1_ex <- PANC1B_R1[PANC1B_R1$chrom == "chr6",]

PANC1B_R1_ex <- PANC1B_R1_ex[PANC1B_R1_ex$start > 100000001 & PANC1B_R1_ex$start < 110000000,]

PANC1B_R1_ex[["sign"]] = ifelse(PANC1B_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("PANC1B_chr6_100Mb-110Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(PANC1B_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 6", y = "Directionality Index", size=12 )
       + theme_bw()
       + geom_vline(xintercept = 105000001, linetype=2, colour="blue", size = 1)
       + xlim(100000000, 110000000)
       + ylim(-35, 50)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()



G401A_R1_ex <- G401A_R1[G401A_R1$chrom == "chr6",]

G401A_R1_ex <- G401A_R1_ex[G401A_R1_ex$start > 100000000 & G401A_R1_ex$start < 110000000,]

G401A_R1_ex[["sign"]] = ifelse(G401A_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("G401A_chr6_100Mb-110Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(G401A_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 6", y = "Directionality Index", size=12 )
       + geom_vline(xintercept = 105000001, linetype=2, colour="blue", size = 1)
       + theme_bw()
       + xlim(100000000, 110000000)
       + ylim(-35, 50)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()

