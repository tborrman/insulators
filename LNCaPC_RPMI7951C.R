library(ggplot2)
suffix <- "__genome__C-40000-iced__chr_allchrom.scaled-1000000_dixonDI_w-2000000"

LNCaPC_R1 <-read.table(paste('ENCODE3-LNCaPC-HindIII-R1__hg19', suffix, sep=""), header=TRUE, sep="\t")
RPMI7951C_R1 <- read.table(paste('ENCODE3-RPMI7951C-HindIII-R1__hg19', suffix, sep=""), header=TRUE, sep="\t")

diff.df <- data.frame(LNCaPC_R1$chrom, LNCaPC_R1$start, LNCaPC_R1$end, abs(LNCaPC_R1$DI - RPMI7951C_R1$DI))
colnames(diff.df) = c("chrom","start","end", "diff")


sorted.diff.df <- diff.df[order(diff.df$diff, decreasing=TRUE),]

LNCaPC_R1_ex <- LNCaPC_R1[LNCaPC_R1$chrom == "chr21",]

LNCaPC_R1_ex <- LNCaPC_R1_ex[LNCaPC_R1_ex$start > 40000000 & LNCaPC_R1_ex$start < 46000000,]

# options(scipen=100)
# png("LNCaPC_chr21_40Mb-50Mb.png", width=5000, height=1000, res=300)
# plot(LNCaPC_R1_ex$start, LNCaPC_R1_ex$DI, col="blue", pch=20, type="o", xlab="chr21",
#      ylab="DI",ylim=c(-650, 600))
# abline(v= 42880001, col="red", lty=2)
# dev.off()


LNCaPC_R1_ex[["sign"]] = ifelse(LNCaPC_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("LNCaPC_chr21_40Mb-50Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(LNCaPC_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 21", y = "Directionality Index", size=12 )
       + theme_bw()
       + ylim(-750, 625)
       + xlim(40000000, 46000000)
       + geom_vline(xintercept = 42880001, linetype=2, colour="blue", size = 1)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()

RPMI7951C_R1_ex <- RPMI7951C_R1[RPMI7951C_R1$chrom == "chr21",]

RPMI7951C_R1_ex <- RPMI7951C_R1_ex[RPMI7951C_R1_ex$start > 40000000 & RPMI7951C_R1_ex$start < 46000000,]

# options(scipen=100)
# png("RPMI7951C_chr21_40Mb-50Mb.png", width=5000, height=1500, res=300)
# plot(RPMI7951C_R1_ex$start, RPMI7951C_R1_ex$DI, col="blue", pch=20, type="o", xlab="chr21",
#      ylab="DI", ylim=c(-650, 600))
# abline(v= 42880001, col="red", lty=2)
# dev.off()


RPMI7951C_R1_ex[["sign"]] = ifelse(RPMI7951C_R1_ex[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("RPMI7951C_chr21_40Mb-50Mb.png", width=5000, height=1000, res=300)
p1 <- (ggplot(RPMI7951C_R1_ex, aes(x= start, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 21", y = "Directionality Index", size=12 )
       + theme_bw()
       + ylim(-750, 625)
       + xlim(40000000, 46000000)
       + geom_vline(xintercept = 42880001, linetype=2, colour="blue", size = 1)
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()
