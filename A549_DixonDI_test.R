library(ggplot2)
A549_R1 <-read.table("ENCODE3-A549C-HindIII-R1__hg19__genome__C-40000-iced__chr2_dixonDI_w-2000000", comment.char="#", header=TRUE, sep="\t")
A549_R1_win <- A549_R1[A549_R1$start > 137000000 & A549_R1$start < 141000000,]
midpoint <-  (A549_R1_win$start + A549_R1_win$end)/2
A549_R1_win <- cbind(A549_R1_win, midpoint)

A549_R1_win[["sign"]] = ifelse(A549_R1_win[["DI"]] >= 0, "positive", "negative")

options(scipen=100)
png("A549_chr2_137Mb-141Mb_DixonDI.png", width=5000, height=1000, res=300)
p1 <- (ggplot(A549_R1_win, aes(x= midpoint, y=DI, fill=sign)) 
       + geom_bar(stat='identity')
       + scale_fill_manual(values = c("positive" = "red", "negative" = "forestgreen"))
       + labs(x= "chr 2", y = "Directionality Index", size=12 )
       + theme_bw()
       + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour="black", size=0.5),
               axis.text.x = element_text(colour="black", size=14),
               axis.text.y = element_text(colour="black", size=14))
)
p1
dev.off()






stop()
diff.df <- data.frame(A549_R1$header, A549_R1$start, abs(A549_R1$rawInsulationScore - T470A_R1$rawInsulationScore))
colnames(diff.df) = c("loc","start", "diff")


sorted.diff.df <- diff.df[order(diff.df$diff, decreasing=TRUE),]

A549_R1_ex <- A549_R1[grep("chr21", A549_R1$header),]

A549_R1_ex <- A549_R1_ex[A549_R1$start > 12000000 & A549_R1$start < 18000000,]

options(scipen=100)
png("A549_chr21_12Mb-18Mb.png", width=5000, height=1500, res=300)
plot(A549_R1_ex$midpoint, A549_R1_ex$rawInsulationScore, col="blue", pch=20, type="o", xlab="chr21",
     ylab="Insulation Score")
dev.off()


T470A_R1_ex <- T470A_R1[grep("chr21", T470A_R1$header),]

T470A_R1_ex <- T470A_R1_ex[T470A_R1$start > 12000000 & T470A_R1$start < 18000000,]

options(scipen=100)
png("T470A_chr21_12Mb-18Mb.png", width=5000, height=1500, res=300)
plot(T470A_R1_ex$midpoint, T470A_R1_ex$rawInsulationScore, col="blue", pch=20, type="o", xlab="chr21",
     ylab="Insulation Score")
dev.off()

A549_R1_ex <- A549_R1[grep("chr22", A549_R1$header),]

A549_R1_ex <- A549_R1_ex[A549_R1$start > 21000000 & A549_R1$start < 29000000,]

options(scipen=100)
png("A549_chr22_21Mb-29Mb.png", width=5000, height=1500, res=300)
plot(A549_R1_ex$midpoint, A549_R1_ex$rawInsulationScore, col="blue", pch=20, type="o", xlab="chr22",
     ylab="Insulation Score",ylim=c(10,35))
abline(v= 24460000, col="red", lty=2)
dev.off()


T470A_R1_ex <- T470A_R1[grep("chr22", T470A_R1$header),]

T470A_R1_ex <- T470A_R1_ex[T470A_R1$start > 21000000 & T470A_R1$start < 29000000,]

options(scipen=100)
png("T470A_chr22_21Mb-29Mb.png", width=5000, height=1500, res=300)
plot(T470A_R1_ex$midpoint, T470A_R1_ex$rawInsulationScore, col="blue", pch=20, type="o", xlab="chr22",
     ylab="Insulation Score", ylim=c(10,35))
abline(v= 24460000, col="red", lty=2)
dev.off()