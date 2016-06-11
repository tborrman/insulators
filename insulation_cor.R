
samples = c('ENCODE3-A549C-HindIII-R1__hg19', 
           'ENCODE3-NCIH460B-HindIII-R2__hg19',
           'ENCODE3-SKMEL5B-HindIII-R2__hg19',		
           'ENCODE3-A549D-HindIII-R2__hg19',         
           'ENCODE3-PANC1B-HindIII-R1__hg19',         
           'ENCODE3-SKNDZA-HindIII-R1__hg19',
           'ENCODE3-CAK12B-HindIII-R2__hg19',        
           'ENCODE3-PANCIC-HindIII-R2__hg19',         
           'ENCODE3-SKNDZB-HindIII-R2__hg19',
           'ENCODE3-Caki2A-HindIII-R1__hg19',        
           'ENCODE3-RPMI7951C-HindIII-R1__hg19',      
           'ENCODE3-SKNMCC-HindIII-R1__hg19',
           'ENCODE3-G401A-HindIII-R1__hg19',         
           'ENCODE3-RPMI7951D-HindIII-R2__hg19',      
           'ENCODE3-SKNMCD-HindII-R2__hg19',
           'ENCODE3-G401B-HindIII-R2__hg19',         
           'ENCODE3-SJCRH30A-HindIII-R1__hg19',       
           'ENCODE3-T470A-HindIII-R1__hg19',
           'ENCODE3-LNCaP-HindII-R2__hg19',          
           'ENCODE3-SJCRH30B-HindIII-R2__hg19',       
           'ENCODE3-T470B-HindIII-R2__hg19',
           'ENCODE3-NCIH460A-HindIII-R1__hg19',      
           'ENCODE3-SKMEL5A-HindIII-R1__hg19',
           'ENCODE3-LNCaPC-HindIII-R1__hg19')

suffix <- "__genome__C-40000-iced__chr_allchrom_dixonDI_w-2000000"


A549_R1 <-read.table(paste('ENCODE3-A549C-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
A549_R2 <-read.table(paste('ENCODE3-A549D-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
NCIH460A_R1 <-read.table(paste('ENCODE3-NCIH460A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
NCIH460B_R2 <-read.table(paste('ENCODE3-NCIH460B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
SKMEL5A_R1 <-read.table(paste('ENCODE3-SKMEL5A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
SKMEL5B_R2 <- read.table(paste('ENCODE3-SKMEL5B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
PANC1B_R1 <- read.table(paste('ENCODE3-PANC1B-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
PANC1C_R2 <- read.table(paste('ENCODE3-PANCIC-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
SKNDZA_R1 <- read.table(paste('ENCODE3-SKNDZA-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
SKNDZB_R2 <- read.table(paste('ENCODE3-SKNDZB-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
CAKi2A_R1 <- read.table(paste('ENCODE3-Caki2A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
CAKi2B_R2 <- read.table(paste('ENCODE3-CAK12B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
RPMI7951C_R1 <- read.table(paste('ENCODE3-RPMI7951C-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
RPMI7951D_R2 <- read.table(paste('ENCODE3-RPMI7951D-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
SKNMCC_R1 <- read.table(paste('ENCODE3-SKNMCC-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
SKNMCD_R2 <- read.table(paste('ENCODE3-SKNMCD-HindII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
G401A_R1 <- read.table(paste('ENCODE3-G401A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
G401B_R2 <- read.table(paste('ENCODE3-G401B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
SJCRH30A_R1 <- read.table(paste('ENCODE3-SJCRH30A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
SJCRH30B_R2 <- read.table(paste('ENCODE3-SJCRH30B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
T470A_R1 <- read.table(paste('ENCODE3-T470A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
T470B_R2 <- read.table(paste('ENCODE3-T470B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
LNCaPC_R1 <- read.table(paste('ENCODE3-LNCaPC-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI
LNCaP_R2 <- read.table(paste('ENCODE3-LNCaP-HindII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")$DI

png("rep_cors.png", height=2500, width=2000, res=300)
par(mfrow=c(4,3), mar=c(4,4,1,1) + 0.5)
plot(A549_R1, A549_R2, pch=21, col="black", bg="blue")
text(5000,-7000, labels= paste("r =", round(cor(A549_R1, A549_R2, use="complete.obs"), 2)))
plot(NCIH460A_R1, NCIH460B_R2, pch=21, col="black", bg="blue")
text(5000,-7000, labels= paste("r =", round(cor(NCIH460A_R1, NCIH460B_R2, use="complete.obs"), 2)))
plot(SKMEL5A_R1, SKMEL5B_R2, pch=21, col="black", bg="blue")
text(5000,-15000, labels= paste("r =", round(cor(SKMEL5A_R1, SKMEL5B_R2, use="complete.obs"), 2)))
plot(PANC1B_R1, PANC1C_R2, pch=21, col="black", bg="blue")
text(5000,-7000, labels= paste("r =", round(cor(PANC1B_R1, PANC1C_R2, use="complete.obs"), 2)))
plot(SKNDZA_R1, SKNDZB_R2, pch=21, col="black", bg="blue")
text(5000,-7000, labels= paste("r =", round(cor(SKNDZA_R1, SKNDZB_R2, use="complete.obs"), 2)))
plot(CAKi2A_R1, CAKi2B_R2, pch=21, col="black", bg="blue")
text(5000,-7000, labels= paste("r =", round(cor(CAKi2A_R1, CAKi2B_R2, use="complete.obs"), 2)))
plot(RPMI7951C_R1, RPMI7951D_R2, pch=21, col="black", bg="blue")
text(5000,-7000, labels= paste("r =", round(cor(RPMI7951C_R1, RPMI7951D_R2, use="complete.obs"), 2)))
plot(SKNMCC_R1, SKNMCD_R2, pch=21, col="black", bg="blue")
text(5000,-4000, labels= paste("r =", round(cor(SKNMCC_R1, SKNMCD_R2, use="complete.obs"), 2)))
plot(G401A_R1, G401B_R2, pch=21, col="black", bg="blue")
text(5000,-7000, labels= paste("r =", round(cor(G401A_R1, G401B_R2, use="complete.obs"), 2)))
plot(SJCRH30A_R1, SJCRH30B_R2, pch=21, col="black", bg="blue")
text(10000,-10000, labels= paste("r =", round(cor(SJCRH30A_R1, SJCRH30B_R2, use="complete.obs"), 2)))
plot(T470A_R1, T470B_R2, pch=21, col="black", bg="blue")
text(5000,-7000, labels= paste("r =", round(cor(T470A_R1, T470B_R2, use="complete.obs"), 2)))
plot(LNCaPC_R1, LNCaP_R2, pch=21, col="black", bg="blue")
text(5000,-5000, labels= paste("r =", round(cor(LNCaPC_R1, LNCaP_R2, use="complete.obs"), 2)))
dev.off()

df <- cbind(A549_R1,NCIH460A_R1,SKMEL5A_R1,PANC1B_R1,SKNDZA_R1, CAKi2A_R1, RPMI7951C_R1, SKNMCC_R1, G401A_R1, SJCRH30A_R1, T470A_R1, LNCaPC_R1) 


cor_matrix <- cor(df, use="pairwise.complete.obs");
# png("rep1_cors.png", height=3000, width=3000, res=300)
# pairs(df)
# dev.off()

library(gplots)


png("rep1_cors_heatmap.png", height=3000, width=3000, res=300)
par(oma=c(3,3,3,3) + 0.5)
heatmap.2(cor_matrix, trace="none")
dev.off()

