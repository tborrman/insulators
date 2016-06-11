
DI_cor <- function(rep1, rep2, chrom) {
  return(cor(rep1[rep1$chrom == paste("chr",chrom, sep=""),]$DI, rep2[rep2$chrom == paste("chr",chrom, sep=""),]$DI, use="complete.obs"))
  
}



suffix <- "__genome__C-40000-iced__chr_allchrom_dixonDI_w-2000000"


A549_R1 <-read.table(paste('ENCODE3-A549C-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
A549_R2 <-read.table(paste('ENCODE3-A549D-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
NCIH460A_R1 <-read.table(paste('ENCODE3-NCIH460A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
NCIH460B_R2 <-read.table(paste('ENCODE3-NCIH460B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
SKMEL5A_R1 <-read.table(paste('ENCODE3-SKMEL5A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
SKMEL5B_R2 <- read.table(paste('ENCODE3-SKMEL5B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
PANC1B_R1 <- read.table(paste('ENCODE3-PANC1B-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
PANC1C_R2 <- read.table(paste('ENCODE3-PANCIC-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
SKNDZA_R1 <- read.table(paste('ENCODE3-SKNDZA-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
SKNDZB_R2 <- read.table(paste('ENCODE3-SKNDZB-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
CAKi2A_R1 <- read.table(paste('ENCODE3-Caki2A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
CAKi2B_R2 <- read.table(paste('ENCODE3-CAK12B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
RPMI7951C_R1 <- read.table(paste('ENCODE3-RPMI7951C-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
RPMI7951D_R2 <- read.table(paste('ENCODE3-RPMI7951D-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
SKNMCC_R1 <- read.table(paste('ENCODE3-SKNMCC-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
SKNMCD_R2 <- read.table(paste('ENCODE3-SKNMCD-HindII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
G401A_R1 <- read.table(paste('ENCODE3-G401A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
G401B_R2 <- read.table(paste('ENCODE3-G401B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
SJCRH30A_R1 <- read.table(paste('ENCODE3-SJCRH30A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
SJCRH30B_R2 <- read.table(paste('ENCODE3-SJCRH30B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
T470A_R1 <- read.table(paste('ENCODE3-T470A-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
T470B_R2 <- read.table(paste('ENCODE3-T470B-HindIII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")
LNCaPC_R1 <- read.table(paste('ENCODE3-LNCaPC-HindIII-R1__hg19', suffix, sep=""),  header=TRUE, sep="\t")
LNCaP_R2 <- read.table(paste('ENCODE3-LNCaP-HindII-R2__hg19', suffix, sep=""),  header=TRUE, sep="\t")

df <- data.frame(matrix(ncol=23, nrow=12))
chroms <- c(seq(1:22), "X")
colnames(df) <- paste(rep("chr", 23), chroms, sep="")
for (chrom in chroms) {
  corrs <- c()
  corrs <- c(DI_cor(A549_R1, A549_R2, chrom),
             DI_cor(NCIH460A_R1, NCIH460B_R2, chrom),
             DI_cor(SKMEL5A_R1, SKMEL5B_R2, chrom),
             DI_cor(PANC1B_R1, PANC1C_R2, chrom),
             DI_cor(SKNDZA_R1, SKNDZB_R2, chrom),
             DI_cor(CAKi2A_R1, CAKi2B_R2, chrom),
             DI_cor(RPMI7951C_R1, RPMI7951D_R2, chrom),
             DI_cor(SKNMCC_R1, SKNMCD_R2, chrom),
             DI_cor(G401A_R1, G401B_R2, chrom),
             DI_cor(SJCRH30A_R1, SJCRH30B_R2, chrom),
             DI_cor(T470A_R1, T470B_R2, chrom),
             DI_cor(LNCaPC_R1, LNCaP_R2, chrom)
             )
  df[paste("chr", chrom, sep="")] <- corrs
}
png("cor_chrom_boxplots.png", width= 3000, height = 1500, res = 300)
  boxplot(df, col=rainbow(23), ylab = "DI correlation between replicates")
dev.off()


