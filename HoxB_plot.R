A549_df <- read.table("ENCODE3-A549C-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000", header=TRUE, sep="\t")
Caki2_df <- read.table("ENCODE3-Caki2A-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000", header=TRUE, sep="\t")
G401_df <- read.table("ENCODE3-G401A-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000", header=TRUE, sep="\t")
LNCaP_df <- read.table("ENCODE3-LNCaPC-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000", header=TRUE, sep="\t")
NCIH460_df <- read.table("ENCODE3-NCIH460A-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000", header=TRUE, sep="\t")
PANC1_df <- read.table("ENCODE3-PANC1B-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000",  header=TRUE, sep="\t")
RPMI7951_df <- read.table("ENCODE3-RPMI7951C-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000", header=TRUE, sep="\t")
SKMEL5_df <- read.table("ENCODE3-SKMEL5A-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000", header=TRUE, sep="\t")
SKNDZ_df <- read.table("ENCODE3-SKNDZA-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000", header=TRUE, sep="\t")
SKNMC_df <- read.table("ENCODE3-SKNMCC-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000", header=TRUE, sep="\t")
T47D_df <- read.table("ENCODE3-T470A-HindIII-R1__hg19__genome__C-40000-iced__chr17.scaled-1000000_dixonDI_w-2000000", header=TRUE, sep="\t")


HoxB <- data.frame()
HoxB[1,1] <- "A549"
HoxB[1,2] <- A549_df[A549_df$start == 46640001, "DI"]
HoxB[1,3] <- A549_df[A549_df$start == 46800001, "DI"]
HoxB[2,1] <- "Caki2"
HoxB[2,2] <- Caki2_df[Caki2_df$start == 46640001, "DI"]
HoxB[2,3] <- Caki2_df[Caki2_df$start == 46800001, "DI"]
HoxB[3,1] <- "G401"
HoxB[3,2] <- G401_df[G401_df$start == 46640001, "DI"]
HoxB[3,3] <- G401_df[G401_df$start == 46800001, "DI"]
HoxB[4,1] <- "LNCaP"
HoxB[4,2] <- LNCaP_df[LNCaP_df$start == 46640001, "DI"]
HoxB[4,3] <- LNCaP_df[LNCaP_df$start == 46800001, "DI"]
HoxB[5,1] <- "NCIH460"
HoxB[5,2] <- NCIH460_df[NCIH460_df$start == 46640001, "DI"]
HoxB[5,3] <- NCIH460_df[NCIH460_df$start == 46800001, "DI"]
HoxB[6,1] <- "PANC1"
HoxB[6,2] <- PANC1_df[PANC1_df$start == 46640001, "DI"]
HoxB[6,3] <- PANC1_df[PANC1_df$start == 46800001, "DI"]
HoxB[7,1] <- "RPMI7951"
HoxB[7,2] <- RPMI7951_df[RPMI7951_df$start == 46640001, "DI"]
HoxB[7,3] <- RPMI7951_df[RPMI7951_df$start == 46800001, "DI"]
HoxB[8,1] <- "SKMEL5"
HoxB[8,2] <- SKMEL5_df[SKMEL5_df$start == 46640001, "DI"]
HoxB[8,3] <- SKMEL5_df[SKMEL5_df$start == 46800001, "DI"]
HoxB[9,1] <- "SKNDZ"
HoxB[9,2] <- SKNDZ_df[SKNDZ_df$start == 46640001, "DI"]
HoxB[9,3] <- SKNDZ_df[SKNDZ_df$start == 46800001, "DI"]
HoxB[10,1] <- "SKNMC"
HoxB[10,2] <- SKNMC_df[SKNMC_df$start == 46640001, "DI"]
HoxB[10,3] <- SKNMC_df[SKNMC_df$start == 46800001, "DI"]
HoxB[11,1] <- "T47D"
HoxB[11,2] <- T47D_df[T47D_df$start == 46640001, "DI"]
HoxB[11,3] <- T47D_df[T47D_df$start == 46800001, "DI"]

colnames(HoxB) <- c("cell", "B3B6", "B13")

plot(HoxB$B3B6, HoxB$B13)
text(HoxB$B3B6, HoxB$B13, HoxB$cell)

