C_df = read.csv("t5C10best.csv")
Mn_df = read.table("Mn_thr5_top10.csv")
Cu_df = read.table("Cu_thr5_top10.csv")
Zn_df = read.table("Zn_thr5_top10.csv")
Fe_df = read.table("Fe_thr5_top10.csv")
N_df = read.csv("t5N10best.csv")
Na_df = read.csv("t5Na10best.csv")
Mg_df = read.table("Mg_thr5_top10.csv")
Comp1_df = read.csv("t10Comp1.csv")


library(dplyr)


# Add an id col and reorder columns
C_df["SNP"] = paste(C_df$pos,"_",C_df$chr)
C_df = C_df[,c("SNP","chr","pos")]
C_df["C"] = 1

Mn_df["SNP"] = paste(Mn_df$pos,"_",Mn_df$chr)
Mn_df["Mn"] = 1
Mn_df = Mn_df[,c("SNP","chr","pos","Mn")]

Cu_df["SNP"] = paste(Cu_df$pos,"_",Cu_df$chr)
Cu_df["Cu"] = 1
Cu_df = Cu_df[,c("SNP","chr","pos", "Cu")]

Zn_df["SNP"] = paste(Zn_df$pos,"_",Zn_df$chr)
Zn_df["Zn"] = 1
Zn_df = Zn_df[,c("SNP","chr","pos", "Zn")]

Fe_df["SNP"] = paste(Fe_df$pos,"_",Fe_df$chr)
Fe_df["Fe"] = 1
Fe_df = Fe_df[,c("SNP","chr","pos", "Fe")]

N_df["SNP"] = paste(N_df$pos,"_",N_df$chr)
N_df["N"] = 1
N_df = N_df[,c("SNP","chr","pos", "N")]

Na_df["SNP"] = paste(Na_df$pos,"_",Na_df$chr)
Na_df["Na"] = 1
Na_df = Na_df[,c("SNP","chr","pos", "Na")]

Mg_df["SNP"] = paste(Mg_df$pos,"_",Mg_df$chr)
Mg_df["Mg"] = 1
Mg_df = Mg_df[,c("SNP","chr","pos", "Mg")]

Comp1_df["SNP"] = paste(Comp1_df$pos,"_",Comp1_df$chr)
Comp1_df["Comp1"] = 1
Comp1_df = Comp1_df[,c("SNP","chr","pos", "Comp1")]

total_SNP = full_join(C_df,Mn_df, by = c("SNP","chr","pos"))
total_SNP = full_join(total_SNP, Cu_df, by = c("SNP","chr","pos"))
total_SNP = full_join(total_SNP, Zn_df, by = c("SNP","chr","pos"))
total_SNP = full_join(total_SNP, Fe_df, by = c("SNP","chr","pos"))
total_SNP = full_join(total_SNP, N_df, by = c("SNP","chr","pos"))
total_SNP = full_join(total_SNP, Na_df, by = c("SNP","chr","pos"))
total_SNP = full_join(total_SNP, Mg_df, by = c("SNP","chr","pos"))
total_SNP = full_join(total_SNP, Comp1_df, by = c("SNP","chr","pos"))

write.table(total_SNP,"SNP_allelements_thr5.csv")

any(duplicated(total_SNP$SNP))
