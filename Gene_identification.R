library(dplyr)


get_region <- function(chromosome,position,bed,type="all", distance = 5000){
  d <- filter(bed, bed$chr == paste0("Chr",chromosome) & position > bed$start - distance & position < bed$end +distance)
  d$name <- stringr::str_split_fixed(stringr::str_split_fixed(d$name, "ID=", 2)[,2], ";", 2)[,1]
  if(type == "gene")
    d <- d[d$type == "gene",]
  d$snp_pos <- rep(position, nrow(d))
  return(d)
}

# ne r?cup?re pas les snp proches des d?buts et fins de g?nes

describe_snps <- function(t, bed, type = "gene", distance = 5000){
  res <- get_region(as.numeric(t[1, "chr"]),as.numeric(t[1, "pos"]), bed,
                                type = type, distance = distance)
  for(i in 2:nrow(t)){
    res <- rbind.data.frame(res, 
                            get_region(as.numeric(t[i, "chr"]),
                                                   as.numeric(t[i, "pos"]), bed,
                                                   type = type, distance = distance))
  }
  return(res)
}


# result !!! now the genes found (column "name")
# can be searched on TAIR10 to see their functions, etc
genesCu <- describe_snps(Cu_df, bed)
genesCu
genesCu <- na.omit(genesCu)
genesCu["Cu"] = 1
genesCu_dedup <- subset(genesCu, select = -c(snp_pos))
genesCu_dedup <- genesCu_dedup[!duplicated(genesCu_dedup), ]

genesC <- describe_snps(C_df,bed)
genesC <- na.omit(genesC)
genesC["C"] = 1
genesC_dedup <- subset(genesC, select = -c(snp_pos))
genesC_dedup <- genesC_dedup[!duplicated(genesC_dedup), ]

genesFe <- describe_snps(Fe_df,bed)
genesFe <- na.omit(genesFe)
genesFe["Fe"] = 1
genesFe_dedup <- subset(genesFe, select = -c(snp_pos))
genesFe_dedup <- genesFe_dedup[!duplicated(genesFe_dedup), ]

genesZn <- describe_snps(Zn_df,bed)
genesZn <- na.omit(genesZn)
genesZn["Zn"] = 1
genesZn_dedup <- subset(genesZn, select = -c(snp_pos))
genesZn_dedup <- genesZn_dedup[!duplicated(genesZn_dedup), ]


genesMn <- describe_snps(Mn_df,bed)
genesMn <- na.omit(genesMn)
genesMn["Mn"] = 1
genesMn_dedup <- subset(genesMn, select = -c(snp_pos))
genesMn_dedup <- genesMn_dedup[!duplicated(genesMn_dedup), ]

genesN <- describe_snps(N_df,bed)
genesN["N"] = 1
genesN <- na.omit(genesN)
genesN_dedup <- subset(genesN, select = -c(snp_pos))
genesN_dedup <- genesN_dedup[!duplicated(genesN_dedup), ]

genesMg <- describe_snps(Mg_df,bed)
genesMg <- na.omit(genesMg)
genesMg["Mg"] = 1
genesMg_dedup <- subset(genesMg, select = -c(snp_pos))
genesMg_dedup <- genesMg_dedup[!duplicated(genesMg_dedup), ]

genesNa <- describe_snps(Na_df,bed)
genesNa <- na.omit(genesNa)
genesNa["Na"] = 1
genesNa_dedup <- subset(genesNa, select = -c(snp_pos))
genesNa_dedup <- genesNa_dedup[!duplicated(genesNa_dedup), ]

genesComp1 <- describe_snps(Comp1_df,bed)
genesComp1 <- na.omit(genesComp1)
genesComp1["Comp1"] = 1
genesComp1_dedup <- subset(genesComp1, select = -c(snp_pos))
genesComp1_dedup <- genesComp1_dedup[!duplicated(genesComp1_dedup), ]


total_gene <- full_join(genesFe_dedup ,genesCu_dedup, by = c("chr","annot","type","start","end","name"))
total_gene <- full_join(total_gene ,genesZn_dedup, by = c("chr","annot","type","start","end","name"))
total_gene <- full_join(total_gene ,genesMn_dedup, by = c("chr","annot","type","start","end","name"))
total_gene <- full_join(total_gene ,genesN_dedup, by = c("chr","annot","type","start","end","name"))
total_gene <- full_join(total_gene ,genesC_dedup, by = c("chr","annot","type","start","end","name"))
total_gene <- full_join(total_gene ,genesMg_dedup, by = c("chr","annot","type","start","end","name"))
total_gene <- full_join(total_gene ,genesNa_dedup, by = c("chr","annot","type","start","end","name"))
total_gene <- full_join(total_gene ,genesComp1_dedup, by = c("chr","annot","type","start","end","name"))

any(duplicated(total_gene$name))

write.table(total_SNP,"genes_5kb.csv")

