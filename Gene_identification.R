library(dplyr)


get_region <- function(chromosome,position,bed,type="all", distance = 0){
  d <- filter(bed, bed$chr == paste0("Chr",chromosome) & position > bed$start - distance & position < bed$end +distance)
  d$name <- stringr::str_split_fixed(stringr::str_split_fixed(d$name, "ID=", 2)[,2], ";", 2)[,1]
  if(type == "gene")
    d <- d[d$type == "gene",]
  d$snp_pos <- rep(position, nrow(d))
  return(d)
}

# ne r?cup?re pas les snp proches des d?buts et fins de g?nes

describe_snps <- function(t, bed, type = "gene", distance = 0){
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

genesC <- describe_snps(C_df,bed)
genesC <- na.omit(genesC)
genesC["C"] = 1

genesFe <- describe_snps(Fe_df,bed)
genesFe <- na.omit(genesFe)
genesFe["Fe"] = 1

genesZn <- describe_snps(Zn_df,bed)
genesZn <- na.omit(genesZn)
genesZn["Zn"] = 1


genesMn <- describe_snps(Mn_df,bed)
genesMn <- na.omit(genesMn)
genesMn["Mn"] = 1

genesN <- describe_snps(N_df,bed)
genesN["N"] = 1


total_gene <- full_join(genesFe ,genesCu, by = c("chr","annot","type","start","end","name"))
total_gene <- full_join(total_gene ,genesZn, by = c("chr","annot","type","start","end","name"))
total_gene <- full_join(total_gene ,genesMn, by = c("chr","annot","type","start","end","name"))
total_gene <- full_join(total_gene ,genesN, by = c("chr","annot","type","start","end","name"))

write.table(total_SNP,"genes_almost_all_elements.csv")
