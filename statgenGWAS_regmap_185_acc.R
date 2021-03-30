# loads phenotype data
load("rdata/phenotypes_regmap.RData")

# loads genotype data (snp)
load("rdata/inferred_SNP_matrix_regmap_185.RData")

# loads snp locations in the genome
load("rdata/inferred_snp_chromosomic_map_185.RData")

library(statgenGWAS)

# individuals for which we have snp data
acc <- intersect(regmap$ecotype, colnames(snp))


# chromosomic positions
colnames(map) <- c("chr", "pos")
head(map)

# phenotype data are the changes caused by elevated CO2
Y <- regmap[match(acc, regmap$ecotype),stringr::str_detect(colnames(regmap), "change")]
Y$genotype <- acc
Y <- Y[,c("genotype", colnames(Y)[colnames(Y) != "genotype"])]

# snps
X <- t(snp)

# kinship matrix : computes genetic relatedness between individuals
kin <- kinship(X, method = c("astle"))
heatmap(kin)


# object storing all info for gwas analysis with statgenGWAS
gData <- createGData(geno = X, map = map, pheno = Y, kin = kin)
summary(gData)

# running mixed models for all snps and getting their pvalues
gwas <- runSingleTraitGwas(gData = gData,
                           traits = c("N_change"),
                           GLSMethod = "single",
                           remlAlgo = "EMMA",
                           thrType = "fixed",
                           LODThr = 6)

# run for all elements at the same time (run independently,
# but the results are merged)
# gwas <- runSingleTraitGwas(gData = gData,
#                            GLSMethod = "single",
#                            remlAlgo = "EMMA",
#                            nCores = 6)

# get significant SNPs
t <- gwas$signSnp$Y[,c("pValue", "chr", "pos")]
t
#capture.output(t, file = "t50CVanRaden.csv")
write.table(t, "tZnAstle.csv", row.names=FALSE, sep=",")

# peut prendre du temps
plot(gwas, plotType = "manhattan", trait = "N_change", lod = 3)


# bed <- read.csv("data/TAIR10_GFF3_genes.gff", sep = '\t', h = F)[,c(1:5,9) ]
# colnames(bed) <- c("chr", "annot", "type", "start", "end", "name")
#save(bed, file = "rdata/tair10_annotation.Rdata")




# bed variable with gene coordinates
load("rdata/tair10_annotation.Rdata")

# functions that find which genes are around
get_overlapping_region <- function(chr, pos, bed, type = "all"){
  d <- bed[bed$chr == paste0("Chr",chr),]
  d <- bed[pos > bed$start & pos < bed$end, ]
  d$name <- stringr::str_split_fixed(stringr::str_split_fixed(d$name, "ID=", 2)[,2], ";", 2)[,1]
  if(type == "gene")
    d <- d[d$type == "gene",]
  d$snp_pos <- rep(pos, nrow(d))
  return(d)
}
# ne r?cup?re pas les snp proches des d?buts et fins de g?nes


t1 = read.table(file = "tZnAstle.csv")


describe_snps <- function(t, bed, type = "gene"){
  res <- get_overlapping_region(as.numeric(t[1, "chr"]),as.numeric(t[1, "pos"]), bed,
                                type = type)
  for(i in 2:nrow(t)){
    res <- rbind.data.frame(res, 
                            get_overlapping_region(as.numeric(t[i, "chr"]),
                                                   as.numeric(t[i, "pos"]), bed,
                                                   type = type))
  }
  return(res)
}


# result !!! now the genes found (column "name")
# can be searched on TAIR10 to see their functions, etc
genes <- describe_snps(t1, bed)
genes



#infos <- DIANE::get_gene_information(unique(genes$name), "Arabidopsis thaliana")
#View(infos)