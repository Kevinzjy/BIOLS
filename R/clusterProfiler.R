# "Enrichment analysis of differentially expression genes"
# library(BiocManager)
# BiocManager::install('clusterProfiler')
# BiocManager::install("org.Hs.eg.db")
rm(list=ls())

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))

setwd('/home/zhangjy/04.XueY_Autophagy')
help_info <- keytypes(org.Hs.eg.db)

# Parse parameters
parser <- OptionParser()
parser <- add_option(parser, c("--input"), action="store", default=NA, type='character',
                     help="input gene list")
parser <- add_option(parser, c("--output"), action="store", default=NA, type='character',
                     help="output csv file")
parser <- add_option(parser, c("--type"), action="store", default="SYMBOL", type='character',
                     help=paste(help_info, collapse = '; '))
opt <- parse_args(parser)

# main point of program is here, do this whether or not "verbose" is set
if (is.na(opt$input) || is.na(opt$output)) {
    cat("Please specify --in/--out, refer to the manual for detailed instruction!\n",
        file=stderr())
    quit()
}

# Load gene list
gene_list <- as.character(read.csv(opt$input, header = FALSE)$V1)
id_converted <- bitr(gene_list, fromType=opt$type, toType="ENTREZID", OrgDb="org.Hs.eg.db")
go <- enrichGO(id_converted$ENTREZID, OrgDb = org.Hs.eg.db, ont="ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05,
               qvalueCutoff = 0.05, keyType = 'ENTREZID')
write.csv(go, file = opt$output)
