# "Enrichment analysis of differentially expression genes"
# library(BiocManager)
# BiocManager::install('clusterProfiler')
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("org.Hs.eg.db")
rm(list=ls())

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(org.Hs.eg.db))

help_info <- c("ENSEMBL", "ENTREZID", "GENENAME", "REFSEQ", "SYMBOL", "UNIPROT")

# Parse parameters
parser <- OptionParser()
parser <- add_option(parser, c("--input"), action="store", default=NA, type='character',
                     help="input gene list")
parser <- add_option(parser, c("--output"), action="store", default=NA, type='character',
                     help="output csv file")
parser <- add_option(parser, c("--db"), action="store", default='Hg', type='character',
                     help="reference database")
parser <- add_option(parser, c("--type"), action="store", default="SYMBOL", type='character',
                     help=paste(help_info, collapse = '; '))
opt <- parse_args(parser)

# main point of program is here, do this whether or not "verbose" is set
if (is.na(opt$input) || is.na(opt$output)) {
    cat("Please specify --in/--out, refer to the manual for detailed instruction!\n",
        file=stderr())
    quit()
}
if (!(opt$db %in% c("Hs", "Mm"))) {
    cat("Please specify --in/--out, refer to the manual for detailed instruction!\n",
        file=stderr())
    quit() 
}

# Load gene list
db = paste(c("org", opt$db, "eg", "db"), collapse=".")
gene_list <- as.character(read.csv(opt$input, header = FALSE)$V1)
id_converted <- bitr(gene_list, fromType=opt$type, toType="ENTREZID", OrgDb=db)

if (opt$db == "Hs") {
    bp <- enrichGO(id_converted$ENTREZID, OrgDb = org.Hs.eg.db, ont="BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.05, keyType = 'ENTREZID')
    mf <- enrichGO(id_converted$ENTREZID, OrgDb = org.Hs.eg.db, ont="MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.05, keyType = 'ENTREZID')
    cc <- enrichGO(id_converted$ENTREZID, OrgDb = org.Hs.eg.db, ont="CC", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.05, keyType = 'ENTREZID')
} else if (opt$db == "Mm") {
    bp <- enrichGO(id_converted$ENTREZID, OrgDb = org.Mm.eg.db, ont="BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.05, keyType = 'ENTREZID')
    mf <- enrichGO(id_converted$ENTREZID, OrgDb = org.Mm.eg.db, ont="MF", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.05, keyType = 'ENTREZID')
    cc <- enrichGO(id_converted$ENTREZID, OrgDb = org.Mm.eg.db, ont="CC", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.05, keyType = 'ENTREZID')
} else {
    # Pass
}

bp <- data.frame(simplify(bp))
bp$ONTOLOGY = "BP"

mf <- data.frame(simplify(mf))
mf$ONTOLOGY = "MF"

cc <- data.frame(simplify(cc))
cc$ONTOLOGY = "CC"

res <- rbind(bp, mf, cc)

write.csv(res, file = opt$output)
