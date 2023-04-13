library(ape)
library(nlme)
library(RERconverge)
library(optparse)
source("glsMC.R")
source("permulation_functions.R")

option_list = list(make_option(c("-m", "--elementscoresfolder"), type="character", default=NULL, help="folder path to motif scores of all elements", metavar="character"),
		   make_option(c("-t", "--treepath"), type="character", default=NULL, help="file path to Newick neutral tree", metavar="character"),
		   make_option(c("-p", "--traitpath"), type="character", default=NULL, help="file (.RDS) path to trait values", metavar="character"),
		   make_option(c("-n", "--numperms"), type="integer", default=100, help="number of permulations", metavar="integer"),
		   make_option(c("-x", "--nulltraits"), type="character", default=NULL, help="file (.RDS) path to precomputed null traits list object", metavar="character"),
		   make_option(c("-o", "--outputfolder"), type="character", default=NULL, help="output folder", metavar="character"),
		   make_option(c("-s", "--allstatsfolder"), type="character", default=NULL, help="all statistics folder", metavar="character")
)


print('Parsing arguments')
opt_parser = OptionParser(option_list=option_list)
opt=parse_args(opt_parser)

elementscoresfolder = opt$elementscoresfolder
treepath = opt$treepath
traitpath = opt$traitpath
if (!is.null(opt$numperms)){
	numperms = opt$numperms
}
if (!is.null(opt$nulltraits)){
	nulltraitspath = opt$nulltraits
}
outputfolder = opt$outputfolder
allstatsfolder = opt$allstatsfolder

element_idx = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

all_elements = list.files(elementscoresfolder)


print('Loading inputs')
masterTree=ape::read.tree(file=treepath)
phenoVals = readRDS(traitpath)
nulltraits = readRDS(nulltraitspath)

thistree=masterTree # we may actually want to use a different tree

use.sp=intersect(thistree$tip.label, names(phenoVals)) #retain only the species we have a phenotype for
thistree=drop.tip(thistree, setdiff(thistree$tip.label, use.sp))
phenoVals=phenoVals[use.sp]

elementscorespath = all_elements[element_idx]
featureMat = readDataFast(paste0(elementscoresfolder, elementscorespath))

outputpath = paste0(outputfolder, sub("motif_scores_", "permulated_scores_", sub(".bed", ".RDS",elementscorespath)))
allstatspath = paste0(allstatsfolder, sub("motif_scores_", "all_statistics_", sub(".bed", ".RDS",elementscorespath)))

use.sp=intersect(names(phenoVals),rownames(featureMat)) #retain only the species we have a phenotype for that also have the promoter

print(paste("Using", length(use.sp), "species"))

#make the subset tree
thistree=drop.tip(thistree, setdiff(thistree$tip.label, use.sp))
featureMat=featureMat[use.sp,, drop=F]
phenoVals=phenoVals[use.sp]

print('Computing permulated scores')
coeff_obs = cor(featureMat, phenoVals, method='spearman')


# permulations
coeff_null = lapply(nulltraits, getCorrelationCoefficient, featureMat=featureMat, method='spearman')
coeff_null_df = as.data.frame(coeff_null)
colnames(coeff_null_df) = paste0('null', seq(1,ncol(coeff_null_df), 1))
# combine with coeff_obs
coeff_all_df = cbind(coeff_obs, coeff_null_df)

permulation_output = lapply(1:length(coeff_obs), computePermulatedScoresApplyV2, coeff_obs_all=coeff_obs, coeff_null_all=coeff_null_df)
permulation_output_df = as.data.frame(do.call(rbind,permulation_output))
rownames(permulation_output_df) = rownames(coeff_null_df)

print('Saving')
saveRDS(permulation_output_df, outputpath)
saveRDS(coeff_all_df, allstatspath)

