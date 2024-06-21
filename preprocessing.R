# This script was provided by Phil McCown, describing the data composition and annotation procedure within the NEPTUNE consortium
# Comments by Phil McCown
# 
# Update by Keon Arbabi (Maze, Consultant) Jan 2024
# 
# > added code to download relevant files for this analysis from Maze AWS
# > added roundToInt = TRUE to SoupX::adjustCounts() to prevent the count being converted to floats 
# > removed dropping of doublets, instead added `IsDoublet` and `pANN` to metadata
# > changed threshold for what percentage of total cells are assumed to be doublets from 7.5% to 3%, this is rather arbritary but 
# >> adjusts for preceding filtering out of high `percent.mt` and `nFeature_RNA` cells, which are usually doublets, and the fact that 
# >> current fluidics technology is better than 6 years ago when 7.5% doublets was assumed 

Download data ###########################################################################################################################
# Make the 10x data on AWS available locally
library(aws.s3)
# aws configure sso-session
Sys.setenv(AWS_PROFILE = "mazetx-sso-data-analyst-169945233738")

# Get required files
s3_path = 's3://maze-data-deposit-vendors/E6A4BDFF-7EDD-458D-BD17-9658D8CD0CA4/neptune-consortium/'
folders = paste0(s3_path, gsub('^(.+) ', '', system(paste0('aws s3 ls ', s3_path), intern = TRUE)))
folders = grep("-EO/", folders, value = TRUE)

data_paths = c()

for (folder in folders) {
  analysis_folder = grep("10x_analysis", gsub('^(.+) ', '', system(paste0('aws s3 ls ', folder), intern = TRUE)), value = TRUE)
  if (length(analysis_folder) > 0) {
    analysis_subfolders = gsub('^(.+) ', '', system(paste0('aws s3 ls ', folder, analysis_folder), intern = TRUE))
    for (subfolder in analysis_subfolders) {
      data_paths = c(data_paths, paste0(folder, analysis_folder, subfolder, "analysis/"))
      data_paths = c(data_paths, paste0(folder, analysis_folder, subfolder, "molecule_info.h5"))
      raw_matrix_path = paste0(folder, analysis_folder, subfolder, "raw_feature_bc_matrix/")
      data_paths = c(data_paths,
                     paste0(raw_matrix_path, "barcodes.tsv.gz"),
                     paste0(raw_matrix_path, "features.tsv.gz"),
                     paste0(raw_matrix_path, "matrix.mtx.gz"))
      filtered_matrix_path = paste0(folder, analysis_folder, subfolder, "filtered_feature_bc_matrix/")
      data_paths = c(data_paths,
                     paste0(filtered_matrix_path, "barcodes.tsv.gz"),
                     paste0(filtered_matrix_path, "features.tsv.gz"),
                     paste0(filtered_matrix_path, "matrix.mtx.gz"))
    }
  }
}
# Pull down files from AWS
system("cd ~/workspace/MZ-167-NEPTUNE-scratch")
for (data in data_paths) {
  command = paste("maze analysis input assign", data)
  system(command)
}
# Extract each file
data_paths_local = gsub("s3://", "/data/", data_paths)
for (data in data_paths_local) {
  if (grepl("\\.gz$", data)) {
    gunzip_command = paste("gunzip -k", data)
    system(gunzip_command)
  }
}

# Pre-process each sample #############################################################################################################

# These packages are necessary to run the base Seurat pipeline
library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(umap)

# Thankfully, the only package you need to run SoupX is SoupX. Make sure you have the most recent version installed,
# as there was an issue that I worked out with the package developer regarding the new cellranger file directory formatting
library(SoupX)

#These packages are needed to run DoubletFinder; some of these need to be installed via github
library(fields)
library(Matrix)
library(KernSmooth)
library(ROCR)
library(parallel)
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

setwd("/data/maze-data-deposit-vendors/E6A4BDFF-7EDD-458D-BD17-9658D8CD0CA4/neptune-consortium")

# 6801-EO-1 sample ####################################################################################################################
# I'll give a more lengthy commentary here for this sample, but the same general thoughts are incorporated with the rest of the samples.
# The order of my QC processing goes SoupX --> Seurat --> DoubletFinder. For the archived samples, I always use DoubletFinder as a QC step. 
# For relatively fresh samples and/or scRNA-seq, DoubletFinder isn't "as" necessary.
# SoupX to start
sc=load10X("6801-EO/10x_analysis_6801-EO/Sample_6801-EO-1-GEX_TGCGCGGTTT-TTTATCCTTG")
# The next three steps are allowing SoupX to run on autopilot. Thankfully, the samples that were selected didn't require a bit more 
# manual control to get high quality nuclei.
sc = autoEstCont(sc)
out = adjustCounts(sc, roundToInt = TRUE)
Patient68011 = CreateSeuratObject(out)
# The above step created a Seurat object from the SoupX-cleaned data. Note, SoupX does NOT remove nuclei 
# it only scrubs transcripts that are determined to be ambient contamination.
# Meta data loading...
Patient68011$EdgarID = "Patient68011"
Patient68011$ID = "Maze1"
Patient68011$Cohort = "FSGS"
Patient68011$Sex = "Male"
Patient68011$Age = "15"
Patient68011$APOL_Allele_Number = "2"
Patient68011$APOL_Alleles = "G2G2"
Patient68011$N264K = "No"
Patient68011$eGFR_Bx = "101.727"
Patient68011$UPCR_Bx = "7.6627"
# The next line is from base Seurat. This identifies the mitochondrial-encoded genes in the nuclei.
Patient68011$percent.mt = PercentageFeatureSet(Patient68011, pattern="^MT-")
# The next line actually removes the nuclei that have a lot of mitochondrial contamination.
# For nuclei from archived samples, I try to be a bit more generous with the mitochondrial contamination acceptance. 
# However, the features criteria are pretty standard and further minimize doublets (anything over 5k tends to be a doublet),
# while also removing cells that are possibly dead/dying (anything under 300 tends to be dead/dying, but I set it to 500 due to past issues with archived samples). 
Patient68011 = subset(Patient68011, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
# The next few lines are pretty normal for Seurat. 
# The parameters that I have set here are what have (historically) been reliable at ensuring good resolution with these data.
Patient68011 = NormalizeData(Patient68011, normalization.method="LogNormalize",scale.factor = 10000)
Patient68011 = FindVariableFeatures(Patient68011, selection.method="vst",nfeatures= 2000)
Patient68011 = ScaleData(Patient68011)
Patient68011 = RunPCA(Patient68011, features = VariableFeatures(object = Patient68011))
Patient68011 = FindNeighbors(Patient68011, dims = 1:30)
Patient68011 = FindClusters(Patient68011, resolution = 0.6)
Patient68011 = RunUMAP(Patient68011, dims = 1:30)
# The above is needed to run the data to a point where DoubletFinder can take over. 
# The remaining steps for this sample identify doublets and then I remove them with the subset command below.
# I also clean up the resulting matrix with the last two lines in this sample.
# I set the PCs = 1:10 to keep things simple from a computational processor resource perspective, 
# while also keeping in mind that doublets tend to appear in the larger clusters, hence the 1 through 10 variable usage. 
# The PCs are counted by Linux, not R, hence why 1 is selected as the first PC and not 0. 
# The next few lines are a bit of variable optimization based on the default DoubletFinder parameters so that you can copy/paste the script (up to a point).
sweep.res.list_Patient68011 = paramSweep(Patient68011, PCs = 1:10, sct = FALSE)
sweep.stats_Patient68011 = summarizeSweep(sweep.res.list_Patient68011, GT = FALSE)
bcmvn_Patient68011 = find.pK(sweep.stats_Patient68011)
pK_Patient68011 = bcmvn_Patient68011 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68011 = as.numeric(as.character(pK_Patient68011[[1]]))
annotations_Patient68011 = Patient68011@meta.data$seurat_clusters
homotypic.prop_Patient68011 = modelHomotypic(annotations_Patient68011)
nExp_poi_Patient68011 = round(0.05*nrow(Patient68011@meta.data))
nExp_poi_adj_Patient68011 = round(nExp_poi_Patient68011*(1-homotypic.prop_Patient68011))
Patient68011 = doubletFinder(Patient68011, PCs = 1:10, pN = 0.25, pK = pK_Patient68011, nExp = nExp_poi_adj_Patient68011, reuse.pANN = FALSE, sct = FALSE)
Patient68011@meta.data$IsDoublet = Patient68011@meta.data[, grep("DF.classifications", names(Patient68011@meta.data), value = TRUE)] != "Singlet"
Patient68011@meta.data$ProbDoublet = Patient68011@meta.data[, grep("pANN", names(Patient68011@meta.data), value = TRUE)] 
Patient68011@meta.data = Patient68011@meta.data[, !grepl("DF.classifications|pANN", names(Patient68011@meta.data))]

# # I ran table(Patient68011@meta.data$DF.classifications_0.25 (tabbed)) to determine the number of possible doublets that were identified.
# # I always proceed in removing doublets, but you don't necessarily have to.
# # Patient68011 = subset(x=Patient68011, subset = DF.classifications_0.25_0.005_709 == "Singlet")
# # The last two lines here are helpful in removing extraneous matrix columns that no longer have any relevance to the rest of the sample's processing in Seurat. 
# # Removal of these two columns now minimizes downstream matrix size, which doesn't matter now, but not having 2(N of samples) less of defunct matrices in the integrated object really frees up a lot of memory later on.
# # Rinse and repeat for the remaining samples. I won't add any additional comments for the rest of the individual sample processing.

# 6801-EO-2 sample ####################################################################################################################
sc1=load10X("6801-EO/10x_analysis_6801-EO/Sample_6801-EO-2-GEX_TGCAATGTTC-TTCGACAAGC")
sc1 = autoEstCont(sc1)
out1 = adjustCounts(sc1, roundToInt = TRUE)
Patient68012 = CreateSeuratObject(out1)
Patient68012$EdgarID = "Patient68012"
Patient68012$ID = "Maze2"
Patient68012$Cohort = "FSGS"
Patient68012$Sex = "Male"
Patient68012$Age = "11"
Patient68012$APOL_Allele_Number = "1"
Patient68012$APOL_Alleles = "G0G1"
Patient68012$N264K = "No"
Patient68012$eGFR_Bx = "84.839"
Patient68012$UPCR_Bx = "2.5058"
Patient68012$percent.mt = PercentageFeatureSet(Patient68012, pattern="^MT-")
Patient68012 = subset(Patient68012, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68012 = NormalizeData(Patient68012, normalization.method="LogNormalize",scale.factor = 10000)
Patient68012 = FindVariableFeatures(Patient68012, selection.method="vst",nfeatures= 2000)
Patient68012 = ScaleData(Patient68012)
Patient68012 = RunPCA(Patient68012, features = VariableFeatures(object = Patient68012))
Patient68012 = FindNeighbors(Patient68012, dims = 1:30)
Patient68012 = FindClusters(Patient68012, resolution = 0.6)
Patient68012 = RunUMAP(Patient68012, dims = 1:30)
sweep.res.list_Patient68012 = paramSweep(Patient68012, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68012 = summarizeSweep(sweep.res.list_Patient68012, GT = FALSE)
bcmvn_Patient68012 = find.pK(sweep.stats_Patient68012)
pK_Patient68012 = bcmvn_Patient68012 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68012 = as.numeric(as.character(pK_Patient68012[[1]]))
annotations_Patient68012 = Patient68012@meta.data$seurat_clusters
homotypic.prop_Patient68012 = modelHomotypic(annotations_Patient68012)
nExp_poi_Patient68012 = round(0.05*nrow(Patient68012@meta.data))
nExp_poi_adj_Patient68012 = round(nExp_poi_Patient68012*(1-homotypic.prop_Patient68012))
Patient68012 = doubletFinder(Patient68012, PCs = 1:10, pN = 0.25, pK = pK_Patient68012, nExp = nExp_poi_adj_Patient68012, reuse.pANN = FALSE, sct = FALSE)
Patient68012@meta.data$IsDoublet = Patient68012@meta.data[, grep("DF.classifications", names(Patient68012@meta.data), value = TRUE)] != "Singlet"
Patient68012@meta.data$ProbDoublet = Patient68012@meta.data[, grep("pANN", names(Patient68012@meta.data), value = TRUE)]
Patient68012@meta.data = Patient68012@meta.data[, !grepl("DF.classifications|pANN", names(Patient68012@meta.data))]


# 6801-EO-3 sample ####################################################################################################################
sc2=load10X("6801-EO/10x_analysis_6801-EO/Sample_6801-EO-3-GEX_TTATTCGAGG-AGCAGGACAG")
sc2 = autoEstCont(sc2)
out2 = adjustCounts(sc2, roundToInt = TRUE)
Patient68013 = CreateSeuratObject(out2)
Patient68013$EdgarID = "Patient68013"
Patient68013$ID = "Maze3"
Patient68013$Cohort = "FSGS"
Patient68013$Sex = "Male"
Patient68013$Age = "7"
Patient68013$APOL_Allele_Number = "0"
Patient68013$APOL_Alleles = "G0G0"
Patient68013$N264K = "No"
Patient68013$eGFR_Bx = "112.467"
Patient68013$UPCR_Bx = "6.64"
Patient68013$percent.mt = PercentageFeatureSet(Patient68013, pattern="^MT-")
Patient68013 = subset(Patient68013, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68013 = NormalizeData(Patient68013, normalization.method="LogNormalize",scale.factor = 10000)
Patient68013 = FindVariableFeatures(Patient68013, selection.method="vst",nfeatures= 2000)
Patient68013 = ScaleData(Patient68013)
Patient68013 = RunPCA(Patient68013, features = VariableFeatures(object = Patient68013))
Patient68013 = FindNeighbors(Patient68013, dims = 1:30)
Patient68013 = FindClusters(Patient68013, resolution = 0.6)
Patient68013 = RunUMAP(Patient68013, dims = 1:30)
sweep.res.list_Patient68013 = paramSweep(Patient68013, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68013 = summarizeSweep(sweep.res.list_Patient68013, GT = FALSE)
bcmvn_Patient68013 = find.pK(sweep.stats_Patient68013)
pK_Patient68013 = bcmvn_Patient68013 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68013 = as.numeric(as.character(pK_Patient68013[[1]]))
annotations_Patient68013 = Patient68013@meta.data$seurat_clusters
homotypic.prop_Patient68013 = modelHomotypic(annotations_Patient68013)
nExp_poi_Patient68013 = round(0.05*nrow(Patient68013@meta.data))
nExp_poi_adj_Patient68013 = round(nExp_poi_Patient68013*(1-homotypic.prop_Patient68013))
Patient68013 = doubletFinder(Patient68013, PCs = 1:10, pN = 0.25, pK = pK_Patient68013, nExp = nExp_poi_adj_Patient68013, reuse.pANN = FALSE, sct = FALSE)
Patient68013@meta.data$IsDoublet = Patient68013@meta.data[, grep("DF.classifications", names(Patient68013@meta.data), value = TRUE)] != "Singlet"
Patient68013@meta.data$ProbDoublet = Patient68013@meta.data[, grep("pANN", names(Patient68013@meta.data), value = TRUE)]
Patient68013@meta.data = Patient68013@meta.data[, !grepl("DF.classifications|pANN", names(Patient68013@meta.data))]


# 6801-EO-4 sample ####################################################################################################################
sc3=load10X("6801-EO/10x_analysis_6801-EO/Sample_6801-EO-4-GEX_AAGATTGGAT-AAATCCCGCT")
sc3 = autoEstCont(sc3)
out3 = adjustCounts(sc3, roundToInt = TRUE)
Patient68014 = CreateSeuratObject(out3)
Patient68014$EdgarID = "Patient68014"
Patient68014$ID = "Maze4"
Patient68014$Cohort = "MCD"
Patient68014$Sex = "Male"
Patient68014$Age = "15"
Patient68014$APOL_Allele_Number = "0"
Patient68014$APOL_Alleles = "G0G0"
Patient68014$N264K = "No"
Patient68014$eGFR_Bx = "149.864"
Patient68014$UPCR_Bx = "0.0388"
Patient68014$percent.mt = PercentageFeatureSet(Patient68014, pattern="^MT-")
Patient68014 = subset(Patient68014, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68014 = NormalizeData(Patient68014, normalization.method="LogNormalize",scale.factor = 10000)
Patient68014 = FindVariableFeatures(Patient68014, selection.method="vst",nfeatures= 2000)
Patient68014 = ScaleData(Patient68014)
Patient68014 = RunPCA(Patient68014, features = VariableFeatures(object = Patient68014))
Patient68014 = FindNeighbors(Patient68014, dims = 1:30)
Patient68014 = FindClusters(Patient68014, resolution = 0.6)
Patient68014 = RunUMAP(Patient68014, dims = 1:30)
sweep.res.list_Patient68014 = paramSweep(Patient68014, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68014 = summarizeSweep(sweep.res.list_Patient68014, GT = FALSE)
bcmvn_Patient68014 = find.pK(sweep.stats_Patient68014)
pK_Patient68014 = bcmvn_Patient68014 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68014 = as.numeric(as.character(pK_Patient68014[[1]]))
annotations_Patient68014 = Patient68014@meta.data$seurat_clusters
homotypic.prop_Patient68014 = modelHomotypic(annotations_Patient68014)
nExp_poi_Patient68014 = round(0.05*nrow(Patient68014@meta.data))
nExp_poi_adj_Patient68014 = round(nExp_poi_Patient68014*(1-homotypic.prop_Patient68014))
Patient68014 = doubletFinder(Patient68014, PCs = 1:10, pN = 0.25, pK = pK_Patient68014, nExp = nExp_poi_adj_Patient68014, reuse.pANN = FALSE, sct = FALSE)
Patient68014@meta.data$IsDoublet = Patient68014@meta.data[, grep("DF.classifications", names(Patient68014@meta.data), value = TRUE)] != "Singlet"
Patient68014@meta.data$ProbDoublet = Patient68014@meta.data[, grep("pANN", names(Patient68014@meta.data), value = TRUE)]
Patient68014@meta.data = Patient68014@meta.data[, !grepl("DF.classifications|pANN", names(Patient68014@meta.data))]


# 6817-EO-1 sample ####################################################################################################################
sc4=load10X("6817-EO/10x_analysis_6817-EO/Sample_6817-EO-1-GEX_TGTAGTCATT-TACGATCAAG")
sc4 = autoEstCont(sc4)
out4 = adjustCounts(sc4, roundToInt = TRUE)
Patient68171 = CreateSeuratObject(out4)
Patient68171$EdgarID = "Patient68171"
Patient68171$ID = "Maze5"
Patient68171$Cohort = "FSGS"
Patient68171$Sex = "Male"
Patient68171$Age = "16"
Patient68171$APOL_Allele_Number = "2"
Patient68171$APOL_Alleles = "G1G2"
Patient68171$N264K = "G2"
Patient68171$eGFR_Bx = "114.61"
Patient68171$UPCR_Bx = "0.5"
Patient68171$percent.mt = PercentageFeatureSet(Patient68171, pattern="^MT-")
Patient68171 = subset(Patient68171, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68171 = NormalizeData(Patient68171, normalization.method="LogNormalize",scale.factor = 10000)
Patient68171 = FindVariableFeatures(Patient68171, selection.method="vst",nfeatures= 2000)
Patient68171 = ScaleData(Patient68171)
Patient68171 = RunPCA(Patient68171, features = VariableFeatures(object = Patient68171))
Patient68171 = FindNeighbors(Patient68171, dims = 1:30)
Patient68171 = FindClusters(Patient68171, resolution = 0.6)
Patient68171 = RunUMAP(Patient68171, dims = 1:30)
sweep.res.list_Patient68171 = paramSweep(Patient68171, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68171 = summarizeSweep(sweep.res.list_Patient68171, GT = FALSE)
bcmvn_Patient68171 = find.pK(sweep.stats_Patient68171)
pK_Patient68171 = bcmvn_Patient68171 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68171 = as.numeric(as.character(pK_Patient68171[[1]]))
annotations_Patient68171 = Patient68171@meta.data$seurat_clusters
homotypic.prop_Patient68171 = modelHomotypic(annotations_Patient68171)
nExp_poi_Patient68171 = round(0.05*nrow(Patient68171@meta.data))
nExp_poi_adj_Patient68171 = round(nExp_poi_Patient68171*(1-homotypic.prop_Patient68171))
Patient68171 = doubletFinder(Patient68171, PCs = 1:10, pN = 0.25, pK = pK_Patient68171, nExp = nExp_poi_adj_Patient68171, reuse.pANN = FALSE, sct = FALSE)
Patient68171@meta.data$IsDoublet = Patient68171@meta.data[, grep("DF.classifications", names(Patient68171@meta.data), value = TRUE)] != "Singlet"
Patient68171@meta.data$ProbDoublet = Patient68171@meta.data[, grep("pANN", names(Patient68171@meta.data), value = TRUE)]
Patient68171@meta.data = Patient68171@meta.data[, !grepl("DF.classifications|pANN", names(Patient68171@meta.data))]


# 6817-EO-2 sample ####################################################################################################################
sc5=load10X("6817-EO/10x_analysis_6817-EO/Sample_6817-EO-2-GEX_ACAATGTGAA-TAACGGTACG")
sc5 = autoEstCont(sc5)
out5 = adjustCounts(sc5, roundToInt = TRUE)
Patient68172 = CreateSeuratObject(out5)
Patient68172$EdgarID = "Patient68172"
Patient68172$ID = "Maze6"
Patient68172$Cohort = "MCD"
Patient68172$Sex = "Female"
Patient68172$Age = "36"
Patient68172$APOL_Allele_Number = "1"
Patient68172$APOL_Alleles = "G0G2"
Patient68172$N264K = "G0"
Patient68172$eGFR_Bx = "90.385"
Patient68172$UPCR_Bx = "0.66"
Patient68172$percent.mt = PercentageFeatureSet(Patient68172, pattern="^MT-")
Patient68172 = subset(Patient68172, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68172 = NormalizeData(Patient68172, normalization.method="LogNormalize",scale.factor = 10000)
Patient68172 = FindVariableFeatures(Patient68172, selection.method="vst",nfeatures= 2000)
Patient68172 = ScaleData(Patient68172)
Patient68172 = RunPCA(Patient68172, features = VariableFeatures(object = Patient68172))
Patient68172 = FindNeighbors(Patient68172, dims = 1:30)
Patient68172 = FindClusters(Patient68172, resolution = 0.6)
Patient68172 = RunUMAP(Patient68172, dims = 1:30)
sweep.res.list_Patient68172 = paramSweep(Patient68172, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68172 = summarizeSweep(sweep.res.list_Patient68172, GT = FALSE)
bcmvn_Patient68172 = find.pK(sweep.stats_Patient68172)
pK_Patient68172 = bcmvn_Patient68172 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68172 = as.numeric(as.character(pK_Patient68172[[1]]))
annotations_Patient68172 = Patient68172@meta.data$seurat_clusters
homotypic.prop_Patient68172 = modelHomotypic(annotations_Patient68172)
nExp_poi_Patient68172 = round(0.05*nrow(Patient68172@meta.data))
nExp_poi_adj_Patient68172 = round(nExp_poi_Patient68172*(1-homotypic.prop_Patient68172))
Patient68172 = doubletFinder(Patient68172, PCs = 1:10, pN = 0.25, pK = pK_Patient68172, nExp = nExp_poi_adj_Patient68172, reuse.pANN = FALSE, sct = FALSE)
Patient68172@meta.data$IsDoublet = Patient68172@meta.data[, grep("DF.classifications", names(Patient68172@meta.data), value = TRUE)] != "Singlet"
Patient68172@meta.data$ProbDoublet = Patient68172@meta.data[, grep("pANN", names(Patient68172@meta.data), value = TRUE)]
Patient68172@meta.data = Patient68172@meta.data[, !grepl("DF.classifications|pANN", names(Patient68172@meta.data))]


# 6817-EO-3 sample ####################################################################################################################
sc6=load10X("6817-EO/10x_analysis_6817-EO/Sample_6817-EO-3-GEX_GTGGATCAAA-CAGGGTTGGC")
sc6 = autoEstCont(sc6)
out6 = adjustCounts(sc6, roundToInt = TRUE)
Patient68173 = CreateSeuratObject(out6)
Patient68173$EdgarID = "Patient68173"
Patient68173$ID = "Maze7"
Patient68173$Cohort = "FSGS"
Patient68173$Sex = "Male"
Patient68173$Age = "16"
Patient68173$APOL_Allele_Number = "2"
Patient68173$APOL_Alleles = "G1G2"
Patient68173$N264K = "No"
Patient68173$eGFR_Bx = "79.794"
Patient68173$UPCR_Bx = "0.7982"
Patient68173$percent.mt = PercentageFeatureSet(Patient68173, pattern="^MT-")
Patient68173 = subset(Patient68173, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68173 = NormalizeData(Patient68173, normalization.method="LogNormalize",scale.factor = 10000)
Patient68173 = FindVariableFeatures(Patient68173, selection.method="vst",nfeatures= 2000)
Patient68173 = ScaleData(Patient68173)
Patient68173 = RunPCA(Patient68173, features = VariableFeatures(object = Patient68173))
Patient68173 = FindNeighbors(Patient68173, dims = 1:30)
Patient68173 = FindClusters(Patient68173, resolution = 0.6)
Patient68173 = RunUMAP(Patient68173, dims = 1:30)
#The above is needed to run the data to a point where DoubletFinder can take over. the remaining steps for this sample identify doublets and then I remove them with the subset command below. I also clean up the resulting matrix with the last two lines in this sample. To save on space, I won't repeat this particular comment, but just kinda use this template down below.
sweep.res.list_Patient68173 = paramSweep(Patient68173, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68173 = summarizeSweep(sweep.res.list_Patient68173, GT = FALSE)
bcmvn_Patient68173 = find.pK(sweep.stats_Patient68173)
pK_Patient68173 = bcmvn_Patient68173 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68173 = as.numeric(as.character(pK_Patient68173[[1]]))
annotations_Patient68173 = Patient68173@meta.data$seurat_clusters
homotypic.prop_Patient68173 = modelHomotypic(annotations_Patient68173)
nExp_poi_Patient68173 = round(0.05*nrow(Patient68173@meta.data))
nExp_poi_adj_Patient68173 = round(nExp_poi_Patient68173*(1-homotypic.prop_Patient68173))
Patient68173 = doubletFinder(Patient68173, PCs = 1:10, pN = 0.25, pK = pK_Patient68173, nExp = nExp_poi_adj_Patient68173, reuse.pANN = FALSE, sct = FALSE)
#I ran table(Patient68173@meta.data$DF.classifications_0.25 (tabbed)) to determine the number of possible doublets that were identified.
Patient68173@meta.data$IsDoublet = Patient68173@meta.data[, grep("DF.classifications", names(Patient68173@meta.data), value = TRUE)] != "Singlet"
Patient68173@meta.data$ProbDoublet = Patient68173@meta.data[, grep("pANN", names(Patient68173@meta.data), value = TRUE)]
Patient68173@meta.data = Patient68173@meta.data[, !grepl("DF.classifications|pANN", names(Patient68173@meta.data))]


# 6817-EO-4 sample ####################################################################################################################
sc7=load10X("6817-EO/10x_analysis_6817-EO/Sample_6817-EO-4-GEX_TCTACCATTT-GACTCTCCCG")
sc7 = autoEstCont(sc7)
out7 = adjustCounts(sc7, roundToInt = TRUE)
Patient68174 = CreateSeuratObject(out7)
Patient68174$EdgarID = "Patient68174"
Patient68174$ID = "Maze8"
Patient68174$Cohort = "MCD"
Patient68174$Sex = "Female"
Patient68174$Age = "12"
Patient68174$APOL_Allele_Number = "1"
Patient68174$APOL_Alleles = "G0G2"
Patient68174$N264K = "No"
Patient68174$eGFR_Bx = "116.015"
Patient68174$UPCR_Bx = "0.0035"
Patient68174$percent.mt = PercentageFeatureSet(Patient68174, pattern="^MT-")
Patient68174 = subset(Patient68174, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68174 = NormalizeData(Patient68174, normalization.method="LogNormalize",scale.factor = 10000)
Patient68174 = FindVariableFeatures(Patient68174, selection.method="vst",nfeatures= 2000)
Patient68174 = ScaleData(Patient68174)
Patient68174 = RunPCA(Patient68174, features = VariableFeatures(object = Patient68174))
Patient68174 = FindNeighbors(Patient68174, dims = 1:30)
Patient68174 = FindClusters(Patient68174, resolution = 0.6)
Patient68174 = RunUMAP(Patient68174, dims = 1:30)
sweep.res.list_Patient68174 = paramSweep(Patient68174, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68174 = summarizeSweep(sweep.res.list_Patient68174, GT = FALSE)
bcmvn_Patient68174 = find.pK(sweep.stats_Patient68174)
pK_Patient68174 = bcmvn_Patient68174 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68174 = as.numeric(as.character(pK_Patient68174[[1]]))
annotations_Patient68174 = Patient68174@meta.data$seurat_clusters
homotypic.prop_Patient68174 = modelHomotypic(annotations_Patient68174)
nExp_poi_Patient68174 = round(0.05*nrow(Patient68174@meta.data))
nExp_poi_adj_Patient68174 = round(nExp_poi_Patient68174*(1-homotypic.prop_Patient68174))
Patient68174 = doubletFinder(Patient68174, PCs = 1:10, pN = 0.25, pK = pK_Patient68174, nExp = nExp_poi_adj_Patient68174, reuse.pANN = FALSE, sct = FALSE)
Patient68174@meta.data$IsDoublet = Patient68174@meta.data[, grep("DF.classifications", names(Patient68174@meta.data), value = TRUE)] != "Singlet"
Patient68174@meta.data$ProbDoublet = Patient68174@meta.data[, grep("pANN", names(Patient68174@meta.data), value = TRUE)]
Patient68174@meta.data = Patient68174@meta.data[, !grepl("DF.classifications|pANN", names(Patient68174@meta.data))]


# 6831-EO-1 sample  ####################################################################################################################
sc8=load10X("6831-EO/10x_analysis_6831-EO/Sample_6831-EO-1-GEX_CACGGTGAAT-TGTGACGAAC")
sc8 = autoEstCont(sc8)
out8 = adjustCounts(sc8, roundToInt = TRUE)
Patient68311 = CreateSeuratObject(out8)
Patient68311$EdgarID = "Patient68311"
Patient68311$ID = "Maze9"
Patient68311$Cohort = "FSGS"
Patient68311$Sex = "Male"
Patient68311$Age = "66"
Patient68311$APOL_Allele_Number = "0"
Patient68311$APOL_Alleles = "G0G0"
Patient68311$N264K = "No"
Patient68311$eGFR_Bx = "96.773"
Patient68311$UPCR_Bx = "0.2881"
Patient68311$percent.mt = PercentageFeatureSet(Patient68311, pattern="^MT-")
Patient68311 = subset(Patient68311, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68311 = NormalizeData(Patient68311, normalization.method="LogNormalize",scale.factor = 10000)
Patient68311 = FindVariableFeatures(Patient68311, selection.method="vst",nfeatures= 2000)
Patient68311 = ScaleData(Patient68311)
Patient68311 = RunPCA(Patient68311, features = VariableFeatures(object = Patient68311))
Patient68311 = FindNeighbors(Patient68311, dims = 1:30)
Patient68311 = FindClusters(Patient68311, resolution = 0.6)
Patient68311 = RunUMAP(Patient68311, dims = 1:30)
sweep.res.list_Patient68311 = paramSweep(Patient68311, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68311 = summarizeSweep(sweep.res.list_Patient68311, GT = FALSE)
bcmvn_Patient68311 = find.pK(sweep.stats_Patient68311)
pK_Patient68311 = bcmvn_Patient68311 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68311 = as.numeric(as.character(pK_Patient68311[[1]]))
annotations_Patient68311 = Patient68311@meta.data$seurat_clusters
homotypic.prop_Patient68311 = modelHomotypic(annotations_Patient68311)
nExp_poi_Patient68311 = round(0.05*nrow(Patient68311@meta.data))
nExp_poi_adj_Patient68311 = round(nExp_poi_Patient68311*(1-homotypic.prop_Patient68311))
Patient68311 = doubletFinder(Patient68311, PCs = 1:10, pN = 0.25, pK = pK_Patient68311, nExp = nExp_poi_adj_Patient68311, reuse.pANN = FALSE, sct = FALSE)
Patient68311@meta.data$IsDoublet = Patient68311@meta.data[, grep("DF.classifications", names(Patient68311@meta.data), value = TRUE)] != "Singlet"
Patient68311@meta.data$ProbDoublet = Patient68311@meta.data[, grep("pANN", names(Patient68311@meta.data), value = TRUE)]
Patient68311@meta.data = Patient68311@meta.data[, !grepl("DF.classifications|pANN", names(Patient68311@meta.data))]


# 6831-EO-2 sample ####################################################################################################################
sc9=load10X("6831-EO/10x_analysis_6831-EO/Sample_6831-EO-2-GEX_ATGGCTTGTG-CACAACATTC")
sc9 = autoEstCont(sc9)
out9 = adjustCounts(sc9, roundToInt = TRUE)
Patient68312 = CreateSeuratObject(out9)
Patient68312$EdgarID = "Patient68312"
Patient68312$ID = "Maze10"
Patient68312$Cohort = "MCD"
Patient68312$Sex = "Male"
Patient68312$Age = "50"
Patient68312$APOL_Allele_Number = "2"
Patient68312$APOL_Alleles = "G2G2"
Patient68312$N264K = "No"
Patient68312$eGFR_Bx = "23.285"
Patient68312$UPCR_Bx = "2.5272"
Patient68312$percent.mt = PercentageFeatureSet(Patient68312, pattern="^MT-")
Patient68312 = subset(Patient68312, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68312 = NormalizeData(Patient68312, normalization.method="LogNormalize",scale.factor = 10000)
Patient68312 = FindVariableFeatures(Patient68312, selection.method="vst",nfeatures= 2000)
Patient68312 = ScaleData(Patient68312)
Patient68312 = RunPCA(Patient68312, features = VariableFeatures(object = Patient68312))
Patient68312 = FindNeighbors(Patient68312, dims = 1:30)
Patient68312 = FindClusters(Patient68312, resolution = 0.6)
Patient68312 = RunUMAP(Patient68312, dims = 1:30)
sweep.res.list_Patient68312 = paramSweep(Patient68312, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68312 = summarizeSweep(sweep.res.list_Patient68312, GT = FALSE)
bcmvn_Patient68312 = find.pK(sweep.stats_Patient68312)
pK_Patient68312 = bcmvn_Patient68312 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68312 = as.numeric(as.character(pK_Patient68312[[1]]))
annotations_Patient68312 = Patient68312@meta.data$seurat_clusters
homotypic.prop_Patient68312 = modelHomotypic(annotations_Patient68312)
nExp_poi_Patient68312 = round(0.05*nrow(Patient68312@meta.data))
nExp_poi_adj_Patient68312 = round(nExp_poi_Patient68312*(1-homotypic.prop_Patient68312))
Patient68312 = doubletFinder(Patient68312, PCs = 1:10, pN = 0.25, pK = pK_Patient68312, nExp = nExp_poi_adj_Patient68312, reuse.pANN = FALSE, sct = FALSE)
Patient68312@meta.data$IsDoublet = Patient68312@meta.data[, grep("DF.classifications", names(Patient68312@meta.data), value = TRUE)] != "Singlet"
Patient68312@meta.data$ProbDoublet = Patient68312@meta.data[, grep("pANN", names(Patient68312@meta.data), value = TRUE)]
Patient68312@meta.data = Patient68312@meta.data[, !grepl("DF.classifications|pANN", names(Patient68312@meta.data))]


# 6831-EO-3 sample  ####################################################################################################################
sc10=load10X("6831-EO/10x_analysis_6831-EO/Sample_6831-EO-3-GEX_CCTTCTAGAG-TCGTTGTATT")
sc10 = autoEstCont(sc10)
out10 = adjustCounts(sc10, roundToInt = TRUE)
Patient68313 = CreateSeuratObject(out10)
Patient68313$EdgarID = "Patient68313"
Patient68313$ID = "Maze11"
Patient68313$Cohort = "MCD"
Patient68313$Sex = "Male"
Patient68313$Age = "4"
Patient68313$APOL_Allele_Number = "N/A"
Patient68313$APOL_Alleles = "N/A"
Patient68313$N264K = "N/A"
Patient68313$eGFR_Bx = "123.799"
Patient68313$UPCR_Bx = "3.4929"
Patient68313$percent.mt = PercentageFeatureSet(Patient68313, pattern="^MT-")
Patient68313 = subset(Patient68313, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68313 = NormalizeData(Patient68313, normalization.method="LogNormalize",scale.factor = 10000)
Patient68313 = FindVariableFeatures(Patient68313, selection.method="vst",nfeatures= 2000)
Patient68313 = ScaleData(Patient68313)
Patient68313 = RunPCA(Patient68313, features = VariableFeatures(object = Patient68313))
Patient68313 = FindNeighbors(Patient68313, dims = 1:30)
Patient68313 = FindClusters(Patient68313, resolution = 0.6)
Patient68313 = RunUMAP(Patient68313, dims = 1:30)
sweep.res.list_Patient68313 = paramSweep(Patient68313, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68313 = summarizeSweep(sweep.res.list_Patient68313, GT = FALSE)
bcmvn_Patient68313 = find.pK(sweep.stats_Patient68313)
pK_Patient68313 = bcmvn_Patient68313 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68313 = as.numeric(as.character(pK_Patient68313[[1]]))
annotations_Patient68313 = Patient68313@meta.data$seurat_clusters
homotypic.prop_Patient68313 = modelHomotypic(annotations_Patient68313)
nExp_poi_Patient68313 = round(0.05*nrow(Patient68313@meta.data))
nExp_poi_adj_Patient68313 = round(nExp_poi_Patient68313*(1-homotypic.prop_Patient68313))
Patient68313 = doubletFinder(Patient68313, PCs = 1:10, pN = 0.25, pK = pK_Patient68313, nExp = nExp_poi_adj_Patient68313, reuse.pANN = FALSE, sct = FALSE)
Patient68313@meta.data$IsDoublet = Patient68313@meta.data[, grep("DF.classifications", names(Patient68313@meta.data), value = TRUE)] != "Singlet"
Patient68313@meta.data$ProbDoublet = Patient68313@meta.data[, grep("pANN", names(Patient68313@meta.data), value = TRUE)]
Patient68313@meta.data = Patient68313@meta.data[, !grepl("DF.classifications|pANN", names(Patient68313@meta.data))]


# 6831-EO-4 sample ####################################################################################################################
sc11=load10X("6831-EO/10x_analysis_6831-EO/Sample_6831-EO-4-GEX_ACCAGACAAC-CCTAGTTCCT")
sc11 = autoEstCont(sc11)
out11 = adjustCounts(sc11, roundToInt = TRUE)
Patient68314 = CreateSeuratObject(out11)
Patient68314$EdgarID = "Patient68314"
Patient68314$ID = "Maze12"
Patient68314$Cohort = "FSGS"
Patient68314$Sex = "Female"
Patient68314$Age = "15"
Patient68314$APOL_Allele_Number = "0"
Patient68314$APOL_Alleles = "G0G0"
Patient68314$N264K = "No"
Patient68314$eGFR_Bx = "107.653"
Patient68314$UPCR_Bx = "11.1714"
Patient68314$percent.mt = PercentageFeatureSet(Patient68314, pattern="^MT-")
Patient68314 = subset(Patient68314, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68314 = NormalizeData(Patient68314, normalization.method="LogNormalize",scale.factor = 10000)
Patient68314 = FindVariableFeatures(Patient68314, selection.method="vst",nfeatures= 2000)
Patient68314 = ScaleData(Patient68314)
Patient68314 = RunPCA(Patient68314, features = VariableFeatures(object = Patient68314))
Patient68314 = FindNeighbors(Patient68314, dims = 1:30)
Patient68314 = FindClusters(Patient68314, resolution = 0.6)
Patient68314 = RunUMAP(Patient68314, dims = 1:30)
sweep.res.list_Patient68314 = paramSweep(Patient68314, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68314 = summarizeSweep(sweep.res.list_Patient68314, GT = FALSE)
bcmvn_Patient68314 = find.pK(sweep.stats_Patient68314)
pK_Patient68314 = bcmvn_Patient68314 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68314 = as.numeric(as.character(pK_Patient68314[[1]]))
annotations_Patient68314 = Patient68314@meta.data$seurat_clusters
homotypic.prop_Patient68314 = modelHomotypic(annotations_Patient68314)
nExp_poi_Patient68314 = round(0.05*nrow(Patient68314@meta.data))
nExp_poi_adj_Patient68314 = round(nExp_poi_Patient68314*(1-homotypic.prop_Patient68314))
Patient68314 = doubletFinder(Patient68314, PCs = 1:10, pN = 0.25, pK = pK_Patient68314, nExp = nExp_poi_adj_Patient68314, reuse.pANN = FALSE, sct = FALSE)
Patient68314@meta.data$IsDoublet = Patient68314@meta.data[, grep("DF.classifications", names(Patient68314@meta.data), value = TRUE)] != "Singlet"
Patient68314@meta.data$ProbDoublet = Patient68314@meta.data[, grep("pANN", names(Patient68314@meta.data), value = TRUE)]
Patient68314@meta.data = Patient68314@meta.data[, !grepl("DF.classifications|pANN", names(Patient68314@meta.data))]


# 6334-EO-1 sample ####################################################################################################################
sc12=load10X("6334-EO/10x_analysis_6334-EO/Sample_6334-EO-1-GEX_GCACTGAGAA-TTCACGCATA")
sc12 = autoEstCont(sc12)
out12 = adjustCounts(sc12, roundToInt = TRUE)
Patient63341 = CreateSeuratObject(out12)
Patient63341$EdgarID = "Patient63341"
Patient63341$ID = "Maze13"
Patient63341$Cohort = "N/A"
Patient63341$Sex = "Male"
Patient63341$Age = "60"
Patient63341$APOL_Allele_Number = "N/A"
Patient63341$APOL_Alleles = "N/A"
Patient63341$N264K = "N/A"
Patient63341$eGFR_Bx = "28.63"
Patient63341$UPCR_Bx = "5.1901"
Patient63341$percent.mt = PercentageFeatureSet(Patient63341, pattern="^MT-")
Patient63341 = subset(Patient63341, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient63341 = NormalizeData(Patient63341, normalization.method="LogNormalize",scale.factor = 10000)
Patient63341 = FindVariableFeatures(Patient63341, selection.method="vst",nfeatures= 2000)
Patient63341 = ScaleData(Patient63341)
Patient63341 = RunPCA(Patient63341, features = VariableFeatures(object = Patient63341))
Patient63341 = FindNeighbors(Patient63341, dims = 1:30)
Patient63341 = FindClusters(Patient63341, resolution = 0.6)
Patient63341 = RunUMAP(Patient63341, dims = 1:30)
sweep.res.list_Patient63341 = paramSweep(Patient63341, PCs = 1:10, sct=FALSE)
sweep.stats_Patient63341 = summarizeSweep(sweep.res.list_Patient63341, GT = FALSE)
bcmvn_Patient63341 = find.pK(sweep.stats_Patient63341)
pK_Patient63341 = bcmvn_Patient63341 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient63341 = as.numeric(as.character(pK_Patient63341[[1]]))
annotations_Patient63341 = Patient63341@meta.data$seurat_clusters
homotypic.prop_Patient63341 = modelHomotypic(annotations_Patient63341)
nExp_poi_Patient63341 = round(0.05*nrow(Patient63341@meta.data))
nExp_poi_adj_Patient63341 = round(nExp_poi_Patient63341*(1-homotypic.prop_Patient63341))
Patient63341 = doubletFinder(Patient63341, PCs = 1:10, pN = 0.25, pK = pK_Patient63341, nExp = nExp_poi_adj_Patient63341, reuse.pANN = FALSE, sct = FALSE)
Patient63341@meta.data$IsDoublet = Patient63341@meta.data[, grep("DF.classifications", names(Patient63341@meta.data), value = TRUE)] != "Singlet"
Patient63341@meta.data$ProbDoublet = Patient63341@meta.data[, grep("pANN", names(Patient63341@meta.data), value = TRUE)]
Patient63341@meta.data = Patient63341@meta.data[, !grepl("DF.classifications|pANN", names(Patient63341@meta.data))]


# 6334-EO-2 sample ####################################################################################################################
sc13=load10X("6334-EO/10x_analysis_6334-EO/Sample_6334-EO-2-GEX_GCTACAAAGC-AGGGCACGTG")
sc13 = autoEstCont(sc13)
out13 = adjustCounts(sc13, roundToInt = TRUE)
Patient63342 = CreateSeuratObject(out13)
Patient63342$EdgarID = "Patient63342"
Patient63342$ID = "Maze14"
Patient63342$Cohort = "N/A"
Patient63342$Sex = "Male"
Patient63342$Age = "N/A"
Patient63342$APOL_Allele_Number = "N/A"
Patient63342$APOL_Alleles = "N/A"
Patient63342$N264K = "N/A"
Patient63342$eGFR_Bx = "60.537"
Patient63342$UPCR_Bx = "3.641"
Patient63342$percent.mt = PercentageFeatureSet(Patient63342, pattern="^MT-")
Patient63342 = subset(Patient63342, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient63342= NormalizeData(Patient63342, normalization.method="LogNormalize",scale.factor = 10000)
Patient63342 = FindVariableFeatures(Patient63342, selection.method="vst",nfeatures= 2000)
Patient63342 = ScaleData(Patient63342)
Patient63342 = RunPCA(Patient63342, features = VariableFeatures(object = Patient63342))
Patient63342 = FindNeighbors(Patient63342, dims = 1:30)
Patient63342 = FindClusters(Patient63342, resolution = 0.6)
Patient63342 = RunUMAP(Patient63342, dims = 1:30)
sweep.res.list_Patient63342 = paramSweep(Patient63342, PCs = 1:10, sct=FALSE)
sweep.stats_Patient63342 = summarizeSweep(sweep.res.list_Patient63342, GT = FALSE)
bcmvn_Patient63342 = find.pK(sweep.stats_Patient63342)
pK_Patient63342 = bcmvn_Patient63342 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient63342 = as.numeric(as.character(pK_Patient63342[[1]]))
annotations_Patient63342 = Patient63342@meta.data$seurat_clusters
homotypic.prop_Patient63342 = modelHomotypic(annotations_Patient63342)
nExp_poi_Patient63342 = round(0.05*nrow(Patient63342@meta.data))
nExp_poi_adj_Patient63342 = round(nExp_poi_Patient63342*(1-homotypic.prop_Patient63342))
Patient63342 = doubletFinder(Patient63342, PCs = 1:10, pN = 0.25, pK = pK_Patient63342, nExp = nExp_poi_adj_Patient63342, reuse.pANN = FALSE, sct = FALSE)
Patient63342@meta.data$IsDoublet = Patient63342@meta.data[, grep("DF.classifications", names(Patient63342@meta.data), value = TRUE)] != "Singlet"
Patient63342@meta.data$ProbDoublet = Patient63342@meta.data[, grep("pANN", names(Patient63342@meta.data), value = TRUE)]
Patient63342@meta.data = Patient63342@meta.data[, !grepl("DF.classifications|pANN", names(Patient63342@meta.data))]


# 6348-EO-1 sample ####################################################################################################################
sc14=load10X("6348-EO/10x_analysis_6348-EO/Sample_6348-EO-1-GEX_ATAAGGATAC-CCCTATCTAT")
sc14 = autoEstCont(sc14)
out14 = adjustCounts(sc14, roundToInt = TRUE)
Patient63481 = CreateSeuratObject(out14)
Patient63481$EdgarID = "Patient63481"
Patient63481$ID = "Maze15"
Patient63481$Cohort = "AIN/DKD"
Patient63481$Sex = "Male"
Patient63481$Age = "N/A"
Patient63481$APOL_Allele_Number = "N/A"
Patient63481$APOL_Alleles = "N/A"
Patient63481$N264K = "N/A"
Patient63481$eGFR_Bx = "10.366"
Patient63481$UPCR_Bx = "12.0241"
Patient63481$percent.mt = PercentageFeatureSet(Patient63481, pattern="^MT-")
Patient63481= subset(Patient63481, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient63481= NormalizeData(Patient63481, normalization.method="LogNormalize",scale.factor = 10000)
Patient63481 = FindVariableFeatures(Patient63481, selection.method="vst",nfeatures= 2000)
Patient63481 = ScaleData(Patient63481)
Patient63481 = RunPCA(Patient63481, features = VariableFeatures(object = Patient63481))
Patient63481 = FindNeighbors(Patient63481, dims = 1:30)
Patient63481 = FindClusters(Patient63481, resolution = 0.6)
Patient63481 = RunUMAP(Patient63481, dims = 1:30)
sweep.res.list_Patient63481 = paramSweep(Patient63481, PCs = 1:10, sct=FALSE)
sweep.stats_Patient63481 = summarizeSweep(sweep.res.list_Patient63481, GT = FALSE)
bcmvn_Patient63481 = find.pK(sweep.stats_Patient63481)
pK_Patient63481 = bcmvn_Patient63481 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient63481 = as.numeric(as.character(pK_Patient63481[[1]]))
annotations_Patient63481 = Patient63481@meta.data$seurat_clusters
homotypic.prop_Patient63481 = modelHomotypic(annotations_Patient63481)
nExp_poi_Patient63481 = round(0.05*nrow(Patient63481@meta.data))
nExp_poi_adj_Patient63481 = round(nExp_poi_Patient63481*(1-homotypic.prop_Patient63481))
Patient63481 = doubletFinder(Patient63481, PCs = 1:10, pN = 0.25, pK = pK_Patient63481, nExp = nExp_poi_adj_Patient63481, reuse.pANN = FALSE, sct = FALSE)
Patient63481@meta.data$IsDoublet = Patient63481@meta.data[, grep("DF.classifications", names(Patient63481@meta.data), value = TRUE)] != "Singlet"
Patient63481@meta.data$ProbDoublet = Patient63481@meta.data[, grep("pANN", names(Patient63481@meta.data), value = TRUE)]
Patient63481@meta.data = Patient63481@meta.data[, !grepl("DF.classifications|pANN", names(Patient63481@meta.data))]


# 6348-EO-2 sample ####################################################################################################################
sc15=load10X("6348-EO/10x_analysis_6348-EO/Sample_6348-EO-2-GEX_AAGTGGAGAG-GTAACAGGAA")
sc15 = autoEstCont(sc15)
out15 = adjustCounts(sc15, roundToInt = TRUE)
Patient63482 = CreateSeuratObject(out15)
Patient63482$EdgarID = "Patient63482"
Patient63482$ID = "Maze16"
Patient63482$Cohort = "MGN"
Patient63482$Sex = "Female"
Patient63482$Age = "N/A"
Patient63482$APOL_Allele_Number = "N/A"
Patient63482$APOL_Alleles = "N/A"
Patient63482$N264K = "N/A"
Patient63482$eGFR_Bx = "77.826"
Patient63482$UPCR_Bx = "5.9386"
Patient63482$percent.mt = PercentageFeatureSet(Patient63482, pattern="^MT-")
Patient63482= subset(Patient63482, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient63482 = NormalizeData(Patient63482, normalization.method="LogNormalize",scale.factor = 10000)
Patient63482 = FindVariableFeatures(Patient63482, selection.method="vst",nfeatures= 2000)
Patient63482 = ScaleData(Patient63482)
Patient63482 = RunPCA(Patient63482, features = VariableFeatures(object = Patient63482))
Patient63482 = FindNeighbors(Patient63482, dims = 1:30)
Patient63482 = FindClusters(Patient63482, resolution = 0.6)
Patient63482 = RunUMAP(Patient63482, dims = 1:30)
sweep.res.list_Patient63482 = paramSweep(Patient63482, PCs = 1:10, sct=FALSE)
sweep.stats_Patient63482 = summarizeSweep(sweep.res.list_Patient63482, GT = FALSE)
bcmvn_Patient63482 = find.pK(sweep.stats_Patient63482)
pK_Patient63482 = bcmvn_Patient63482 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient63482 = as.numeric(as.character(pK_Patient63482[[1]]))
annotations_Patient63482 = Patient63482@meta.data$seurat_clusters
homotypic.prop_Patient63482 = modelHomotypic(annotations_Patient63482)
nExp_poi_Patient63482 = round(0.05*nrow(Patient63482@meta.data))
nExp_poi_adj_Patient63482 = round(nExp_poi_Patient63482*(1-homotypic.prop_Patient63482))
Patient63482 = doubletFinder(Patient63482, PCs = 1:10, pN = 0.25, pK = pK_Patient63482, nExp = nExp_poi_adj_Patient63482, reuse.pANN = FALSE, sct = FALSE)
Patient63482@meta.data$IsDoublet = Patient63482@meta.data[, grep("DF.classifications", names(Patient63482@meta.data), value = TRUE)] != "Singlet"
Patient63482@meta.data$ProbDoublet = Patient63482@meta.data[, grep("pANN", names(Patient63482@meta.data), value = TRUE)]
Patient63482@meta.data = Patient63482@meta.data[, !grepl("DF.classifications|pANN", names(Patient63482@meta.data))]


# 6793-EO-1 sample ####################################################################################################################
sc16=load10X("6793-EO/10x_analysis_6793-EO/Sample_6793-EO-1-GEX_CACCGCACCA-ATTGACAGTC")
sc16 = autoEstCont(sc16)
out16 = adjustCounts(sc16, roundToInt = TRUE)
Patient67931 = CreateSeuratObject(out16)
Patient67931$EdgarID = "Patient67931"
Patient67931$ID = "Maze17"
Patient67931$Cohort = "MCD"
Patient67931$Sex = "Female"
Patient67931$Age = "57"
Patient67931$APOL_Allele_Number = "0"
Patient67931$APOL_Alleles = "G0G0"
Patient67931$N264K = "No"
Patient67931$eGFR_Bx = "78.557"
Patient67931$UPCR_Bx = "2.5582"
Patient67931$percent.mt = PercentageFeatureSet(Patient67931, pattern="^MT-")
Patient67931= subset(Patient67931, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient67931 = NormalizeData(Patient67931, normalization.method="LogNormalize",scale.factor = 10000)
Patient67931 = FindVariableFeatures(Patient67931, selection.method="vst",nfeatures= 2000)
Patient67931 = ScaleData(Patient67931)
Patient67931 = RunPCA(Patient67931, features = VariableFeatures(object = Patient67931))
Patient67931 = FindNeighbors(Patient67931, dims = 1:30)
Patient67931 = FindClusters(Patient67931, resolution = 0.6)
Patient67931 = RunUMAP(Patient67931, dims = 1:30)
sweep.res.list_Patient67931 = paramSweep(Patient67931, PCs = 1:10, sct=FALSE)
sweep.stats_Patient67931 = summarizeSweep(sweep.res.list_Patient67931, GT = FALSE)
bcmvn_Patient67931 = find.pK(sweep.stats_Patient67931)
pK_Patient67931 = bcmvn_Patient67931 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient67931 = as.numeric(as.character(pK_Patient67931[[1]]))
annotations_Patient67931 = Patient67931@meta.data$seurat_clusters
homotypic.prop_Patient67931 = modelHomotypic(annotations_Patient67931)
nExp_poi_Patient67931 = round(0.05*nrow(Patient67931@meta.data))
nExp_poi_adj_Patient67931 = round(nExp_poi_Patient67931*(1-homotypic.prop_Patient67931))
Patient67931 = doubletFinder(Patient67931, PCs = 1:10, pN = 0.25, pK = pK_Patient67931, nExp = nExp_poi_adj_Patient67931, reuse.pANN = FALSE, sct = FALSE)
Patient67931@meta.data$IsDoublet = Patient67931@meta.data[, grep("DF.classifications", names(Patient67931@meta.data), value = TRUE)] != "Singlet"
Patient67931@meta.data$ProbDoublet = Patient67931@meta.data[, grep("pANN", names(Patient67931@meta.data), value = TRUE)]
Patient67931@meta.data = Patient67931@meta.data[, !grepl("DF.classifications|pANN", names(Patient67931@meta.data))]


# 6793-EO-2 sample ####################################################################################################################
sc17=load10X("6793-EO/10x_analysis_6793-EO/Sample_6793-EO-2-GEX_CGTCAAGGGC-GAGTGACCTA")
sc17 = autoEstCont(sc17)
out17 = adjustCounts(sc17, roundToInt = TRUE)
Patient67932 = CreateSeuratObject(out17)
Patient67932$EdgarID = "Patient67932"
Patient67932$ID = "Maze18"
Patient67932$Cohort = "FSGS"
Patient67932$Sex = "Male"
Patient67932$Age = "17"
Patient67932$APOL_Allele_Number = "2"
Patient67932$APOL_Alleles = "G1G1"
Patient67932$N264K = "No"
Patient67932$eGFR_Bx = "68.52"
Patient67932$UPCR_Bx = "0.2652"
Patient67932$percent.mt = PercentageFeatureSet(Patient67932, pattern="^MT-")
Patient67932= subset(Patient67932, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient67932 = NormalizeData(Patient67932, normalization.method="LogNormalize",scale.factor = 10000)
Patient67932 = FindVariableFeatures(Patient67932, selection.method="vst",nfeatures= 2000)
Patient67932 = ScaleData(Patient67932)
Patient67932 = RunPCA(Patient67932, features = VariableFeatures(object = Patient67932))
Patient67932 = FindNeighbors(Patient67932, dims = 1:30)
Patient67932 = FindClusters(Patient67932, resolution = 0.6)
Patient67932 = RunUMAP(Patient67932, dims = 1:30)
sweep.res.list_Patient67932 = paramSweep(Patient67932, PCs = 1:10, sct=FALSE)
sweep.stats_Patient67932 = summarizeSweep(sweep.res.list_Patient67932, GT = FALSE)
bcmvn_Patient67932 = find.pK(sweep.stats_Patient67932)
pK_Patient67932 = bcmvn_Patient67932 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient67932 = as.numeric(as.character(pK_Patient67932[[1]]))
annotations_Patient67932 = Patient67932@meta.data$seurat_clusters
homotypic.prop_Patient67932 = modelHomotypic(annotations_Patient67932)
nExp_poi_Patient67932 = round(0.05*nrow(Patient67932@meta.data))
nExp_poi_adj_Patient67932 = round(nExp_poi_Patient67932*(1-homotypic.prop_Patient67932))
Patient67932 = doubletFinder(Patient67932, PCs = 1:10, pN = 0.25, pK = pK_Patient67932, nExp = nExp_poi_adj_Patient67932, reuse.pANN = FALSE, sct = FALSE)
Patient67932@meta.data$IsDoublet = Patient67932@meta.data[, grep("DF.classifications", names(Patient67932@meta.data), value = TRUE)] != "Singlet"
Patient67932@meta.data$ProbDoublet = Patient67932@meta.data[, grep("pANN", names(Patient67932@meta.data), value = TRUE)]
Patient67932@meta.data = Patient67932@meta.data[, !grepl("DF.classifications|pANN", names(Patient67932@meta.data))]


# 6793-EO-3 sample ####################################################################################################################
sc18=load10X("6793-EO/10x_analysis_6793-EO/Sample_6793-EO-3-GEX_TCGTCAAGAT-CCTGAGTTGC")
sc18 = autoEstCont(sc18)
out18 = adjustCounts(sc18, roundToInt = TRUE)
Patient67933 = CreateSeuratObject(out18)
Patient67933$EdgarID = "Patient67933"
Patient67933$ID = "Maze19"
Patient67933$Cohort = "FSGS"
Patient67933$Sex = "Female"
Patient67933$Age = "4"
Patient67933$APOL_Allele_Number = "1"
Patient67933$APOL_Alleles = "G0G1"
Patient67933$N264K = "No"
Patient67933$eGFR_Bx = "64.046"
Patient67933$UPCR_Bx = "1.0247"
Patient67933$percent.mt = PercentageFeatureSet(Patient67933, pattern="^MT-")
Patient67933= subset(Patient67933, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient67933 = NormalizeData(Patient67933, normalization.method="LogNormalize",scale.factor = 10000)
Patient67933 = FindVariableFeatures(Patient67933, selection.method="vst",nfeatures= 2000)
Patient67933 = ScaleData(Patient67933)
Patient67933 = RunPCA(Patient67933, features = VariableFeatures(object = Patient67933))
Patient67933 = FindNeighbors(Patient67933, dims = 1:30)
Patient67933 = FindClusters(Patient67933, resolution = 0.6)
Patient67933 = RunUMAP(Patient67933, dims = 1:30)
sweep.res.list_Patient67933 = paramSweep(Patient67933, PCs = 1:10, sct=FALSE)
sweep.stats_Patient67933 = summarizeSweep(sweep.res.list_Patient67933, GT = FALSE)
bcmvn_Patient67933 = find.pK(sweep.stats_Patient67933)
pK_Patient67933 = bcmvn_Patient67933 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient67933 = as.numeric(as.character(pK_Patient67933[[1]]))
annotations_Patient67933 = Patient67933@meta.data$seurat_clusters
homotypic.prop_Patient67933 = modelHomotypic(annotations_Patient67933)
nExp_poi_Patient67933 = round(0.05*nrow(Patient67933@meta.data))
nExp_poi_adj_Patient67933 = round(nExp_poi_Patient67933*(1-homotypic.prop_Patient67933))
Patient67933 = doubletFinder(Patient67933, PCs = 1:10, pN = 0.25, pK = pK_Patient67933, nExp = nExp_poi_adj_Patient67933, reuse.pANN = FALSE, sct = FALSE)
Patient67933@meta.data$IsDoublet = Patient67933@meta.data[, grep("DF.classifications", names(Patient67933@meta.data), value = TRUE)] != "Singlet"
Patient67933@meta.data$ProbDoublet = Patient67933@meta.data[, grep("pANN", names(Patient67933@meta.data), value = TRUE)]
Patient67933@meta.data = Patient67933@meta.data[, !grepl("DF.classifications|pANN", names(Patient67933@meta.data))]


# 6839-EO-1 sample ####################################################################################################################
sc19=load10X("6839-EO/10x_analysis_6839-EO/Sample_6839-EO-1-GEX_TAAGCAAC-TTGAGTAT")
sc19 = autoEstCont(sc19)
out19 = adjustCounts(sc19, roundToInt = TRUE)
Patient68391 = CreateSeuratObject(out19)
Patient68391$EdgarID = "Patient68391"
Patient68391$ID = "Maze20"
Patient68391$Cohort = "MCD"
Patient68391$Sex = "Male"
Patient68391$Age = "24"
Patient68391$APOL_Allele_Number = "0"
Patient68391$APOL_Alleles = "G0G0"
Patient68391$N264K = "No"
Patient68391$eGFR_Bx = "61.602"
Patient68391$UPCR_Bx = "6.1707"
Patient68391$percent.mt = PercentageFeatureSet(Patient68391, pattern="^MT-")
Patient68391 = subset(Patient68391, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68391 = NormalizeData(Patient68391, normalization.method="LogNormalize",scale.factor = 10000)
Patient68391 = FindVariableFeatures(Patient68391, selection.method="vst",nfeatures= 2000)
Patient68391 = ScaleData(Patient68391)
Patient68391 = RunPCA(Patient68391, features = VariableFeatures(object = Patient68391))
Patient68391 = FindNeighbors(Patient68391, dims = 1:30)
Patient68391 = FindClusters(Patient68391, resolution = 0.6)
Patient68391 = RunUMAP(Patient68391, dims = 1:30)
sweep.res.list_Patient68391 = paramSweep(Patient68391, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68391 = summarizeSweep(sweep.res.list_Patient68391, GT = FALSE)
bcmvn_Patient68391 = find.pK(sweep.stats_Patient68391)
pK_Patient68391 = bcmvn_Patient68391 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68391 = as.numeric(as.character(pK_Patient68391[[1]]))
annotations_Patient68391 = Patient68391@meta.data$seurat_clusters
homotypic.prop_Patient68391 = modelHomotypic(annotations_Patient68391)
nExp_poi_Patient68391 = round(0.05*nrow(Patient68391@meta.data))
nExp_poi_adj_Patient68391 = round(nExp_poi_Patient68391*(1-homotypic.prop_Patient68391))
Patient68391 = doubletFinder(Patient68391, PCs = 1:10, pN = 0.25, pK = pK_Patient68391, nExp = nExp_poi_adj_Patient68391, reuse.pANN = FALSE, sct = FALSE)
Patient68391@meta.data$IsDoublet = Patient68391@meta.data[, grep("DF.classifications", names(Patient68391@meta.data), value = TRUE)] != "Singlet"
Patient68391@meta.data$ProbDoublet = Patient68391@meta.data[, grep("pANN", names(Patient68391@meta.data), value = TRUE)]
Patient68391@meta.data = Patient68391@meta.data[, !grepl("DF.classifications|pANN", names(Patient68391@meta.data))]


# 6839-EO-2 sample ####################################################################################################################
sc20=load10X("6839-EO/10x_analysis_6839-EO/Sample_6839-EO-2-GEX_ATAAGGAT-CCCTATCT")
sc20 = autoEstCont(sc20)
out20 = adjustCounts(sc20, roundToInt = TRUE)
Patient68392 = CreateSeuratObject(out20)
Patient68392$EdgarID = "Patient68392"
Patient68392$ID = "Maze21"
Patient68392$Cohort = "FSGS"
Patient68392$Sex = "Male"
Patient68392$Age = "77"
Patient68392$APOL_Allele_Number = "1"
Patient68392$APOL_Alleles = "G0G1"
Patient68392$N264K = "No"
Patient68392$eGFR_Bx = "33.666"
Patient68392$UPCR_Bx = "4.5162"
Patient68392$percent.mt = PercentageFeatureSet(Patient68392, pattern="^MT-")
Patient68392 = subset(Patient68392, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68392 = NormalizeData(Patient68392, normalization.method="LogNormalize",scale.factor = 10000)
Patient68392 = FindVariableFeatures(Patient68392, selection.method="vst",nfeatures= 2000)
Patient68392 = ScaleData(Patient68392)
Patient68392 = RunPCA(Patient68392, features = VariableFeatures(object = Patient68392))
Patient68392 = FindNeighbors(Patient68392, dims = 1:30)
Patient68392 = FindClusters(Patient68392, resolution = 0.6)
Patient68392 = RunUMAP(Patient68392, dims = 1:30)
sweep.res.list_Patient68392 = paramSweep(Patient68392, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68392 = summarizeSweep(sweep.res.list_Patient68392, GT = FALSE)
bcmvn_Patient68392 = find.pK(sweep.stats_Patient68392)
pK_Patient68392 = bcmvn_Patient68392 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68392 = as.numeric(as.character(pK_Patient68392[[1]]))
annotations_Patient68392 = Patient68392@meta.data$seurat_clusters
homotypic.prop_Patient68392 = modelHomotypic(annotations_Patient68392)
nExp_poi_Patient68392 = round(0.05*nrow(Patient68392@meta.data))
nExp_poi_adj_Patient68392 = round(nExp_poi_Patient68392*(1-homotypic.prop_Patient68392))
Patient68392 = doubletFinder(Patient68392, PCs = 1:10, pN = 0.25, pK = pK_Patient68392, nExp = nExp_poi_adj_Patient68392, reuse.pANN = FALSE, sct = FALSE)
Patient68392@meta.data$IsDoublet = Patient68392@meta.data[, grep("DF.classifications", names(Patient68392@meta.data), value = TRUE)] != "Singlet"
Patient68392@meta.data$ProbDoublet = Patient68392@meta.data[, grep("pANN", names(Patient68392@meta.data), value = TRUE)]
Patient68392@meta.data = Patient68392@meta.data[, !grepl("DF.classifications|pANN", names(Patient68392@meta.data))]


# 6839-EO-3 sample ####################################################################################################################
sc21=load10X("6839-EO/10x_analysis_6839-EO/Sample_6839-EO-3-GEX_AAGTGGAG-GTAACAGG")
sc21 = autoEstCont(sc21)
out21 = adjustCounts(sc21, roundToInt = TRUE)
Patient68393 = CreateSeuratObject(out21)
Patient68393$EdgarID = "Patient68393"
Patient68393$ID = "Maze22"
Patient68393$Cohort = "FSGS"
Patient68393$Sex = "Male"
Patient68393$Age = "20"
Patient68393$APOL_Allele_Number = "2"
Patient68393$APOL_Alleles = "G1G1"
Patient68393$N264K = "No"
Patient68393$eGFR_Bx = "19.367"
Patient68393$UPCR_Bx = "3.121"
Patient68393$percent.mt = PercentageFeatureSet(Patient68393, pattern="^MT-")
Patient68393 = subset(Patient68393, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68393 = NormalizeData(Patient68393, normalization.method="LogNormalize",scale.factor = 10000)
Patient68393 = FindVariableFeatures(Patient68393, selection.method="vst",nfeatures= 2000)
Patient68393 = ScaleData(Patient68393)
Patient68393 = RunPCA(Patient68393, features = VariableFeatures(object = Patient68393))
Patient68393 = FindNeighbors(Patient68393, dims = 1:30)
Patient68393 = FindClusters(Patient68393, resolution = 0.6)
Patient68393 = RunUMAP(Patient68393, dims = 1:30)
sweep.res.list_Patient68393 = paramSweep(Patient68393, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68393 = summarizeSweep(sweep.res.list_Patient68393, GT = FALSE)
bcmvn_Patient68393 = find.pK(sweep.stats_Patient68393)
pK_Patient68393 = bcmvn_Patient68393 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68393 = as.numeric(as.character(pK_Patient68393[[1]]))
annotations_Patient68393 = Patient68393@meta.data$seurat_clusters
homotypic.prop_Patient68393 = modelHomotypic(annotations_Patient68393)
nExp_poi_Patient68393 = round(0.05*nrow(Patient68393@meta.data))
nExp_poi_adj_Patient68393 = round(nExp_poi_Patient68393*(1-homotypic.prop_Patient68393))
Patient68393 = doubletFinder(Patient68393, PCs = 1:10, pN = 0.25, pK = pK_Patient68393, nExp = nExp_poi_adj_Patient68393, reuse.pANN = FALSE, sct = FALSE)
Patient68393@meta.data$IsDoublet = Patient68393@meta.data[, grep("DF.classifications", names(Patient68393@meta.data), value = TRUE)] != "Singlet"
Patient68393@meta.data$ProbDoublet = Patient68393@meta.data[, grep("pANN", names(Patient68393@meta.data), value = TRUE)]
Patient68393@meta.data = Patient68393@meta.data[, !grepl("DF.classifications|pANN", names(Patient68393@meta.data))]


# 6839-EO-4 sample  ####################################################################################################################
sc22=load10X("6839-EO/10x_analysis_6839-EO/Sample_6839-EO-4-GEX_TATTGAGG-CACTTACC")
sc22 = autoEstCont(sc22)
out22 = adjustCounts(sc22, roundToInt = TRUE)
Patient68394 = CreateSeuratObject(out22)
Patient68394$EdgarID = "Patient68394"
Patient68394$ID = "Maze23"
Patient68394$Cohort = "FSGS"
Patient68394$Sex = "Male"
Patient68394$Age = "51"
Patient68394$APOL_Allele_Number = "2"
Patient68394$APOL_Alleles = "G1G1"
Patient68394$N264K = "No"
Patient68394$eGFR_Bx = "21.972"
Patient68394$UPCR_Bx = "1.2121"
Patient68394$percent.mt = PercentageFeatureSet(Patient68394, pattern="^MT-")
Patient68394 = subset(Patient68394, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68394 = NormalizeData(Patient68394, normalization.method="LogNormalize",scale.factor = 10000)
Patient68394 = FindVariableFeatures(Patient68394, selection.method="vst",nfeatures= 2000)
Patient68394 = ScaleData(Patient68394)
Patient68394 = RunPCA(Patient68394, features = VariableFeatures(object = Patient68394))
Patient68394 = FindNeighbors(Patient68394, dims = 1:30)
Patient68394 = FindClusters(Patient68394, resolution = 0.6)
Patient68394 = RunUMAP(Patient68394, dims = 1:30)
sweep.res.list_Patient68394 = paramSweep(Patient68394, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68394 = summarizeSweep(sweep.res.list_Patient68394, GT = FALSE)
bcmvn_Patient68394 = find.pK(sweep.stats_Patient68394)
pK_Patient68394 = bcmvn_Patient68394 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68394 = as.numeric(as.character(pK_Patient68394[[1]]))
annotations_Patient68394 = Patient68394@meta.data$seurat_clusters
homotypic.prop_Patient68394 = modelHomotypic(annotations_Patient68394)
nExp_poi_Patient68394 = round(0.05*nrow(Patient68394@meta.data))
nExp_poi_adj_Patient68394 = round(nExp_poi_Patient68394*(1-homotypic.prop_Patient68394))
Patient68394 = doubletFinder(Patient68394, PCs = 1:10, pN = 0.25, pK = pK_Patient68394, nExp = nExp_poi_adj_Patient68394, reuse.pANN = FALSE, sct = FALSE)
Patient68394@meta.data$IsDoublet = Patient68394@meta.data[, grep("DF.classifications", names(Patient68394@meta.data), value = TRUE)] != "Singlet"
Patient68394@meta.data$ProbDoublet = Patient68394@meta.data[, grep("pANN", names(Patient68394@meta.data), value = TRUE)]
Patient68394@meta.data = Patient68394@meta.data[, !grepl("DF.classifications|pANN", names(Patient68394@meta.data))]


# 6871-EO-1 sample ####################################################################################################################
sc23=load10X("6871-EO/10x_analysis_6871-EO/Sample_6871-EO-1-GEX_CACCGCACCA-ATTGACAGTC")
sc23 = autoEstCont(sc23)
out23 = adjustCounts(sc23, roundToInt = TRUE)
Patient68711 = CreateSeuratObject(out23)
Patient68711$EdgarID = "Patient68711"
Patient68711$ID = "Maze24"
Patient68711$Cohort = "MCD"
Patient68711$Sex = "Female"
Patient68711$Age = "13"
Patient68711$APOL_Allele_Number = "1"
Patient68711$APOL_Alleles = "G0G1"
Patient68711$N264K = "No"
Patient68711$eGFR_Bx = "80.348"
Patient68711$UPCR_Bx = "11.7202"
Patient68711$percent.mt = PercentageFeatureSet(Patient68711, pattern="^MT-")
Patient68711 = subset(Patient68711, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68711 = NormalizeData(Patient68711, normalization.method="LogNormalize",scale.factor = 10000)
Patient68711 = FindVariableFeatures(Patient68711, selection.method="vst",nfeatures= 2000)
Patient68711 = ScaleData(Patient68711)
Patient68711 = RunPCA(Patient68711, features = VariableFeatures(object = Patient68711))
Patient68711 = FindNeighbors(Patient68711, dims = 1:30)
Patient68711 = FindClusters(Patient68711, resolution = 0.6)
Patient68711 = RunUMAP(Patient68711, dims = 1:30)
sweep.res.list_Patient68711 = paramSweep(Patient68711, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68711 = summarizeSweep(sweep.res.list_Patient68711, GT = FALSE)
bcmvn_Patient68711 = find.pK(sweep.stats_Patient68711)
pK_Patient68711 = bcmvn_Patient68711 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68711 = as.numeric(as.character(pK_Patient68711[[1]]))
annotations_Patient68711 = Patient68711@meta.data$seurat_clusters
homotypic.prop_Patient68711 = modelHomotypic(annotations_Patient68711)
nExp_poi_Patient68711 = round(0.05*nrow(Patient68711@meta.data))
nExp_poi_adj_Patient68711 = round(nExp_poi_Patient68711*(1-homotypic.prop_Patient68711))
Patient68711 = doubletFinder(Patient68711, PCs = 1:10, pN = 0.25, pK = pK_Patient68711, nExp = nExp_poi_adj_Patient68711, reuse.pANN = FALSE, sct = FALSE)
Patient68711@meta.data$IsDoublet = Patient68711@meta.data[, grep("DF.classifications", names(Patient68711@meta.data), value = TRUE)] != "Singlet"
Patient68711@meta.data$ProbDoublet = Patient68711@meta.data[, grep("pANN", names(Patient68711@meta.data), value = TRUE)]
Patient68711@meta.data = Patient68711@meta.data[, !grepl("DF.classifications|pANN", names(Patient68711@meta.data))]


# 6871-EO-2 sample ####################################################################################################################
sc24=load10X("6871-EO/10x_analysis_6871-EO/Sample_6871-EO-2-GEX_CGTCAAGGGC-GAGTGACCTA")
sc24 = autoEstCont(sc24)
out24 = adjustCounts(sc24, roundToInt = TRUE)
Patient68712 = CreateSeuratObject(out24)
Patient68712$EdgarID = "Patient68712"
Patient68712$ID = "Maze25"
Patient68712$Cohort = "FSGS"
Patient68712$Sex = "Male"
Patient68712$Age = "50"
Patient68712$APOL_Allele_Number = "N/A"
Patient68712$APOL_Alleles = "N/A"
Patient68712$N264K = "N/A"
Patient68712$eGFR_Bx = "36.354"
Patient68712$UPCR_Bx = "1.5229"
Patient68712$percent.mt = PercentageFeatureSet(Patient68712, pattern="^MT-")
Patient68712 = subset(Patient68712, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68712 = NormalizeData(Patient68712, normalization.method="LogNormalize",scale.factor = 10000)
Patient68712 = FindVariableFeatures(Patient68712, selection.method="vst",nfeatures= 2000)
Patient68712 = ScaleData(Patient68712)
Patient68712 = RunPCA(Patient68712, features = VariableFeatures(object = Patient68712))
Patient68712 = FindNeighbors(Patient68712, dims = 1:30)
Patient68712 = FindClusters(Patient68712, resolution = 0.6)
Patient68712 = RunUMAP(Patient68712, dims = 1:30)
sweep.res.list_Patient68712 = paramSweep(Patient68712, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68712 = summarizeSweep(sweep.res.list_Patient68712, GT = FALSE)
bcmvn_Patient68712 = find.pK(sweep.stats_Patient68712)
pK_Patient68712 = bcmvn_Patient68712 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68712 = as.numeric(as.character(pK_Patient68712[[1]]))
annotations_Patient68712 = Patient68712@meta.data$seurat_clusters
homotypic.prop_Patient68712 = modelHomotypic(annotations_Patient68712)
nExp_poi_Patient68712 = round(0.05*nrow(Patient68712@meta.data))
nExp_poi_adj_Patient68712 = round(nExp_poi_Patient68712*(1-homotypic.prop_Patient68712))
Patient68712 = doubletFinder(Patient68712, PCs = 1:10, pN = 0.25, pK = pK_Patient68712, nExp = nExp_poi_adj_Patient68712, reuse.pANN = FALSE, sct = FALSE)
Patient68712@meta.data$IsDoublet = Patient68712@meta.data[, grep("DF.classifications", names(Patient68712@meta.data), value = TRUE)] != "Singlet"
Patient68712@meta.data$ProbDoublet = Patient68712@meta.data[, grep("pANN", names(Patient68712@meta.data), value = TRUE)]
Patient68712@meta.data = Patient68712@meta.data[, !grepl("DF.classifications|pANN", names(Patient68712@meta.data))]


# 6871-EO-3 sample ####################################################################################################################
sc25=load10X("6871-EO/10x_analysis_6871-EO/Sample_6871-EO-3-GEX_TCGTCAAGAT-CCTGAGTTGC")
sc25 = autoEstCont(sc25)
out25 = adjustCounts(sc25, roundToInt = TRUE)
Patient68713 = CreateSeuratObject(out25)
Patient68713$EdgarID = "Patient68713"
Patient68713$ID = "Maze26"
Patient68713$Cohort = "MCD"
Patient68713$Sex = "Female"
Patient68713$Age = "5"
Patient68713$APOL_Allele_Number = "N/A"
Patient68713$APOL_Alleles = "N/A"
Patient68713$N264K = "N/A"
Patient68713$eGFR_Bx = "120.046"
Patient68713$UPCR_Bx = "11.6545"
Patient68713$percent.mt = PercentageFeatureSet(Patient68713, pattern="^MT-")
Patient68713 = subset(Patient68713, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68713 = NormalizeData(Patient68713, normalization.method="LogNormalize",scale.factor = 10000)
Patient68713 = FindVariableFeatures(Patient68713, selection.method="vst",nfeatures= 2000)
Patient68713 = ScaleData(Patient68713)
Patient68713 = RunPCA(Patient68713, features = VariableFeatures(object = Patient68713))
Patient68713 = FindNeighbors(Patient68713, dims = 1:30)
Patient68713 = FindClusters(Patient68713, resolution = 0.6)
Patient68713 = RunUMAP(Patient68713, dims = 1:30)
sweep.res.list_Patient68713 = paramSweep(Patient68713, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68713 = summarizeSweep(sweep.res.list_Patient68713, GT = FALSE)
bcmvn_Patient68713 = find.pK(sweep.stats_Patient68713)
pK_Patient68713 = bcmvn_Patient68713 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68713 = as.numeric(as.character(pK_Patient68713[[1]]))
annotations_Patient68713 = Patient68713@meta.data$seurat_clusters
homotypic.prop_Patient68713 = modelHomotypic(annotations_Patient68713)
nExp_poi_Patient68713 = round(0.05*nrow(Patient68713@meta.data))
nExp_poi_adj_Patient68713 = round(nExp_poi_Patient68713*(1-homotypic.prop_Patient68713))
Patient68713 = doubletFinder(Patient68713, PCs = 1:10, pN = 0.25, pK = pK_Patient68713, nExp = nExp_poi_adj_Patient68713, reuse.pANN = FALSE, sct = FALSE)
Patient68713@meta.data$IsDoublet = Patient68713@meta.data[, grep("DF.classifications", names(Patient68713@meta.data), value = TRUE)] != "Singlet"
Patient68713@meta.data$ProbDoublet = Patient68713@meta.data[, grep("pANN", names(Patient68713@meta.data), value = TRUE)]
Patient68713@meta.data = Patient68713@meta.data[, !grepl("DF.classifications|pANN", names(Patient68713@meta.data))]


# 6886-EO-1 sample ####################################################################################################################
sc26=load10X("6886-EO/10x_analysis_6886-EO/Sample_6886-EO-1-GEX_TCCGGGAC-TGGCATTC")
sc26 = autoEstCont(sc26)
out26 = adjustCounts(sc26, roundToInt = TRUE)
Patient68861 = CreateSeuratObject(out26)
Patient68861$EdgarID = "Patient68861"
Patient68861$ID = "Maze27"
Patient68861$Cohort = "MCD"
Patient68861$Sex = "Female"
Patient68861$Age = "57"
Patient68861$APOL_Allele_Number = "N/A"
Patient68861$APOL_Alleles = "N/A"
Patient68861$N264K = "N/A"
Patient68861$eGFR_Bx = "48.071"
Patient68861$UPCR_Bx = "0.1245"
Patient68861$percent.mt = PercentageFeatureSet(Patient68861, pattern="^MT-")
Patient68861 = subset(Patient68861, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68861 = NormalizeData(Patient68861, normalization.method="LogNormalize",scale.factor = 10000)
Patient68861 = FindVariableFeatures(Patient68861, selection.method="vst",nfeatures= 2000)
Patient68861 = ScaleData(Patient68861)
Patient68861 = RunPCA(Patient68861, features = VariableFeatures(object = Patient68861))
Patient68861 = FindNeighbors(Patient68861, dims = 1:30)
Patient68861 = FindClusters(Patient68861, resolution = 0.6)
Patient68861 = RunUMAP(Patient68861, dims = 1:30)
sweep.res.list_Patient68861 = paramSweep(Patient68861, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68861 = summarizeSweep(sweep.res.list_Patient68861, GT = FALSE)
bcmvn_Patient68861 = find.pK(sweep.stats_Patient68861)
pK_Patient68861 = bcmvn_Patient68861 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68861 = as.numeric(as.character(pK_Patient68861[[1]]))
annotations_Patient68861 = Patient68861@meta.data$seurat_clusters
homotypic.prop_Patient68861 = modelHomotypic(annotations_Patient68861)
nExp_poi_Patient68861 = round(0.05*nrow(Patient68861@meta.data))
nExp_poi_adj_Patient68861 = round(nExp_poi_Patient68861*(1-homotypic.prop_Patient68861))
Patient68861 = doubletFinder(Patient68861, PCs = 1:10, pN = 0.25, pK = pK_Patient68861, nExp = nExp_poi_adj_Patient68861, reuse.pANN = FALSE, sct = FALSE)
Patient68861@meta.data$IsDoublet = Patient68861@meta.data[, grep("DF.classifications", names(Patient68861@meta.data), value = TRUE)] != "Singlet"
Patient68861@meta.data$ProbDoublet = Patient68861@meta.data[, grep("pANN", names(Patient68861@meta.data), value = TRUE)]
Patient68861@meta.data = Patient68861@meta.data[, !grepl("DF.classifications|pANN", names(Patient68861@meta.data))]


# 6886-EO-2 sample ####################################################################################################################
sc27=load10X("6886-EO/10x_analysis_6886-EO/Sample_6886-EO-2-GEX_TTCACACC-GTGTACAC")
sc27 = autoEstCont(sc27)
out27 = adjustCounts(sc27, roundToInt = TRUE)
Patient68862 = CreateSeuratObject(out27)
Patient68862$EdgarID = "Patient68862"
Patient68862$ID = "Maze28"
Patient68862$Cohort = "FSGS"
Patient68862$Sex = "Male"
Patient68862$Age = "59"
Patient68862$APOL_Allele_Number = "N/A"
Patient68862$APOL_Alleles = "N/A"
Patient68862$N264K = "N/A"
Patient68862$eGFR_Bx = "37.513"
Patient68862$UPCR_Bx = "14.0893"
Patient68862$percent.mt = PercentageFeatureSet(Patient68862, pattern="^MT-")
Patient68862 = subset(Patient68862, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68862 = NormalizeData(Patient68862, normalization.method="LogNormalize",scale.factor = 10000)
Patient68862 = FindVariableFeatures(Patient68862, selection.method="vst",nfeatures= 2000)
Patient68862 = ScaleData(Patient68862)
Patient68862 = RunPCA(Patient68862, features = VariableFeatures(object = Patient68862))
Patient68862 = FindNeighbors(Patient68862, dims = 1:30)
Patient68862 = FindClusters(Patient68862, resolution = 0.6)
Patient68862 = RunUMAP(Patient68862, dims = 1:30)
sweep.res.list_Patient68862 = paramSweep(Patient68862, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68862 = summarizeSweep(sweep.res.list_Patient68862, GT = FALSE)
bcmvn_Patient68862 = find.pK(sweep.stats_Patient68862)
pK_Patient68862 = bcmvn_Patient68862 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68862 = as.numeric(as.character(pK_Patient68862[[1]]))
annotations_Patient68862 = Patient68862@meta.data$seurat_clusters
homotypic.prop_Patient68862 = modelHomotypic(annotations_Patient68862)
nExp_poi_Patient68862 = round(0.05*nrow(Patient68862@meta.data))
nExp_poi_adj_Patient68862 = round(nExp_poi_Patient68862*(1-homotypic.prop_Patient68862))
Patient68862 = doubletFinder(Patient68862, PCs = 1:10, pN = 0.25, pK = pK_Patient68862, nExp = nExp_poi_adj_Patient68862, reuse.pANN = FALSE, sct = FALSE)
Patient68862@meta.data$IsDoublet = Patient68862@meta.data[, grep("DF.classifications", names(Patient68862@meta.data), value = TRUE)] != "Singlet"
Patient68862@meta.data$ProbDoublet = Patient68862@meta.data[, grep("pANN", names(Patient68862@meta.data), value = TRUE)]
Patient68862@meta.data = Patient68862@meta.data[, !grepl("DF.classifications|pANN", names(Patient68862@meta.data))]


# 6886-EO-3 sample ####################################################################################################################
sc28=load10X("6886-EO/10x_analysis_6886-EO/Sample_6886-EO-3-GEX_GATAACCT-GTTTCTAA")
sc28 = autoEstCont(sc28)
out28 = adjustCounts(sc28, roundToInt = TRUE)
Patient68863 = CreateSeuratObject(out28)
Patient68863$EdgarID = "Patient68863"
Patient68863$ID = "Maze29"
Patient68863$Cohort = "MCD"
Patient68863$Sex = "Male"
Patient68863$Age = "58"
Patient68863$APOL_Allele_Number = "N/A"
Patient68863$APOL_Alleles = "N/A"
Patient68863$N264K = "N/A"
Patient68863$eGFR_Bx = "33.003"
Patient68863$UPCR_Bx = "20.61"
Patient68863$percent.mt = PercentageFeatureSet(Patient68863, pattern="^MT-")
Patient68863 = subset(Patient68863, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68863 = NormalizeData(Patient68863, normalization.method="LogNormalize",scale.factor = 10000)
Patient68863 = FindVariableFeatures(Patient68863, selection.method="vst",nfeatures= 2000)
Patient68863 = ScaleData(Patient68863)
Patient68863 = RunPCA(Patient68863, features = VariableFeatures(object = Patient68863))
Patient68863 = FindNeighbors(Patient68863, dims = 1:30)
Patient68863 = FindClusters(Patient68863, resolution = 0.6)
Patient68863 = RunUMAP(Patient68863, dims = 1:30)
sweep.res.list_Patient68863 = paramSweep(Patient68863, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68863 = summarizeSweep(sweep.res.list_Patient68863, GT = FALSE)
bcmvn_Patient68863 = find.pK(sweep.stats_Patient68863)
pK_Patient68863 = bcmvn_Patient68863 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68863 = as.numeric(as.character(pK_Patient68863[[1]]))
annotations_Patient68863 = Patient68863@meta.data$seurat_clusters
homotypic.prop_Patient68863 = modelHomotypic(annotations_Patient68863)
nExp_poi_Patient68863 = round(0.05*nrow(Patient68863@meta.data))
nExp_poi_adj_Patient68863 = round(nExp_poi_Patient68863*(1-homotypic.prop_Patient68863))
Patient68863 = doubletFinder(Patient68863, PCs = 1:10, pN = 0.25, pK = pK_Patient68863, nExp = nExp_poi_adj_Patient68863, reuse.pANN = FALSE, sct = FALSE)
Patient68863@meta.data$IsDoublet = Patient68863@meta.data[, grep("DF.classifications", names(Patient68863@meta.data), value = TRUE)] != "Singlet"
Patient68863@meta.data$ProbDoublet = Patient68863@meta.data[, grep("pANN", names(Patient68863@meta.data), value = TRUE)]
Patient68863@meta.data = Patient68863@meta.data[, !grepl("DF.classifications|pANN", names(Patient68863@meta.data))]


# 6886-EO-4 sample ####################################################################################################################
sc29=load10X("6886-EO/10x_analysis_6886-EO/Sample_6886-EO-4-GEX_ACAATCGA-CATTCCGT")
sc29 = autoEstCont(sc29)
out29 = adjustCounts(sc29, roundToInt = TRUE)
Patient68864 = CreateSeuratObject(out29)
Patient68864$EdgarID = "Patient68864"
Patient68864$ID = "Maze30"
Patient68864$Cohort = "MCD"
Patient68864$Sex = "Male"
Patient68864$Age = "6"
Patient68864$APOL_Allele_Number = "N/A"
Patient68864$APOL_Alleles = "N/A"
Patient68864$N264K = "N/A"
Patient68864$eGFR_Bx = "138.018"
Patient68864$UPCR_Bx = "29.1382"
Patient68864$percent.mt = PercentageFeatureSet(Patient68864, pattern="^MT-")
Patient68864 = subset(Patient68864, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68864 = NormalizeData(Patient68864, normalization.method="LogNormalize",scale.factor = 10000)
Patient68864 = FindVariableFeatures(Patient68864, selection.method="vst",nfeatures= 2000)
Patient68864 = ScaleData(Patient68864)
Patient68864 = RunPCA(Patient68864, features = VariableFeatures(object = Patient68864))
Patient68864 = FindNeighbors(Patient68864, dims = 1:30)
Patient68864 = FindClusters(Patient68864, resolution = 0.6)
Patient68864 = RunUMAP(Patient68864, dims = 1:30)
sweep.res.list_Patient68864 = paramSweep(Patient68864, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68864 = summarizeSweep(sweep.res.list_Patient68864, GT = FALSE)
bcmvn_Patient68864 = find.pK(sweep.stats_Patient68864)
pK_Patient68864 = bcmvn_Patient68864 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68864 = as.numeric(as.character(pK_Patient68864[[1]]))
annotations_Patient68864 = Patient68864@meta.data$seurat_clusters
homotypic.prop_Patient68864 = modelHomotypic(annotations_Patient68864)
nExp_poi_Patient68864 = round(0.05*nrow(Patient68864@meta.data))
nExp_poi_adj_Patient68864 = round(nExp_poi_Patient68864*(1-homotypic.prop_Patient68864))
Patient68864 = doubletFinder(Patient68864, PCs = 1:10, pN = 0.25, pK = pK_Patient68864, nExp = nExp_poi_adj_Patient68864, reuse.pANN = FALSE, sct = FALSE)
Patient68864@meta.data$IsDoublet = Patient68864@meta.data[, grep("DF.classifications", names(Patient68864@meta.data), value = TRUE)] != "Singlet"
Patient68864@meta.data$ProbDoublet = Patient68864@meta.data[, grep("pANN", names(Patient68864@meta.data), value = TRUE)]
Patient68864@meta.data = Patient68864@meta.data[, !grepl("DF.classifications|pANN", names(Patient68864@meta.data))]


# 6936-EO-1 sample ####################################################################################################################
sc30=load10X("6936-EO/10x_analysis_6936-EO/Sample_6936-EO-1-GEX_GTGGATCAAA-CAGGGTTGGC")
sc30 = autoEstCont(sc30)
out30 = adjustCounts(sc30, roundToInt = TRUE)
Patient69361 = CreateSeuratObject(out30)
Patient69361$EdgarID = "Patient69361"
Patient69361$ID = "Maze31"
Patient69361$Cohort = "FSGS"
Patient69361$Sex = "Male"
Patient69361$Age = "16"
Patient69361$APOL_Allele_Number = "2"
Patient69361$APOL_Alleles = "G1G1"
Patient69361$N264K = "No"
Patient69361$eGFR_Bx = "57.139"
Patient69361$UPCR_Bx = "2.3845"
Patient69361$percent.mt = PercentageFeatureSet(Patient69361, pattern="^MT-")
Patient69361 = subset(Patient69361, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient69361 = NormalizeData(Patient69361, normalization.method="LogNormalize",scale.factor = 10000)
Patient69361 = FindVariableFeatures(Patient69361, selection.method="vst",nfeatures= 2000)
Patient69361 = ScaleData(Patient69361)
Patient69361 = RunPCA(Patient69361, features = VariableFeatures(object = Patient69361))
Patient69361 = FindNeighbors(Patient69361, dims = 1:30)
Patient69361 = FindClusters(Patient69361, resolution = 0.6)
Patient69361 = RunUMAP(Patient69361, dims = 1:30)
sweep.res.list_Patient69361 = paramSweep(Patient69361, PCs = 1:10, sct=FALSE)
sweep.stats_Patient69361 = summarizeSweep(sweep.res.list_Patient69361, GT = FALSE)
bcmvn_Patient69361 = find.pK(sweep.stats_Patient69361)
pK_Patient69361 = bcmvn_Patient69361 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient69361 = as.numeric(as.character(pK_Patient69361[[1]]))
annotations_Patient69361 = Patient69361@meta.data$seurat_clusters
homotypic.prop_Patient69361 = modelHomotypic(annotations_Patient69361)
nExp_poi_Patient69361 = round(0.05*nrow(Patient69361@meta.data))
nExp_poi_adj_Patient69361 = round(nExp_poi_Patient69361*(1-homotypic.prop_Patient69361))
Patient69361 = doubletFinder(Patient69361, PCs = 1:10, pN = 0.25, pK = pK_Patient69361, nExp = nExp_poi_adj_Patient69361, reuse.pANN = FALSE, sct = FALSE)
Patient69361@meta.data$IsDoublet = Patient69361@meta.data[, grep("DF.classifications", names(Patient69361@meta.data), value = TRUE)] != "Singlet"
Patient69361@meta.data$ProbDoublet = Patient69361@meta.data[, grep("pANN", names(Patient69361@meta.data), value = TRUE)]
Patient69361@meta.data = Patient69361@meta.data[, !grepl("DF.classifications|pANN", names(Patient69361@meta.data))]


# 6936-EO-2 sample ####################################################################################################################
sc31=load10X("6936-EO/10x_analysis_6936-EO/Sample_6936-EO-2-GEX_TCTACCATTT-GACTCTCCCG")
sc31 = autoEstCont(sc31)
out31 = adjustCounts(sc31, roundToInt = TRUE)
Patient69362 = CreateSeuratObject(out31)
Patient69362$EdgarID = "Patient69362"
Patient69362$ID = "Maze32"
Patient69362$Cohort = "FSGS"
Patient69362$Sex = "Male"
Patient69362$Age = "53"
Patient69362$APOL_Allele_Number = "N/A"
Patient69362$APOL_Alleles = "N/A"
Patient69362$N264K = "N/A"
Patient69362$eGFR_Bx = "14.9625"
Patient69362$UPCR_Bx = "4.46039"
Patient69362$percent.mt = PercentageFeatureSet(Patient69362, pattern="^MT-")
Patient69362 = subset(Patient69362, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient69362 = NormalizeData(Patient69362, normalization.method="LogNormalize",scale.factor = 10000)
Patient69362 = FindVariableFeatures(Patient69362, selection.method="vst",nfeatures= 2000)
Patient69362 = ScaleData(Patient69362)
Patient69362 = RunPCA(Patient69362, features = VariableFeatures(object = Patient69362))
Patient69362 = FindNeighbors(Patient69362, dims = 1:30)
Patient69362 = FindClusters(Patient69362, resolution = 0.6)
Patient69362 = RunUMAP(Patient69362, dims = 1:30)
sweep.res.list_Patient69362 = paramSweep(Patient69362, PCs = 1:10, sct=FALSE)
sweep.stats_Patient69362 = summarizeSweep(sweep.res.list_Patient69362, GT = FALSE)
bcmvn_Patient69362 = find.pK(sweep.stats_Patient69362)
pK_Patient69362 = bcmvn_Patient69362 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient69362 = as.numeric(as.character(pK_Patient69362[[1]]))
annotations_Patient69362 = Patient69362@meta.data$seurat_clusters
homotypic.prop_Patient69362 = modelHomotypic(annotations_Patient69362)
nExp_poi_Patient69362 = round(0.05*nrow(Patient69362@meta.data))
nExp_poi_adj_Patient69362 = round(nExp_poi_Patient69362*(1-homotypic.prop_Patient69362))
Patient69362 = doubletFinder(Patient69362, PCs = 1:10, pN = 0.25, pK = pK_Patient69362, nExp = nExp_poi_adj_Patient69362, reuse.pANN = FALSE, sct = FALSE)
Patient69362@meta.data$IsDoublet = Patient69362@meta.data[, grep("DF.classifications", names(Patient69362@meta.data), value = TRUE)] != "Singlet"
Patient69362@meta.data$ProbDoublet = Patient69362@meta.data[, grep("pANN", names(Patient69362@meta.data), value = TRUE)]
Patient69362@meta.data = Patient69362@meta.data[, !grepl("DF.classifications|pANN", names(Patient69362@meta.data))]


# 6936-EO-3 sample ####################################################################################################################
sc32=load10X("6936-EO/10x_analysis_6936-EO/Sample_6936-EO-3-GEX_CAATCCCGAC-TACTACTCGG")
sc32 = autoEstCont(sc32)
out32 = adjustCounts(sc32, roundToInt = TRUE)
Patient69363 = CreateSeuratObject(out32)
Patient69363$EdgarID = "Patient69363"
Patient69363$ID = "Maze33"
Patient69363$Cohort = "MCD"
Patient69363$Sex = "Male"
Patient69363$Age = "10"
Patient69363$APOL_Allele_Number = "N/A"
Patient69363$APOL_Alleles = "N/A"
Patient69363$N264K = "N/A"
Patient69363$eGFR_Bx = "157.715"
Patient69363$UPCR_Bx = "3.789"
Patient69363$percent.mt = PercentageFeatureSet(Patient69363, pattern="^MT-")
Patient69363 = subset(Patient69363, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient69363 = NormalizeData(Patient69363, normalization.method="LogNormalize",scale.factor = 10000)
Patient69363 = FindVariableFeatures(Patient69363, selection.method="vst",nfeatures= 2000)
Patient69363 = ScaleData(Patient69363)
Patient69363 = RunPCA(Patient69363, features = VariableFeatures(object = Patient69363))
Patient69363 = FindNeighbors(Patient69363, dims = 1:30)
Patient69363 = FindClusters(Patient69363, resolution = 0.6)
Patient69363 = RunUMAP(Patient69363, dims = 1:30)
sweep.res.list_Patient69363 = paramSweep(Patient69363, PCs = 1:10, sct=FALSE)
sweep.stats_Patient69363 = summarizeSweep(sweep.res.list_Patient69363, GT = FALSE)
bcmvn_Patient69363 = find.pK(sweep.stats_Patient69363)
pK_Patient69363 = bcmvn_Patient69363 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient69363 = as.numeric(as.character(pK_Patient69363[[1]]))
annotations_Patient69363 = Patient69363@meta.data$seurat_clusters
homotypic.prop_Patient69363 = modelHomotypic(annotations_Patient69363)
nExp_poi_Patient69363 = round(0.05*nrow(Patient69363@meta.data))
nExp_poi_adj_Patient69363 = round(nExp_poi_Patient69363*(1-homotypic.prop_Patient69363))
Patient69363 = doubletFinder(Patient69363, PCs = 1:10, pN = 0.25, pK = pK_Patient69363, nExp = nExp_poi_adj_Patient69363, reuse.pANN = FALSE, sct = FALSE)
Patient69363@meta.data$IsDoublet = Patient69363@meta.data[, grep("DF.classifications", names(Patient69363@meta.data), value = TRUE)] != "Singlet"
Patient69363@meta.data$ProbDoublet = Patient69363@meta.data[, grep("pANN", names(Patient69363@meta.data), value = TRUE)]
Patient69363@meta.data = Patient69363@meta.data[, !grepl("DF.classifications|pANN", names(Patient69363@meta.data))]


# 6936-EO-4 sample ####################################################################################################################
sc33=load10X("6936-EO/10x_analysis_6936-EO/Sample_6936-EO-4-GEX_TTAATACGCG-ACCCGAGGTG")
sc33 = autoEstCont(sc33)
out33 = adjustCounts(sc33, roundToInt = TRUE)
Patient69364 = CreateSeuratObject(out33)
Patient69364$EdgarID = "Patient69364"
Patient69364$ID = "Maze34"
Patient69364$Cohort = "DKD"
Patient69364$Sex = "Male"
Patient69364$Age = "49"
Patient69364$APOL_Allele_Number = "N/A"
Patient69364$APOL_Alleles = "N/A"
Patient69364$N264K = "N/A"
Patient69364$eGFR_Bx = "50.356"
Patient69364$UPCR_Bx = "10.5455"
Patient69364$percent.mt = PercentageFeatureSet(Patient69364, pattern="^MT-")
Patient69364 = subset(Patient69364, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient69364 = NormalizeData(Patient69364, normalization.method="LogNormalize",scale.factor = 10000)
Patient69364 = FindVariableFeatures(Patient69364, selection.method="vst",nfeatures= 2000)
Patient69364 = ScaleData(Patient69364)
Patient69364 = RunPCA(Patient69364, features = VariableFeatures(object = Patient69364))
Patient69364 = FindNeighbors(Patient69364, dims = 1:30)
Patient69364 = FindClusters(Patient69364, resolution = 0.6)
Patient69364 = RunUMAP(Patient69364, dims = 1:30)
sweep.res.list_Patient69364 = paramSweep(Patient69364, PCs = 1:10, sct=FALSE)
sweep.stats_Patient69364 = summarizeSweep(sweep.res.list_Patient69364, GT = FALSE)
bcmvn_Patient69364 = find.pK(sweep.stats_Patient69364)
pK_Patient69364 = bcmvn_Patient69364 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient69364 = as.numeric(as.character(pK_Patient69364[[1]]))
annotations_Patient69364 = Patient69364@meta.data$seurat_clusters
homotypic.prop_Patient69364 = modelHomotypic(annotations_Patient69364)
nExp_poi_Patient69364 = round(0.05*nrow(Patient69364@meta.data))
nExp_poi_adj_Patient69364 = round(nExp_poi_Patient69364*(1-homotypic.prop_Patient69364))
Patient69364 = doubletFinder(Patient69364, PCs = 1:10, pN = 0.25, pK = pK_Patient69364, nExp = nExp_poi_adj_Patient69364, reuse.pANN = FALSE, sct = FALSE)
Patient69364@meta.data$IsDoublet = Patient69364@meta.data[, grep("DF.classifications", names(Patient69364@meta.data), value = TRUE)] != "Singlet"
Patient69364@meta.data$ProbDoublet = Patient69364@meta.data[, grep("pANN", names(Patient69364@meta.data), value = TRUE)]
Patient69364@meta.data = Patient69364@meta.data[, !grepl("DF.classifications|pANN", names(Patient69364@meta.data))]


# 6896-EO-1 sample ####################################################################################################################
sc34=load10X("6896-EO/10x_analysis_6896-EO/Sample_6896-EO-1-GEX_AACCACGCAT-TAACCTGAAT")
sc34 = autoEstCont(sc34)
out34 = adjustCounts(sc34, roundToInt = TRUE)
Patient68961 = CreateSeuratObject(out34)
Patient68961$EdgarID = "Patient68961"
Patient68961$ID = "Maze35"
Patient68961$Cohort = "FSGS"
Patient68961$Sex = "Female"
Patient68961$Age = "62"
Patient68961$APOL_Allele_Number = "N/A"
Patient68961$APOL_Alleles = "N/A"
Patient68961$N264K = "N/A"
Patient68961$eGFR_Bx = "38.416"
Patient68961$UPCR_Bx = "1.065"
Patient68961$percent.mt = PercentageFeatureSet(Patient68961, pattern="^MT-")
Patient68961 = subset(Patient68961, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68961 = NormalizeData(Patient68961, normalization.method="LogNormalize",scale.factor = 10000)
Patient68961 = FindVariableFeatures(Patient68961, selection.method="vst",nfeatures= 2000)
Patient68961 = ScaleData(Patient68961)
Patient68961 = RunPCA(Patient68961, features = VariableFeatures(object = Patient68961))
Patient68961 = FindNeighbors(Patient68961, dims = 1:30)
Patient68961 = FindClusters(Patient68961, resolution = 0.6)
Patient68961 = RunUMAP(Patient68961, dims = 1:30)
sweep.res.list_Patient68961 = paramSweep(Patient68961, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68961 = summarizeSweep(sweep.res.list_Patient68961, GT = FALSE)
bcmvn_Patient68961 = find.pK(sweep.stats_Patient68961)
pK_Patient68961 = bcmvn_Patient68961 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68961 = as.numeric(as.character(pK_Patient68961[[1]]))
annotations_Patient68961 = Patient68961@meta.data$seurat_clusters
homotypic.prop_Patient68961 = modelHomotypic(annotations_Patient68961)
nExp_poi_Patient68961 = round(0.05*nrow(Patient68961@meta.data))
nExp_poi_adj_Patient68961 = round(nExp_poi_Patient68961*(1-homotypic.prop_Patient68961))
Patient68961 = doubletFinder(Patient68961, PCs = 1:10, pN = 0.25, pK = pK_Patient68961, nExp = nExp_poi_adj_Patient68961, reuse.pANN = FALSE, sct = FALSE)
Patient68961@meta.data$IsDoublet = Patient68961@meta.data[, grep("DF.classifications", names(Patient68961@meta.data), value = TRUE)] != "Singlet"
Patient68961@meta.data$ProbDoublet = Patient68961@meta.data[, grep("pANN", names(Patient68961@meta.data), value = TRUE)]
Patient68961@meta.data = Patient68961@meta.data[, !grepl("DF.classifications|pANN", names(Patient68961@meta.data))]


# 6896-EO-2 sample ####################################################################################################################
sc35=load10X("6896-EO/10x_analysis_6896-EO/Sample_6896-EO-2-GEX_CCCACCACAA-AAGCGGAGGT")
sc35 = autoEstCont(sc35)
out35 = adjustCounts(sc35, roundToInt = TRUE)
Patient68962 = CreateSeuratObject(out35)
Patient68962$EdgarID = "Patient68962"
Patient68962$ID = "Maze36"
Patient68962$Cohort = "MCD"
Patient68962$Sex = "Female"
Patient68962$Age = "34"
Patient68962$APOL_Allele_Number = "N/A"
Patient68962$APOL_Alleles = "N/A"
Patient68962$N264K = "N/A"
Patient68962$eGFR_Bx = "68.085"
Patient68962$UPCR_Bx = "10.5435"
Patient68962$percent.mt = PercentageFeatureSet(Patient68962, pattern="^MT-")
Patient68962 = subset(Patient68962, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68962 = NormalizeData(Patient68962, normalization.method="LogNormalize",scale.factor = 10000)
Patient68962 = FindVariableFeatures(Patient68962, selection.method="vst",nfeatures= 2000)
Patient68962 = ScaleData(Patient68962)
Patient68962 = RunPCA(Patient68962, features = VariableFeatures(object = Patient68962))
Patient68962 = FindNeighbors(Patient68962, dims = 1:30)
Patient68962 = FindClusters(Patient68962, resolution = 0.6)
Patient68962 = RunUMAP(Patient68962, dims = 1:30)
sweep.res.list_Patient68962 = paramSweep(Patient68962, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68962 = summarizeSweep(sweep.res.list_Patient68962, GT = FALSE)
bcmvn_Patient68962 = find.pK(sweep.stats_Patient68962)
pK_Patient68962 = bcmvn_Patient68962 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68962 = as.numeric(as.character(pK_Patient68962[[1]]))
annotations_Patient68962 = Patient68962@meta.data$seurat_clusters
homotypic.prop_Patient68962 = modelHomotypic(annotations_Patient68962)
nExp_poi_Patient68962 = round(0.05*nrow(Patient68962@meta.data))
nExp_poi_adj_Patient68962 = round(nExp_poi_Patient68962*(1-homotypic.prop_Patient68962))
Patient68962 = doubletFinder(Patient68962, PCs = 1:10, pN = 0.25, pK = pK_Patient68962, nExp = nExp_poi_adj_Patient68962, reuse.pANN = FALSE, sct = FALSE)
Patient68962@meta.data$IsDoublet = Patient68962@meta.data[, grep("DF.classifications", names(Patient68962@meta.data), value = TRUE)] != "Singlet"
Patient68962@meta.data$ProbDoublet = Patient68962@meta.data[, grep("pANN", names(Patient68962@meta.data), value = TRUE)]
Patient68962@meta.data = Patient68962@meta.data[, !grepl("DF.classifications|pANN", names(Patient68962@meta.data))]


# 6896-EO-3 sample ####################################################################################################################
sc36=load10X("6896-EO/10x_analysis_6896-EO/Sample_6896-EO-3-GEX_GCGCTTATGG-CTAGCCAGGC")
sc36 = autoEstCont(sc36)
out36 = adjustCounts(sc36, roundToInt = TRUE)
Patient68963 = CreateSeuratObject(out36)
Patient68963$EdgarID = "Patient68963"
Patient68963$ID = "Maze37"
Patient68963$Cohort = "FSGS"
Patient68963$Sex = "Female"
Patient68963$Age = "32"
Patient68963$APOL_Allele_Number = "N/A"
Patient68963$APOL_Alleles = "N/A"
Patient68963$N264K = "N/A"
Patient68963$eGFR_Bx = "33.428"
Patient68963$UPCR_Bx = "17.7314"
Patient68963$percent.mt = PercentageFeatureSet(Patient68963, pattern="^MT-")
Patient68963 = subset(Patient68963, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68963 = NormalizeData(Patient68963, normalization.method="LogNormalize",scale.factor = 10000)
Patient68963 = FindVariableFeatures(Patient68963, selection.method="vst",nfeatures= 2000)
Patient68963 = ScaleData(Patient68963)
Patient68963 = RunPCA(Patient68963, features = VariableFeatures(object = Patient68963))
Patient68963 = FindNeighbors(Patient68963, dims = 1:30)
Patient68963 = FindClusters(Patient68963, resolution = 0.6)
Patient68963 = RunUMAP(Patient68963, dims = 1:30)
sweep.res.list_Patient68963 = paramSweep(Patient68963, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68963 = summarizeSweep(sweep.res.list_Patient68963, GT = FALSE)
bcmvn_Patient68963 = find.pK(sweep.stats_Patient68963)
pK_Patient68963 = bcmvn_Patient68963 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68963 = as.numeric(as.character(pK_Patient68963[[1]]))
annotations_Patient68963 = Patient68963@meta.data$seurat_clusters
homotypic.prop_Patient68963 = modelHomotypic(annotations_Patient68963)
nExp_poi_Patient68963 = round(0.05*nrow(Patient68963@meta.data))
nExp_poi_adj_Patient68963 = round(nExp_poi_Patient68963*(1-homotypic.prop_Patient68963))
Patient68963 = doubletFinder(Patient68963, PCs = 1:10, pN = 0.25, pK = pK_Patient68963, nExp = nExp_poi_adj_Patient68963, reuse.pANN = FALSE, sct = FALSE)
Patient68963@meta.data$IsDoublet = Patient68963@meta.data[, grep("DF.classifications", names(Patient68963@meta.data), value = TRUE)] != "Singlet"
Patient68963@meta.data$ProbDoublet = Patient68963@meta.data[, grep("pANN", names(Patient68963@meta.data), value = TRUE)]
Patient68963@meta.data = Patient68963@meta.data[, !grepl("DF.classifications|pANN", names(Patient68963@meta.data))]


# 6896-EO-4 sample ####################################################################################################################
sc37=load10X("6896-EO/10x_analysis_6896-EO/Sample_6896-EO-4-GEX_AGTTTCCTGG-CTGTGTGGCA")
sc37 = autoEstCont(sc37)
out37 = adjustCounts(sc37, roundToInt = TRUE)
Patient68964 = CreateSeuratObject(out37)
Patient68964$EdgarID = "Patient68964"
Patient68964$ID = "Maze38"
Patient68964$Cohort = "FSGS"
Patient68964$Sex = "Female"
Patient68964$Age = "22"
Patient68964$APOL_Allele_Number = "N/A"
Patient68964$APOL_Alleles = "N/A"
Patient68964$N264K = "N/A"
Patient68964$eGFR_Bx = "100.517"
Patient68964$UPCR_Bx = "1.8313"
Patient68964$percent.mt = PercentageFeatureSet(Patient68964, pattern="^MT-")
Patient68964 = subset(Patient68964, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
Patient68964 = NormalizeData(Patient68964, normalization.method="LogNormalize",scale.factor = 10000)
Patient68964 = FindVariableFeatures(Patient68964, selection.method="vst",nfeatures= 2000)
Patient68964 = ScaleData(Patient68964)
Patient68964 = RunPCA(Patient68964, features = VariableFeatures(object = Patient68964))
Patient68964 = FindNeighbors(Patient68964, dims = 1:30)
Patient68964 = FindClusters(Patient68964, resolution = 0.6)
Patient68964 = RunUMAP(Patient68964, dims = 1:30)
sweep.res.list_Patient68964 = paramSweep(Patient68964, PCs = 1:10, sct=FALSE)
sweep.stats_Patient68964 = summarizeSweep(sweep.res.list_Patient68964, GT = FALSE)
bcmvn_Patient68964 = find.pK(sweep.stats_Patient68964)
pK_Patient68964 = bcmvn_Patient68964 %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK_Patient68964 = as.numeric(as.character(pK_Patient68964[[1]]))
annotations_Patient68964 = Patient68964@meta.data$seurat_clusters
homotypic.prop_Patient68964 = modelHomotypic(annotations_Patient68964)
nExp_poi_Patient68964 = round(0.05*nrow(Patient68964@meta.data))
nExp_poi_adj_Patient68964 = round(nExp_poi_Patient68964*(1-homotypic.prop_Patient68964))
Patient68964 = doubletFinder(Patient68964, PCs = 1:10, pN = 0.25, pK = pK_Patient68964, nExp = nExp_poi_adj_Patient68964, reuse.pANN = FALSE, sct = FALSE)
Patient68964@meta.data$IsDoublet = Patient68964@meta.data[, grep("DF.classifications", names(Patient68964@meta.data), value = TRUE)] != "Singlet"
Patient68964@meta.data$ProbDoublet = Patient68964@meta.data[, grep("pANN", names(Patient68964@meta.data), value = TRUE)]
Patient68964@meta.data = Patient68964@meta.data[, !grepl("DF.classifications|pANN", names(Patient68964@meta.data))]


# RPCA integration ####################################################################################################################

library(Seurat)
library(SeuratDisk)
library(scCustomize)
library(tidyverse)
library(patchwork)
library(biomaRt)
library(Matrix)
library(anndata)

setwd("~/workspace/MZ-167-NEPTUNE-scratch/")

AllMaze0 = merge(Patient68011, y=c(Patient68012,Patient68013,Patient68014,Patient68171,Patient68172,
                                   Patient68173,Patient68174,Patient68311,Patient68312,Patient68313,
                                   Patient68314,Patient63341,Patient63342,Patient63481,Patient63482,
                                   Patient67931,Patient67932,Patient67933,Patient68391,Patient68392,
                                   Patient68393,Patient68394,Patient68711,Patient68712,Patient68713,
                                   Patient68861,Patient68862,Patient68863,Patient68864,Patient69361,
                                   Patient69362,Patient69363,Patient69364,Patient68961,Patient68962,
                                   Patient68963,Patient68964),
                 add.cell.ids=c("68011","68012","68013","68014","68171","68172","68173","68174","68311",
                                "68312","68313","68314","63341","63342","63481","63482","67931","67932",
                                "67933","68391","68392","68393","68394","68711","68712","68713","68861",
                                "68862","68863","68864","69361","69362","69363","69364","68961","68962",
                                "68963","68964"),
                 project="MazeIntegrate")

AllMaze0[["RNA"]] = JoinLayers(AllMaze0[["RNA"]])
anndata  = AnnData(
  X = t(AllMaze0@assays$RNA$counts),
  obs = as.data.frame(AllMaze0@meta.data),
  var = data.frame(rownames(AllMaze0@assays$RNA$counts),
    row.names = rownames(AllMaze@assays$RNA$counts))
)
write_h5ad(anndata, filename = "./output/neptune_10x.h5ad")
