#new attempt of RNAseq analysis of the B.subt data based both on Sandra's script and RNA_seq course from September 2018
setwd("~/Work_Docs/Bacillus_RNA_seq/New_B_subt_analysis/")
library(tidyverse)
library("ensembldb")
library("readr")
library("tximport")
library(tibble)
library("DESeq2")
library("ggplot2")

### Import quant.sf files from salmon
dir <- ("~/Work_Docs/Bacillus_RNA_seq/salmon_mapping_results") # set directory to where the data is
samples <-read.table(file.path(dir,"dirlist_n.txt"),header = TRUE) # import tables containing information about the samples
files <- file.path(dir,samples$sample,"quant.sf") # load all quant.sf files
names(files) <- (samples$sample)
all(file.exists(files)) # Check all files exist

### Building annotation databases from gff files using the ensembldb packages
bsbSq1 <- ensDbFromGff("~/Work_Docs/Bacillus_RNA_seq/Data_analysis/Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.38.chromosome.Chromosome.gff3.gz") 
bsDb <- EnsDb(bsbSq1) # building the annotation databse for Bacillus_subtilis_PY79

### Construct tx2gene table (transcript ID to gene ID)

k <- keys (bsDb,keytype = "TXNAME")
tx2gene <- ensembldb::select(bsDb,k,"GENEID","TXNAME")


tx2gene$TXID <- NULL # I don't why there is a third column "TXID", so I deleted it
head(tx2gene)

### Importing the salmon quantification data using the tximport packages
txi <- tximport(files,type="salmon",tx2gene=tx2gene,dropInfReps=TRUE)
names(txi)
head(txi$counts)

#trying to save txi object to transfer to docker container, since dds object transfer was not wery succesful

save(txi, samples, file="Bacillus_RNA_seq_preprocessing_TXI.RData")


txi_counts <- txi$counts


write.csv(txi_counts, file="RNAseq_Counts_Bacillus_new.csv")


dim(txi_counts)
#filtering rows with low count numbers
keep <- rowSums(txi_counts) > 5
countdata <- txi_counts[keep,]
dim(countdata)

#evaluating lib sizes
librarySizes <- colSums(countdata)
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")
abline(h=3e6, lty=2)

# Get log2 counts per million
logcounts <- log2(countdata + 1)


# make a colour vector
genotypeCol <- as.numeric(factor(samples$genotype)) + 1
# Check distributions of samples using boxplots
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=genotypeCol)

# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(as.matrix(logcounts)), col="blue")

#trying to apply rlog transformation of the count data (from rna_seq course)
#need to convert all elements to integers:


countdata_round <- round(countdata) #converting to integers!

rlog_counts <- rlog(countdata_round) #now rlogs will work!

boxplot(rlog_counts, 
        xlab="", 
        ylab="Rlog(Counts)",
        las=2,
        col=genotypeCol)



#trying to return gene names to rlog_counts matrix!
tmp <- logcounts[,0]
rlog_counts <- cbind(tmp, rlog_counts)



#library(Glimma) #not available!

#Hierarchical clustering with heatmaps

# We estimate the variance for each row in the logcounts matrix
countVar <- apply(rlog_counts, 1, var)
# Get the row numbers for the top 500 most variable genes
highVar <- order(countVar, decreasing=TRUE)[1:100]
# Subset logcounts matrix
hmDat <- rlog_counts[highVar,]

library(gplots)
library(RColorBrewer)
# Get some nicer colours
mypalette <- brewer.pal(11, "RdYlBu")
# http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[samples$treatment]
# Plot the heatmap
heatmap.2(hmDat, 
          col=rev(morecols(50)),
          trace="column", 
          main="Top 10 most variable genes across samples",
          ColSideColors=col.cell,scale="row")

#Creating Deseq2 object for analysis
#checking if data is ok
all(samples$sample == colnames(countdata_round))

# create the design formula
design <- as.formula(~ treatment)
# create the DESeqDataSet object
#NB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#This is nor right way to create ddsObj from salmon mapping! see ddseq2 documnetation, or Sndra's file
#need to use DESeqDataSetFromTximport!

#Finally, we can construct a DESeqDataSet from the txi object and sample information in samples.

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ treatment)


#some low counts filtering
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]



#will try to save and transfer dds object to docker container usage

save(dds, samples, file="Bacillus_RNA_seq_preprocessing_dds_object.RData")
# ddsObj <- DESeqDataSetFromMatrix(countData = countdata_round, #this does not apply to Salmon count data!
#                                  colData = samples,
#                                  design = design)

# Apply normalisation to DDS object
ddsObj <- estimateSizeFactors(ddsTxi)


#Take a look at the normalisation factors for these samples.
sizeFactors(ddsObj)

#can not install limma, so the code below does not work!
# normalizedCounts <- counts(ddsObj, normalized=TRUE) 
# logNormalizedCounts <- log2(normalizedCounts + 1)
# par(mfrow=c(1,2))
# limma::plotMA(logNormalizedCounts, array = 7)
# abline(h=0,col="grey")

#We can save a few data objects to use later so we don’t have to rerun everything:
save(countdata, samples, file="RNA_seq_preprocessing.RData")

#Rna_seq_course materials from day 3
#another design with two factors!
# Use the standard R 'formula' syntax for an additive model
design <- as.formula(~ genotype + treatment)

modelMatrix <- model.matrix(design, data = samples)
modelMatrix

# create the DESeqDataSet object
ddsObj.raw <- DESeqDataSetFromMatrix(countData = countdata_round, #need to use countdata round!
                                     colData = samples,
                                     design = design)

#Let’s plot a PCA from vst transformed data. Can you anticipate if the interaction term will be important?
vstcounts <- vst(ddsObj.raw, blind=TRUE)
plotPCA(vstcounts, intgroup=c("genotype", "treatment", "rep"))

#code above does not work for unclear reasons!
# it works if use only one design factor (~treatment)

#The main DESeq2 work flow is carried out in 3 steps:
#First, Calculate the “median ratio” normalisation size factors…

ddsObj <- estimateSizeFactors(ddsObj.raw)

#… then estimate dispersion …
ddsObj <- estimateDispersions(ddsObj)

#… finally, apply Negative Binomial GLM fitting and calculate Wald statistics
ddsObj <- nbinomWaldTest(ddsObj)


#We can generate a table of differential expression results from the DDS object using the results function of DESeq2.
res <- results(ddsObj, alpha=0.05)
head(res)

#By default, results has returned the contrast encoded by the final column in the model matrix. DESeq2 has the command resultsNames that allows us to view the contrasts that are available directly from the model matrix.

resultsNames(ddsObj)

#Let’s get the top 100 genes by adjusted p-value

topGenes_LC_vs_B <- as.data.frame(res) %>%
  rownames_to_column("GeneID") %>% 
  arrange(padj) %>% 
  head(100)
topGenes_LC_vs_B
#topGenes_LC_vs_B_100 <- dplyr::select(topGenes_LC_vs_B, GeneID, log2FoldChange)


#just trying to make a heatmap of top diff genes based on prevous heatmap
#
#nihera ne vishlo!



#Finally save the results in a new RData object
save(res, ddsObj, samples, file="B_subt_DE.RData")

#day3_part2
# visualisation and annotation

library(biomaRt)


listMarts()
listMarts(mart = NULL, host="http://useast.ensembl.org", path="/biomart/martservice",
          port=80, includeHosts = FALSE, archive = FALSE, ssl.verifypeer = TRUE, 
          ensemblRedirect = NULL, verbose = FALSE)


## set up connection to ensembl database
ensembl=useMart("ENSEMBL_MART_ENSEMBL")

listDatasets(ensembl)

# DOES NOT HAVE ANY BACTERIAL DATA! nEED TO FIND DIFFERENT WAY FOR ANNOTATIONS



#DESeq2 provides a functon called lfcShrink that shrinks log-Fold Change (LFC) estimates towards zero using and empirical Bayes procedure. 

ddsShrink <- lfcShrink(ddsObj, coef="treatment_LC_vs_B")
shrink_LC_vs_B <- as.data.frame(ddsShrink) %>%
  rownames_to_column("GeneID") %>% 
  rename(logFC=log2FoldChange, FDR=padj)

#P-value histogram

hist(shrink_LC_vs_B$pvalue)


#MA plots are a common way to visualize the results of a differential analysis. We met them briefly towards the end of Session 2. This plot shows the log-Fold Change for each gene against its average expression across all samples in the two conditions being contrasted.


plotMA(ddsShrink, alpha=0.05)




# add a column with the names of only the top 10 genes
cutoff <- sort(shrink_LC_vs_B$pvalue)[10]
shrink_LC_vs_B <- shrink_LC_vs_B %>% 
  mutate(TopGeneLabel=ifelse(pvalue<=cutoff, GeneID, ""))
ggplot(shrink_LC_vs_B, aes(x = log2(baseMean), y=logFC)) + 
  geom_point(aes(colour=FDR < 0.05), shape=20, size=0.5) +
  geom_text(aes(label=TopGeneLabel)) +
  labs(x="mean of normalised counts", y="log fold change")



#probably need to install docker, since mane of the libaries do not work here!








### Construct a DESeqDataSet from the txi object
ddsTxi <- DESeqDataSetFromTximport(txi,colData=samples, design = ~ genotype+treatment)

### Set the factors level (control group)
ddsTxi$genotype <- relevel(ddsTxi$genotype, ref="g021")
ddsTxi$treatment <- relevel(ddsTxi$treatment, ref="B")
ddsTxi$rep <- factor(ddsTxi$rep, levels=c("rep1","rep2"))


### Differential epxression analysis
dds <- DESeq(ddsTxi)
