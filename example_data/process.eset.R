# function to process gset data for TransApp
# Apr/09/2015
# Abhaydeep Pandey
# Note name of expression set must be eset !!!!
#library("GEOquery")
#gse <- getGEO("GSE26378",GSEMatrix=TRUE)
#gset <- gse[[1]]

library("hgu133plus2.db")
library("genefilter")
require("org.Hs.eg.db")
process.eset <- function(study.id="GSE26440"){
  cat("Loading eset ",study.id,sep = "")
  cat("\n")
  data.file=paste(study.id,".rda")
  if(file.exists(data.file)){
    cat("loadind data....","\n")
    load(data.file) -> a
  } else {
    cat("Downloadinf data, process will take some time....","\n")
    gset <- getGEO(study.id, GSEMatrix =TRUE)
    eset <- gset[[1]]
    cat("Saving data...","\n")
    save(gset,file=paste(study.id ,".rda"))
  }
  cat("data loaded..","\n")
  annotation(eset) <- "org.Hs.eg.db"

  cat("Creating directory and sub-directories")
  cat("\n")
  
  # Create required directories
  folders <- c("metadata","rawdata","processed_data","result")
  dir.create(study.id)
  parent.dir <- study.id
  for (i in 1:length(folders)){
    dir.create(paste(parent.dir,folders[i], sep="/"))
  } 
  
  datafile <- paste(study.id, "/rawdata/",study.id,"_gexp.txt", sep="")
  feat.datafile <- paste(study.id, "/metadata/",study.id,"_feature_annot.csv", sep="")
  cat("Processing data")
  cat("\n")
  ## Remove genes that have no entrezGene id
  entrezIds <- mget(featureNames(eset), envir=hgu133plus2ENTREZID,ifnotfound = NA)
  haveEntrezId <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]
  numNoEntrezId <- length(featureNames(eset)) - length(haveEntrezId)
  cat(class(haveEntrezId),"\n")
  eset <- eset[haveEntrezId, ]
  
  
  ## Non-specific filtering based on IQR
  iqrCutoff <- 0.05
  esetIqr <- apply(exprs(eset), 1, IQR)
  selected <- esetIqr > iqrCutoff
  nsFiltered <- eset[selected, ]
  
  # nsFiltering-unique
  numNsWithDups <- length(featureNames(nsFiltered))
  nsFilteredIqr <- esetIqr[selected]
  uniqGenes <- findLargest(featureNames(nsFiltered), nsFilteredIqr, 
                           "hgu133plus2")
  eset <- nsFiltered[uniqGenes, ]
  e.Ids <- mget(featureNames(eset), envir=hgu133plus2ENTREZID)
  featureNames(eset) <- e.Ids 
  gexp <- exprs(eset) 
  write.table(gexp,file = datafile,quote = FALSE, sep = "\t")
  feat.annot <- cbind(rownames(gexp),rownames(gexp))
  colnames(feat.annot) <- c("Probe_ID","Entrez_ID")
  write.table(feat.annot,file = feat.datafile,quote = FALSE, sep = "\t",row.names = FALSE)
  cat("Process complete for ",study.id,sep = "")
  cat("\n")
}
