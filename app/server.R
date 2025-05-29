# Note : log-transform data date: 30/apr/2015
library(shiny)
library(shinyIncubator)
knitr::include_graphics("images/fig31_adjust_nodes.png")
#source("trans_pipe.R")
# Define server logic required to summarize and view the selected dataset
shinyServer(function(input, output, session){
  # Return the requested dataset
  # get the tech_id
  # The main script of the automated pipeline for transcriptomic data analysis
  output$trans_pipe <- renderPrint({ # renderText({ #
    techid <- input$platform
    projid <-  input$proj_id
    if(techid != 0 && projid !="GSEXXXX"){
      withProgress(session, min=1, max=15, expr={
        for(i in 1:15) {
          setProgress(message = "Loading required packages" ,
                      value=i)
          Sys.sleep(0.3)
        }
      })
      require("GEOquery")
      require("genefilter")
      require("GSEABase")
      require("KEGG.db")
      require("Category")
      require("lattice")
      require("SPIA")
      require("annotate")
      require("GOstats")
      require("RColorBrewer")
      require("xtable")
      require("Rgraphviz")
      require("hgu133plus2.db") # techid 1
      require("hgu133plus2cdf")
      require("hgu133plus2probe")
      require("illuminaHumanv4.db") # techid 2
      require("limma")
      require("Biobase")
      require("gcrma")
      require("arrayQualityMetrics")
      require("affy") # techid 1
      require("ReportingTools")
      require("hwriter")
      require("KEGGprofile")
      require("DOSE")
      #require("ReactomePA")
      options(digits = 3)
      withProgress(session, min=1, max=15, expr={
        for(i in 1:15) {
          setProgress(message = "Starting the process ..." ,
                      value=i)
          Sys.sleep(0.3)
        }
      })
      
      trim <- function (x) gsub("^\\s+|\\s+$", "", x) # to remove leading or trailing whitespace
      
      sink(paste(paste(projid, "/result/", sep=""),"output.log",sep=''))

      # read the paths
      rawdatapath = paste(projid, "/raw_data/", sep="")
      targetspath = paste(projid, "/metadata/Targets.txt", sep="")
      processeddatapath = paste(projid, "/processed_data/", sep="")
      qcpath = paste(projid, "/QC/", sep="")
      resultpath = paste(projid, "/result/", sep="")
      
      # Read the targets; Targets.txt
      # This file should have atleast two columns: Name  Group  Celfile
      # For Affymetrix (techid 1) there is also Celfile column
      
      # Read the target files for both techids we can read it by readtargets()
      targets <- readTargets(file=targetspath) 
      if(techid == 1) {
        ###########################################################################
        ## For affy data                                                        ###
        ###########################################################################
        eset.rda_file_name <- paste(processeddatapath, projid, "_eset.rda", sep="")
        if(file.exists(eset.rda_file_name)) {
          cat(paste("Reading the file ", eset.rda_file_name, " ...", sep=""))
          load(file=eset.rda_file_name)
          cat("\n")
        } else {
          # Run QC ## from .cel files to expression set
          myCovs = data.frame(targets)
          rownames(myCovs) = myCovs[,1]
          nlev = as.numeric(apply(myCovs, 2, function(x) nlevels(as.factor(x))))
          metadata = data.frame(labelDescription = paste(colnames(myCovs), ": ",
                                                         nlev, " level", ifelse(nlev==1,"","s"), sep=""), 
                                row.names=colnames (myCovs))
          phenoData = new("AnnotatedDataFrame", data=myCovs, varMetadata=metadata)
          dat = ReadAffy(sampleNames=myCovs$Name, 
                         filenames=myCovs$Celfile, 
                         phenoData=phenoData,               
                         celfile.path=rawdatapath)
          
          # Quality assessment. Store the results into a directory called QC.
          cat("QC running")
          # arrayQualityMetrics(expressionset=dat, outdir=qcpath,
          #                     do.logtransform=TRUE, intgroup="Group",
          #                    reporttitle=paste("Quality Control report for", projid,sep=""),force=TRUE)
          
          # Normalization (Background Adjustment Using Sequence Information) 
          # using gcrma; save eset in processed_data file
          eset = gcrma(dat)
          annotation(eset) <- "hgu133plus2.db"
          
          # Non-specific filtering
          eset_filt <- nsFilter(eset, var.filter=F)$eset
          eset <- eset_filt
          featureNames(eset) <- as.character(unlist(mget(featureNames(eset), hgu133plus2ENTREZID)))
          annotation(eset) <- "org.Hs.eg.db"
          
          save(eset, file=eset.rda_file_name)
        }
        
      } else if(techid==2){ ## Normalized data , no need to do normalization again
        ###########################################################################
        ## For illumina data                                                    ###
        ###########################################################################
        eset.rda_file_name <- paste(processeddatapath, projid, "_eset.rda", sep="")
        if(file.exists(eset.rda_file_name)) {
          cat(paste("Reading the file ", eset.rda_file_name, " ...", sep=""))
          load(file=eset.rda_file_name)
          cat("\n")
        } else {
          #  gse <- getGEO(projid,GSEMatrix=TRUE)
          #  gse.ex <- gse[[1]]
          #  annotation(gse.ex) <- "illuminaHumanv4.db"
          #  eset_filt <- nsFilter(gse.ex, var.filter=F)$eset
          #  eset <- eset_filt[,targets$Sample.ID] # selected from targets only
          #  eset$Group  <- targets$Group 
          #  eng.ids <- as.character(unlist(mget(featureNames(eset), illuminaHumanv4ENTREZID)))
          #  save(eset, file=eset.rda_file_name)
          cat("File not found !!" )
        }
        #Quality assessment. Store the results into a directory called QC.
        #cat("QC running")
        #arrayQualityMetrics(expressionset=eset, outdir=qcpath, do.logtransform=TRUE, reporttitle=paste("Quality Control report for", projid,sep=""),force=TRUE)
      }else if(techid==3) {
        # The minimum data provided by user will have matrix with rownames as Probe id, column names as sample_id, and expression values
        ###########################################################################
        ## For customized data                                                        ###
        ###########################################################################
        eset.rda_file_name <- paste(processeddatapath, projid, "_eset.rda", sep="")
	deresultfile <- paste(resultpath, paste(projid,"All_genes.csv", sep=''), sep='')
        if(file.exists(eset.rda_file_name)) {
          cat(paste("Reading the file ", eset.rda_file_name, " ...", sep=""))
          load(file=eset.rda_file_name)
          cat("\n")
        } else {
          withProgress(session, min=1, max=15, expr={
            for(i in 1:15) {
              setProgress(message = "Preparing data for analysis" ,
                          value=i)
              Sys.sleep(0.3)
            }
          })
          targets <- readTargets(file=targetspath) 
          datafile <- paste(projid, "/rawdata/",projid,"_gexp.txt", sep="")
          feat.datafile <- paste(projid, "/metadata/",projid,"_feature_annot.csv", sep="") # manually delete the un usefull ids  
          expr <- as.matrix(read.table(file=datafile, header=TRUE,sep="\t", row.names=1, as.is=TRUE))
          feat_annot <- as.matrix(read.table(file=feat.datafile, header=TRUE,sep="\t", as.is=TRUE))
          rownames(feat_annot) <- feat_annot[,"Probe_ID"]
          
          # Get selected samples on the basis of targets
          expr <- expr[,trim(targets$Sample.ID)]
          pData <- targets[,c("Sample.ID","Group")]
          rownames(pData)<-trim(targets$Sample.ID)
          
          # rownames(expr)<-rownames(as.matrix(expr))
          metadata = data.frame(labelDescription=c("Patient ID","Case/control status"),row.names=c("Sample_ID", "Group"))
          adf = new("AnnotatedDataFrame", data=pData,varMetadata=metadata)
          eset=new("ExpressionSet",exprs=expr,phenoData=adf)

          # Find the probes with non-missing values for at least 3 samples
          nonmissing <- apply(exprs(eset), 1, function(x) sum(!is.na(x))>=3)
          eset <- eset[nonmissing, ]
          
          # Remove duplicate entrez ids
          # a.na <- feat_annot[sel,"Entrez_ID"]
          a.na <- feat_annot[featureNames(eset),"Entrez_ID"]
          
          # Entrez ids without na
          a.na <- a.na[-grep("//", a.na)]
          a <- a.na[!is.na(a.na)]
          a <- trim(a)
          eg.nonmissing <- a
          if(length(a)>0){
            eset <- eset[names(a),]
          } 
          b <- apply(exprs(eset),1,IQR,na.rm=T)
          # Get IQR value without na
          b <- b[!is.na(b)]
          #b <- b[!is.na(a.na)]
          
          # creating named vector
          names(b) <- eg.nonmissing 
          cx <- unique(eg.nonmissing)
          
          maxin <- 0
          if(length(cx) != length(eg.nonmissing)){
            # Get index of single probes and probes with max IQR if duplicated
            for(i in 1:length(cx)){
              index<-which(cx[i]==eg.nonmissing)
              maxin[i]<-index[which.max(b[index])]
            }
          }
          
          # Get final e-set object 
          if(maxin > 0){
            eset <- eset[maxin,]
            featureNames(eset) <- eg.nonmissing[maxin]
          }
          
          null.idx <- which(featureNames(eset)=="")
          if(length(null.idx) > 0){
            eset <- eset[-null.idx,]
          }
          annotation(eset) <- "org.Hs.eg.db"
          Grp <- factor(eset$Group)
          save(eset, file=eset.rda_file_name)
        }
      }
          #################################################
	  # log-transform data                #
	  #################################################	
	if(max(exprs(eset))>20){
	withProgress(session, min=1, max=15, expr={
        	for(i in 1:15) {
        	  setProgress(message = paste("Log-transforming ", projid, "; Max intensity: ", max(exprs(eset)), "...\n", sep=""),
        	              detail = paste("Process may take a while"),
        	              value=i)
        	 # print(i)
        	  Sys.sleep(0.3)
        	}
      	})
#	    cat(paste("Log-transforming ", study.id, "; Max intensity: ", max(exprs(eset)), "...\n", sep=""))

	    edat <- exprs(eset)
	    edat <- apply(edat, 2, function(x) {
	      maxx <- max(x, na.rm=T)
	      if(maxx>20.0) { # log-transform data
	        minx <- min(x, na.rm=T)
	        if(minx < 1) {
	          x <- x+1-minx
	        }
	      }
	      log2(x)
	    })
	    exprs(eset) <- edat
	  }
	  else{
		withProgress(session, min=1, max=15, expr={
		        	for(i in 1:15) {
		        	  setProgress(message = paste("Log-transformed data "),
		        	              value=i)
		#	    print.noquote("Log-transformed data ")
			Sys.sleep(0.3)
		        	}
		      	})
	 }

      # Add group information, Group contains the two groups to be compared. Convert this to a factor.
      eset$Group[which(eset$Group != "Zcontrol")] <- "Case"
      Grp <- factor(eset$Group)
      # run model.matrix to get the design matrix,
      # des = model.matrix(~grp)
      des <- model.matrix(~0+Grp) ###
      colnames(des) <-levels(Grp)  
      
      # lmFit for fitting linear model.
      f1 = lmFit(eset, des)
      my.level <- levels(Grp)
      #print(paste(levels(Grp), collapse="vs"),quote=FALSE)
      contrast.matrix <- makeContrasts(CasevsZcontrol = Case-Zcontrol, levels=des)
      ef1 <- contrasts.fit(f1, contrast.matrix)
      
      # eBayes for regularization of t-statistics.
      ef1 = eBayes(ef1)
      
      # Adding Entrez_ID column to ef1 
      ef1$Entrez_Id <- featureNames(eset) #eng.ids
      
      #topTable for getting the most differentially expressed genes.
      #topTable(ef1, 2)
      # Processing to get allgenes.csv file
      top_table1 <- topTable(ef1, number=length(rownames(ef1)))
      top0 <- featureNames(eset) #rownames(top_table1)  # Entrez Ids
      top1 <- formatC(top_table1[,"logFC"],digits=3) # fold change
      top2 <- formatC(top_table1[,"AveExpr"],digits=3) # AveExpr
      top4 <- formatC(top_table1[,"P.Value"],digits=3) # p-value
      top5 <- formatC(top_table1[,"adj.P.Val"],digits=3) # adj.P.Val
      top6 <- as.character(unlist(mget(featureNames(eset), org.Hs.egGENENAME,ifnotfound=NA))) #gene_name_sel   
      top7 <- as.character(unlist(mget(featureNames(eset), org.Hs.egSYMBOL,ifnotfound=NA))) #Gene SYMBOL  
      
      top_table_all <- cbind( top7, top0, top1, top2, top4, top5, top6)
      colnames(top_table_all) <- c("Gene SYMBOL","EntrezID","Log fold change","AveExpr", "p-value","adj.P.Val", "Gene name")
      no.na <- which(is.na(top0)==FALSE) #remove NA
      top_table_all <- top_table_all[no.na,]
      withProgress(session, min=1, max=15, expr={
        for(i in 1:15) {
          setProgress(message = "Detecting differentially expressed genes",
                      detail = paste("Writing data to", deresultfile, sep=""),
                      value=i)
         # print(i)
          Sys.sleep(0.3)
        }
      })
      write.table(top_table_all,file=deresultfile, col.names=T, row.names=F, quote=F, sep="\t")
      
      ################################################################
      ## Pathway analysis                                           ##
      ################################################################
      
      # Apply rowttests on gset
      rtt <- rowttests(eset, "Group")
      
      ###########################################################################
      # Over-representation analysis                                           ##
      ###########################################################################
      print.noquote("Over-representation analysis: kegg_params")
      withProgress(session, min=1, max=15, expr={
        for(i in 1:15) {
          setProgress(message = 'Detecting up-regulated pathways',value=i)
          Sys.sleep(0.3)
        }
      })

      # the gene universe
      universeGeneIds <- featureNames(eset)
      
      lfc<-rtt$dm
      names(lfc) <- featureNames(eset)
      
      # select the genes of interest, i.e., differentially expressed (upregulated)
      de_egs.up <- rownames(rtt)[which(rtt$p.value<0.01 & rtt$dm>0)] # which.up
      
      # select the genes of interest, i.e., differentially expressed (down regulated)
      de_egs.down <- rownames(rtt)[which(rtt$p.value<0.01 & rtt$dm<0)] # which.down
      
      #### KEGGHyperG test for Up regulated Pathways
      hgCutoff=0.001
      kegg_params.up <- new("KEGGHyperGParams",geneIds=de_egs.up, 
                            universeGeneIds=universeGeneIds, annotation=annotation(eset),
                            pvalueCutoff=hgCutoff,testDirection="over")
      hgOver.up = hyperGTest(kegg_params.up)
hgOver.up.1 <- hgOver.up 
      hgOver.up <- summary(hgOver.up)[summary(hgOver.up)$Count > 5,]

      withProgress(session, min=1, max=15, expr={
        for(i in 1:15) {
          setProgress(message = 'Detecting down-regulated pathways',value=i)
         
          Sys.sleep(0.3)
        }
      })
      
      #### KEGGHyperG test for Down regulated Pathways
      kegg_params.down <- new("KEGGHyperGParams",geneIds=de_egs.down, 
                              universeGeneIds=universeGeneIds, annotation=annotation(eset),
                              pvalueCutoff=hgCutoff,testDirection="under")
      hgOver.down = hyperGTest(kegg_params.down)
      hgOver.down.1 <- hgOver.down 
      hgOver.down <- summary(hgOver.down)[summary(hgOver.down)$Count > 5,]
      mixed_path <- intersect((hgOver.up)[,"KEGGID"],(hgOver.down)[,"KEGGID"])
      
      
      ## Processing top_table to generate html file 
      top_table2 <- top_table_all[1:20,]
      
      ## Processing with keggProfile . 
      rtt_sel <- rtt[,c("dm","p.value")]
      
      rownames(rtt_sel) <- names(lfc)
      colnames(rtt_sel) <- c(paste(projid,"_lfc",sep=""),paste(projid,"_p.value",sep=""))
      rtt_lfc <- as.matrix(rtt[,"dm"])
      rownames(rtt_lfc) <- names(lfc)
      
      # Generate color matrix
      par(mar = rep(2, 4))
      min.range <- min(rtt_lfc)
      max.range <- max(rtt_lfc)
      if(min.range > 0){
        min.range <- 0
      }
      if(max.range < 0){
        max.range <- 0
      }
      col <- col_by_value(rtt_lfc,col = c("green","red"),  range = c(-1,1)) #range = c(min.range,max.range)) #
      #dev.off()
      workdir <- getwd()     
      dev.off()     
      setwd(resultpath) # I got no option to set the output.dir of plot_pathway default is current dir. We want itin Result folder.

      ################################################################
      ## Heatmap (sample level) 
      ################################################################
      png(filename=paste("heatmap.png",sep=""))
      col.var <- eset$Group
      col.var <- replace(eset$Group, which(eset$Group == "Zcontrol"), "green")
      col.var <- replace(col.var, which(col.var != "green"), "red")
      vars <- apply(exprs(eset), 1, var, na.rm=T)
      which.varying <- order(vars, decreasing=T)[1:500]

      #### temp to avoid blak heatmap in working directory
      temp <- heatmap(x=exprs(eset[which.varying,]),na.rm=T,
                      main=paste("Heatmap (using 500 genes with highest variance)", sep=""),
                      labRow=NA, ylab="Entrez_ID", xlab = "Sample IDs",margins = c(10, 5),Colv=NA,ColSideColors = col.var)
      dev.off()
      ########## heatmap <- recordPlot() ##############################
     
      withProgress(session, min=1, max=10, expr={
        for(i in 1:10) {
          setProgress(message = 'Generating the heatmap',value=i)
         # print(i)
          Sys.sleep(0.3)
        }
      })
      numSigKeggUp <- sum(hgOver.up[,"Pvalue"]<hgCutoff)
      numSigKeggDown <- sum(hgOver.down[,"Pvalue"]<hgCutoff)
      
      # Top up regulaed keggids
      if(sum(hgOver.up[,"Pvalue"]<hgCutoff)>1 || sum(hgOver.down[,"Pvalue"]<hgCutoff)>1 ) {
        kegglst.up <- (hgOver.up)[,"KEGGID"]
        for(KEGGID in kegglst.up) {
          # download_KEGGfile(pathway_id = KEGGID, species = "hsa") 
          de_genes.up <- geneIdsByCategory(hgOver.up.1)[[KEGGID]]
          temp <- plot_pathway(as.matrix(rtt_sel)[de_genes.up,], type = "bg", bg_col = col, 
                               text_col = "black", magnify = 1.2, species = "hsa", pathway_id = KEGGID)
          temp2 <- write.table(de_genes.up, file=paste0("de", KEGGID, "_genes.txt"), row.names = FALSE, col.names=FALSE)
        }  
       # withProgress(session, min=1, max=10, expr={
          #for(i in 1:10) {
            #setProgress(message = 'Pathway analysis Running',
            #            detail = 'This may take a while...',
            #            value=i)
           # print(i)
          #  Sys.sleep(0.3)
          #}
        #})
        # Top down regulated keggids
        kegglst.down <- (hgOver.down)[,"KEGGID"]
        for(KEGGID in kegglst.down) {
          #download_KEGGfile(pathway_id = KEGGID, species = "hsa") 
          de_genes.down <- geneIdsByCategory(hgOver.down.1)[[KEGGID]]
          temp <- plot_pathway(as.matrix(rtt_sel)[de_genes.down,], type = "bg", bg_col = col, 
                               text_col = "black", magnify = 1.2, species = "hsa", pathway_id = KEGGID)
          temp2 <- write.table(de_genes.down, file=paste0("de", KEGGID, "_genes.txt"),row.names = FALSE, col.names=FALSE)
        }  
        
        ### Processing for html page 
        setwd(workdir)
        my.mat.up <- (hgOver.up)
        rownames(my.mat.up) <- kegglst.up
        
        my.mat.down <- (hgOver.down)
        rownames(my.mat.down) <- kegglst.down
        
        if(length(mixed_path)>0){
          my.mat.up<- my.mat.up[!rownames(my.mat.up) %in% mixed_path,]
          my.mat.down<- my.mat.down[!rownames(my.mat.down) %in% mixed_path,]
          my.df <- rbind(my.mat.up,my.mat.down)
        } else {
          my.df <- rbind(my.mat.up,my.mat.down)
        }
        #withProgress(session, min=1, max=15, expr={
         # for(i in 1:15) {
            #setProgress(message = 'Calculation in progress',
            #            detail = 'This may take a while...',
            #            value=i)
            
           # Sys.sleep(0.3)
          #}
       # })
        # rownames(my.df) <- c(kegglst.up,kegglst.down)
        my.df <- subset(my.df,select = -c(OddsRatio,ExpCount))
        count.sel <- which(my.df[,"Count"] > 10)
        my.df <- my.df[count.sel,]
        
        # Get top kegg id 
        kegglst <- my.df[,"KEGGID"]
        
        up.reg.path <- intersect(kegglst,(hgOver.up)[,"KEGGID"])
        down.reg.path <- intersect(kegglst,(hgOver.down)[,"KEGGID"])
        
        imagename <- c()
        de_genes_link <- c()   
        imagelink <- c()
        
        for(KEGGID in kegglst) {
          imagename[KEGGID] <- as.character(paste0("hsa", KEGGID))
          imagelink[KEGGID]  <- as.character(paste0("../result/","hsa", KEGGID, "_profile_bg.png"))
          de_genes_link[KEGGID] <- as.character(paste0("../result/","de", KEGGID, "_genes.txt"))
        }
        
        my.df$Mapped_Image <- hwriteImage(imagelink, link = imagelink, table = FALSE, width=100, height=80)
        my.df$KEGGID <- hwrite(as.character(my.df$KEGGID), link = paste("http://www.genome.jp/kegg-bin/show_pathway?hsa",
                                                                        as.character(my.df$KEGGID), sep = ''), table=FALSE)
        my.df$Count <- hwrite(as.character(my.df$Count), link = de_genes_link, table = FALSE)
        my.df$Direction <-c(rep("Up",length(up.reg.path)),rep("Down",length(down.reg.path)))
      } else{
        setwd(workdir)
      }
      
      m.d <- as.data.frame(top_table2)
      m.d$EntrezID <- hwrite(as.character(m.d$EntrezID), 
                             link = paste("http://www.ncbi.nlm.nih.gov/gene/",
                                          as.character(m.d$EntrezID), sep = ''), table=FALSE)
      
      imglink  <- as.character(paste0("../result/","heatmap.png",sep=""))
      img <- hwriteImage(imglink, link = imglink, table = FALSE, width=200, height=180)

      withProgress(session, min=1, max=15, expr={
        for(i in 1:15) {
          setProgress(message = 'Processing data for HTML page ',
                      detail = 'Almost Done...',
                      value=i)
          Sys.sleep(0.3)
        }
      })
      
      cat ("processing for HTML report...")
      ######################################################################
      ### html page for analysis report                                   ##
      ######################################################################
      
      html_report <- HTMLReport(shortName = paste("Gene Expression Data Analysis Report"), 
                                title = paste("Gene Expression Data Analysis Report"),
                                reportDirectory = paste("./",projid,"/report",sep=""))
      
      publish(hwrite(paste("Project ID : ", projid ,sep=''),style='font-size:150%',br=TRUE), html_report)
      publish(" ", html_report)
      if(techid == 1){
        publish(hwrite("Array : Affymetrix Human Genome U133 Plus 2.0 Array  " ,style='font-size:150%',br=TRUE), html_report)
      }else if (techid == 2){
        publish(hwrite("Array : Illumina HumanHT-12 V4.0 expression beadchip " ,style='font-size:150%',br=TRUE), html_report)
      }else if (techid == 3 || techid == 4){
        publish(hwrite("Array : Customised" ,style='font-size:150%',br=TRUE), html_report)
      }
      publish(" ", html_report)
      publish(hwrite("Targets File ", link = paste("../metadata/Targets.txt",sep=''),style='font-size:150%',br=TRUE), html_report)
      publish(" ", html_report)
       if(techid != 3){
        publish(" ", html_report)
        publish(hwrite("QC report ", 
                       link = paste("../QC/index.html",sep=''),
                       style='font-size:150%',br=TRUE), html_report)
      
      publish(" ", html_report)
     publish("HTML report of quality metrics about dataset. The quality metrics are mainly on the per array level, i. e. they can be used to assess the relative quality of different arrays within a dataset. ", html_report)
        }
      publish(" ", html_report)
      publish(hwrite(
        paste("Heatmap of top 500 genes",sep=""),
        style='font-size:150%',br=TRUE), 
              html_report, name = "Heatmap",width=100, height=100)
      publish(" ", html_report)
      publish(img, html_report, name = paste("heatmap.png",sep=""))
      publish(" ", html_report)
      publish(hwrite(paste("Top 20 Differentially expressed genes"),style='font-size:150%',br=TRUE), html_report)
      publish(as.data.frame(m.d), html_report, name = "Differencially expressed genes")
      publish(" ", html_report)
      publish(Link("Click here to get list of all Differencially expressed genes ", paste("../result/", paste(projid,"All_genes.csv", sep=''),sep='')), html_report)
      publish(" ", html_report)
      publish("The column headers are self-explanatory: Gene SYMBOL, EntrezID, Log fold change, AveExpr (Average Gene Expression), p-value, adj.P.Val (Adjusted P-value), Gene name", html_report)
      publish(" ", html_report)
      publish(hwrite(paste(" Pathway Analysis Using KEGG Database ",sep=''),style='font-size:150%',br=TRUE), html_report)
      publish(" ", html_report)
      publish("Explanation of the column header: (1) KEGGID - KEGG Pathway ID, links to www.genome.jp, (2) Pvalue - Enrichment p-value from Fisher's Exact Test: the lower the better, (3) Count - Number of genes in the pathway that are differentially expressed, links to a file with the list of Entrez IDs, (4) Size - Total number of genes in the pathway, (5) Term - Pathway Name, (6) Mapped_Image - Image of the pathway with the nodes painted red (if up-regulated) or green (if down-regulated), (7) Direction - Whether the differentially expressed genes are up- or down-regulated. The number of coloured nodes may be less than the number reported under Count, due to multiple Entrez ids mapping to the same KEGG node.", html_report)        
  if(sum(hgOver.up[,"Pvalue"]<hgCutoff)>0 || sum(hgOver.down[,"Pvalue"]<hgCutoff)>0) {
    publish(as.data.frame(my.df), html_report)
  } else {
    publish(hwrite(paste0("No Pathway found to be significantly moderated"),style='font-size:150%',br=TRUE), html_report)
  }
  publish(" ", html_report)
  publish(hwrite("Session log ", 
               link = paste("../result/output.log",sep=''), 
               style='font-size:150%') , html_report)
  publish(" ", html_report)
  publish("The link above provides information of interest to programmers. ", html_report)        
  publish(" ", html_report)
publish("___________________________________________________________________________________________________________", html_report)
publish("This report was generated with the program 'TransApp' developed by Abhaydeep Pandey and Saroj Kant Mohapatra at the National Institute of Biomedical Genomics, Kalyani, India.", html_report)        
  finish(html_report)
      withProgress(session, min=1, max=10, expr={
        for(i in 1:10) {
          setProgress(message = 'Opening the project report',
                      
                      value=i)
         
          Sys.sleep(0.3)
        }
      })
  browseURL(finish(html_report))
  print(sessionInfo())
  sink()
  rm(list=ls())
  gc()
    }
  })
})

