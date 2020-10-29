setwd("~/work/CostertonBiofilmCenter/Projects/Current/H6paper_RNAseq")

#Get list of count files

counts <- list.files("./generated","*.counts$", full.names = T)
counts

#What conditions we want to compare 

conditions <- list(c("condition","WT_OD1.2_No","WT_OD1.2_Yes"), 
                   c("condition","WT_OD1.5_No","WT_OD1.5_Yes"),
                   c("condition","WT_OD2.5_No","WT_OD2.5_Yes"),
                   c("condition","wspFpelpsl_OD1.2_No","wspFpelpsl_OD1.2_Yes"),
                   c("condition","wspFpelpsl_OD1.5_No","wspFpelpsl_OD1.5_Yes"),
                   c("condition","wspFpelpsl_OD2.5_No","wspFpelpsl_OD2.5_Yes")
                   )

#First we'll define a function for performing the DESeq analysis.
#
# INPUTS:
#
# count_files = a list of count file paths (see above command)
# metadata = tsv file containing sample metadata
# csv_out = If True, writes out results from DESeq
#

run_DEseq <- function(count_files, metadata, conditions, csv_out=FALSE, out_dir=NULL, vst_out=FALSE){
  
  ###########Build DESeq count data frame######################
  count_data <- lapply(count_files, read.delim,
                       header=T,
                       comment= "#")
  
  #Extract counts and gene names from first count_file
  cds <- data.frame(count_data[[1]][7])
  row.names(cds)<-count_data[[1]]$Geneid
  
  #For every other count file, bind genes to 
  
  for (index in 2:length(count_data)){
    name <- colnames(count_data[[index]])[7]
    cds<-cbind(cds,count_data[[index]][7])
    
  }
  
  #Fix colnames in CDS
  
  colnames(cds)<-gsub("..generated.","",colnames(cds))
  colnames(cds)<-gsub(".aligned.bam","",colnames(cds))
  
  ################## Build Metadata ##########################
  
  metadata <- read.delim(metadata,stringsAsFactors = F)
  
  metadata<- subset(metadata, Sample_Name %in% colnames(cds))
  
  cds <- cds[,order(colnames(cds))]
  metadata <- metadata[order(metadata$Sample_Name),]
  
  
  metadata$Strain<-gsub(" PAO1","",metadata$Strain)
  
  metadata$condition<-paste(metadata$Strain,metadata$Growth_time,metadata$Treatment,sep="_")
  
  
  ############ Run DESEQ2 #####################################
  require(DESeq2)
  

  dds <- DESeqDataSetFromMatrix(countData = cds,
                                colData = metadata,
                                design = ~ condition )
  
  if(isTRUE(vst_out)){
    rld <- vst(dds, blind=FALSE)
    write.csv(assay(rld),paste0(out_dir,"VST_normalized_counts.csv"))
  }
  
  
  
  dds <- DESeq(dds)
  
  lapply(conditions, function(x){
    
    print(x)
    
    res <- results(dds, contrast = x )
    
    if(csv_out==TRUE & !is.null(out_dir)){
    
        print(paste0("Writing out: ", paste0(out_dir,x[2],"VS",x[3],"DiffEXP.csv")))
        print(head(subset(data.frame(res),padj<0.1)))
      
        write.csv(res,paste0(out_dir,x[2],"VS",x[3],"DiffEXP.csv"))
      
      }
  })
}




###############################DEfine and run the analyses######################################################


run_DEseq(count_files = counts, 
          metadata = "./sample_metadata.tsv", 
          conditions = conditions, 
          csv_out=TRUE,
          vst_out = TRUE, 
          out_dir = "./DESeq_results_25May20/")


