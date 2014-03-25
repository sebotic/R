#copyright Sebastian Burgstaller-Muehlbacher 2012
#this R code extracts CNVs with zero overlap with DGV from 
#Affymetrix Genotyping Console CNV raw data and generates a list with overlapping genes
#

library("RSQLite")
library("xlsx")


getCNVwZeroOL <- function(index){
  file_path <- paste("~/data/sample_",index,".tsv", sep="")
  zerooverlap <- read.csv(file_path, header=T, sep="\t",stringsAsFactors=F)
  zerooverlap <- subset(zerooverlap, zerooverlap[,4] != "Y" & zerooverlap[,10] == 0)
  
  return(zerooverlap)
}

finemapCNV <- function(chr, cnv_from, cnv_to, gene_list){
  if(chr == "X")
    chr <- 23
  
  if(length(gene_list) == 0){
    return("intergenic")
  }else if (length(gene_list) > 1){
    return("several genes affected")
  }else if(length(gene_list) == 1){
    chr_genes <- switch(as.integer(chr),chr1_refseq,chr2_refseq,chr3_refseq,chr4_refseq,chr5_refseq,chr6_refseq,chr7_refseq,chr8_refseq,chr9_refseq,chr10_refseq,chr11_refseq,chr12_refseq,chr13_refseq,chr14_refseq,chr15_refseq,chr16_refseq,chr17_refseq,chr18_refseq,chr19_refseq,chr20_refseq,chr21_refseq,chr22_refseq,chrX_refseq)
    gene_of_interest <- subset(chr_genes, gene_list[1] == chr_genes[,1])
    #cat(gene_of_interest[1,10],"\n")
    #cat(gene_of_interest[1,11],"\n")
    exon_starts <- unlist(strsplit(gene_of_interest[1,10],","))
    exon_ends <- unlist(strsplit(gene_of_interest[1,11],","))
    
    status <- ""
    for(i in 1:length(exon_starts)){
      start <- as.integer(exon_starts[i])
      end <- as.integer(exon_ends[i])
      if((cnv_from <= start & cnv_to <= end & cnv_to >= start) | (cnv_from >= start & cnv_to >= end & cnv_from <= end) | (cnv_from >= start & cnv_to <= end) | (cnv_from <= start & cnv_to >= end)){
        status <- paste0(status,"exon ",i, ",", collapse="")
      }#else status <- "intron"
    }
    
    if(status == ""){
      status = "intron"
    }
    return(status)
  }
}

#takes a genes list for a certain chromosome and returns a list of genes affected by a CNV
mapToGene <- function(chr, cnv_from, cnv_to){
  if(chr == "X")
    chr <- 23
  chr_genes <- switch(as.integer(chr),chr1_refseq,chr2_refseq,chr3_refseq,chr4_refseq,chr5_refseq,chr6_refseq,chr7_refseq,chr8_refseq,chr9_refseq,chr10_refseq,chr11_refseq,chr12_refseq,chr13_refseq,chr14_refseq,chr15_refseq,chr16_refseq,chr17_refseq,chr18_refseq,chr19_refseq,chr20_refseq,chr21_refseq,chr22_refseq,chrX_refseq)
  
  genes_in_region <- c()
  for(i in 1:nrow(chr_genes)){
    if((cnv_from <= chr_genes[i,5] & cnv_to <= chr_genes[i,6] & cnv_to >= chr_genes[i,5]) | (cnv_from >= chr_genes[i,5] & cnv_to >= chr_genes[i,6] & cnv_from <= chr_genes[i,6]) | (cnv_from >= chr_genes[i,5] & cnv_to <= chr_genes[i,6]) | (cnv_from <= chr_genes[i,5] & cnv_to >= chr_genes[i,6])){
      #genes_in_region <- c(genes_in_region, paste(chr_genes[i,1],"(",chr_genes[i,2],")"))
      genes_in_region <- c(genes_in_region,chr_genes[i,1])
    }
    
  }
  
  
  return(unique(genes_in_region))
}


annot_file <- "~/data/NetAffxGenomicAnnotations.Homo_sapiens.hg19.na32.1.db"

con <- dbConnect(drv="SQLite", dbname=annot_file)
refseq <- dbGetQuery(conn=con, statement=paste("SELECT * FROM '", "browse_refseq", "'", sep=""))

#split into single chromosomes, to reduce number of comparisons. bad code, replace with assign(var.name, value) and regex.
chr1_refseq <- subset(refseq,refseq[,3] == "chr1")
chr2_refseq <- subset(refseq,refseq[,3] == "chr2")
chr3_refseq <- subset(refseq,substring(refseq[,3],1,4) == "chr3")
chr4_refseq <- subset(refseq,substring(refseq[,3],1,4) == "chr4")
chr5_refseq <- subset(refseq,substring(refseq[,3],1,4) == "chr5")
chr6_refseq <- subset(refseq,substring(refseq[,3],1,4) == "chr6")
chr7_refseq <- subset(refseq,substring(refseq[,3],1,4) == "chr7")
chr8_refseq <- subset(refseq,substring(refseq[,3],1,4) == "chr8")
chr9_refseq <- subset(refseq,substring(refseq[,3],1,4) == "chr9")
chr10_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr10")
chr11_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr11")
chr12_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr12")
chr13_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr13")
chr14_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr14")
chr15_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr15")
chr16_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr16")
chr17_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr17")
chr18_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr18")
chr19_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr19")
chr20_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr20")
chr21_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr21")
chr22_refseq <- subset(refseq,substring(refseq[,3],1,5) == "chr22")
chrX_refseq <- subset(refseq,substring(refseq[,3],1,4) == "chrX")


totalNOCNVs <- data.frame()

#load individual CNV sample files and combine it to one large file
for (i in 1:10){
  totalNOCNVs <- rbind(totalNOCNVs,getCNVwZeroOL(i))
  
}

write.csv(totalNOCNVs, "CNVs_not_overlapping_DGV-entries.csv")

#fill list with the genes in region
for (i in 1:nrow(totalNOCNVs)){
#   if(totalNOCNVs[i,4] == "X")
#     chr_nr <- 23
#   else chr_nr <- totalNOCNVs[i,4]
#   
#   chr_genes <- switch(chr_nr,chr1_refseq,chr2_refseq,chr3_refseq,chr4_refseq,chr5_refseq,chr6_refseq,chr7_refseq,chr8_refseq,chr9_refseq,chr101_refseq,chr11_refseq,chr12_refseq,chr13_refseq,chr14_refseq,chr15_refseq,chr16_refseq,chr17_refseq,chr18_refseq,chr19_refseq,chr20_refseq,chr21_refseq,chr22_refseq,chrX_refseq)

  genes_per_cnv <- mapToGene(totalNOCNVs[i,4],totalNOCNVs[i,11], totalNOCNVs[i,12])
  
  totalNOCNVs[i,15] <- paste(genes_per_cnv, collapse = ' ')
  totalNOCNVs[i,16] <- finemapCNV(totalNOCNVs[i,4],totalNOCNVs[i,11], totalNOCNVs[i,12],genes_per_cnv)
}

#write results to excel sheet



noDGV <- totalNOCNVs

group_list <- data.frame()

while(nrow(noDGV) > 0){
  curr_Chr <- noDGV[1,4]
  curr_start <- noDGV[1,11] 
  curr_stop <- noDGV[1,12]
  
  i <- 1
  
  temp_index <- c()
  for(i in 1:nrow(noDGV)){
    if(curr_Chr == noDGV[i,4]){
      if((curr_start <= noDGV[i,11] & curr_stop <= noDGV[i,12] & curr_stop >= noDGV[i,11]) | (curr_start >= noDGV[i,11] & curr_stop >= noDGV[i,12] & curr_start <= noDGV[i,12]) | (curr_start >= noDGV[i,11] & curr_stop <= noDGV[i,12]) | (curr_start <= noDGV[i,11] & curr_stop >= noDGV[i,12])){
        group_list <- rbind(group_list,noDGV[i,])
        #cat(noDGV[i,1],noDGV[i,4],noDGV[i,11],"\n")
        
        temp_index <- c(temp_index,i)
      }
      
    }
  }
  
  noDGV <- noDGV[-temp_index,]
  #cat(noDGV[1,1],noDGV[1,4],"\n")
  group_list <- rbind(group_list,rep(" ",ncol(group_list)))
}



