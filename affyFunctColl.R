require(xlsx)


# conversion of Affy Genotyping Console text export SNP6 bulk genotyping
# data to  vcf format. Required, as AffyCalls are exported unordered, without genomic position 
# and deviating call symbols
affyToVCF <- function(x , ref.calls){
  affy.calls <- x
  total.rows <- nrow(affy.calls)
  total.columns <- (ncol(affy.calls) - 2)/2 + 9
  converted.vcf <- data.frame(matrix(c(rep(0,total.rows * 9)), nrow=total.rows), 
                              stringsAsFactors=F)
  colnames(converted.vcf) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  
  #align snps in affy calls to snps in ref.calls
  #this is ABSOLUTELY CRUCIAL, otherwise, association of ref alleles
  #to SNPs won't be correct.
  affy.calls <- affy.calls[match(ref.calls[,1], affy.calls[,ncol(affy.calls)]),]
  
  #fill converted.vcf with rs#, ref and alt alleles
  converted.vcf[, 3:5] <- ref.calls[, 1:3]
  print(sapply(ref.calls[1,],class))
  #sapply(converted.vcf[,4], print)
  
  for(i in seq(2, ncol(affy.calls) - 2, 2)) {
    #if ref allele is different in affy calls, switch to desired ref allele
    #conversion is done column-wise
    
    is.ref <- (substr(affy.calls[, i + 1], 1, 1) == converted.vcf[, 4] & 
      affy.calls[, i] == "AA") | 
      (substr(affy.calls[, i + 1], 1, 1) == converted.vcf[, 5] & 
         affy.calls[, i] == "BB")
    
    fill.AA <- affy.calls[, i] == "BB" & substr(affy.calls[, i + 1], 1, 1) == converted.vcf[, 4] & !is.ref
    fill.BB <- affy.calls[, i] == "AA" & substr(affy.calls[, i + 1], 1, 1) == converted.vcf[, 5] & !is.ref
    
    affy.calls[fill.AA, i] <- "AA"
    affy.calls[fill.BB, i] <- "BB"
    
    #convert affy calls to vcf format
    affy.calls[,i] <- sapply(affy.calls[,i], function(y) {
      switch(y, "NoCall" = "./.",
             "AA" = "0/0",
             "AB" = "0/1",
             "BB" = "1/1")})
    
    converted.vcf <- cbind(converted.vcf, affy.calls[,i])
    colnames(converted.vcf)[ncol(converted.vcf)] <- colnames(affy.calls)[i]
    converted.vcf[, ncol(converted.vcf)] <- as.character(converted.vcf[, ncol(converted.vcf)])
  }
  
  converted.vcf[,3] <- affy.calls[,ncol(affy.calls)]  
  return(converted.vcf)
}

#takes a vcf file and hapmap population and returns a boolean list were True indicates
#the column/sample belonging to the population. Populations are taken from column name
getPopulationIndex <- function(vcf, pop, pos = F) {
  vcf.headers <- colnames(vcf)
  pop.indices <- c()
  if(pos)
    pop.indices <- grep(paste(pop, "$", sep="", collapse=""), vcf.headers)
  else
    pop.indices <- grepl(paste(pop, "$", sep="", collapse=""), vcf.headers)
  
  #col 1 to 9 are vcf standard fields and don't contain genotyping info, thus set to F
  pop.indices[1:9] <- F
  
  return(pop.indices)
}



# plots a list of snps on the x-axis and the corresponding sample on the y-axis
# takes a vcf format table as input, therefore, genotyping data is from column 10 on
plotSNPs <- function(snp.list, xdim, ydim) {
  
  snp.count <- nrow(snp.list)
  
  # initialize the plot in square mode and write to svg file
  svg(filename = "test.svg", width = xdim, height = ydim)
  par(pty="s", bty="l")
  plot(c(0, snp.count), c(0, ncol(snp.list)), asp=1, col="white")
  axis(side = 1, labels = snp.list[,3], at = (c(1:snp.count) - 0.5), pos = 0, las = 3,
       cex.axis = 0.4)
  axis(side = 2, labels = colnames(snp.list)[10:ncol(snp.list)], at = (c(1:(ncol(snp.list) - 9)) - 0.5), 
       pos = 0, las=1, cex.axis = 0.4)
  
  axis(side = 3, labels = snp.list[, 2], at = (c(1:snp.count) - 0.5), 
       pos = ncol(snp.list) - 9, las = 3, cex.axis = 0.4, col = "black")
  
  
  for(i in 10:ncol(snp.list)){
    
  
    for(y in 1:snp.count) {
      color <- switch(snp.list[y, i], 
                      "./." = "white",
                      "0/0" = "red",
                      "0/1" = "#6495ED",
                      "1/1" = "green")     
      
      rect(xleft = (y - 1), ybottom = (i - 10), xright = y, ytop = (i - 10 + 1), col = color)
    }
  }
  
  dev.off()
}


