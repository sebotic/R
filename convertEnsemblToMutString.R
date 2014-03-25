#set of functions to convert the genomic positions of Ensembl variation/mutation data, to data interpretable by a clinician, 
#namely: mutation position in a transcript and affected amino acid
#deals with indels and SNVs 

#needs code cleanup

#copyright Sebastian Burgstaller, 2014

require(xlsx)

gen.subst <- function(x, coding.start.pos, plus.strand = T, pos.provided = "") {
  position <-  as.integer(unlist(strsplit(x[1, 2], ":"))[3])
  if(plus.strand){
    coding.pos <- position - coding.start.pos + 1
  } else {
    coding.pos <- coding.start.pos - position + 1
  }
  
  if(pos.provided != ""){
    coding.pos <- pos.provided
  }
    
  ref.allele <- unlist(strsplit(x[1, 3], "/"))[1]
  alt.allele <- unlist(strsplit(x[1, 3], "/"))[2]
  
  
  dna.string <- paste0("c.", coding.pos, ref.allele, ">", alt.allele, sep="")
  prot.string <- ""
  
  if(!(x[1,10] == "-" |  x[1,10] == "" | x[1,10] == " ")){
    if(grepl("/", x[1,10])){
      ref.as <- unlist(strsplit(x[1, 10], "/"))[1]
      alt.as <- unlist(strsplit(x[1, 10], "/"))[2]
      if(alt.as == "*") 
        alt.as <- "X"
    } else {
      ref.as <- x[1, 10]
      alt.as <- x[1, 10]
    }
    
    
    prot.string <- paste0("p.", ref.as, x[1, 11], alt.as, sep="")
  }
  
  return(c(dna.string, prot.string))  
}

gen.insertion <- function(x, coding.start.pos, plus.strand = T, pos.provided = "") {
  ins.start <-  as.integer(unlist(strsplit(x[1, 2], " "))[3])
  ins.end <- as.integer(unlist(strsplit(x[1, 2], " "))[5])
  ins.start <- abs(ins.start - coding.start.pos) + 1
  ins.end <- abs(ins.end - coding.start.pos) + 1
  ins.size <- abs(ins.start - ins.end) + 1
  
  if(plus.strand){
    if(pos.provided != 0){
      dna.string <- paste0("c.", pos.provided, "..", pos.provided + ins.size, "_ins")
    } else
      dna.string <- paste0("c.", ins.start, "..", ins.end, "_ins")
  } else {
    if(pos.provided != 0){
      dna.string <- paste0("c.", pos.provided, "..", pos.provided + ins.size, "_ins")
    }else
      dna.string <- paste0("c.", ins.end, "..", ins.start, "_ins")
  }
  
  dna.string <- paste0(dna.string, unlist(strsplit(x[1, 3], "/"))[2])
  
  return(c(dna.string, ""))
}

gen.deletion <- function(x, coding.start.pos, plus.strand = T, pos.provided = "") {
  tmp.string <- unlist(strsplit(x[1, 2], ":"))[3]
  del.start <-  as.integer(unlist(strsplit(tmp.string, "-"))[1])
  del.end <- as.integer(unlist(strsplit(tmp.string, "-"))[2])
  del.start <- abs(del.start - coding.start.pos) + 1
  del.end <- abs(del.end - coding.start.pos) + 1
  del.size <- abs(del.start - del.end) + 1
  
  if(plus.strand){
    if(pos.provided != ""){
      dna.string <- paste0("c.", pos.provided, "..", pos.provided + del.size, "_del")
    }else 
      dna.string <- paste0("c.", del.start, "..", del.end, "_del")
  } else {
    if(pos.provided != ""){
      dna.string <- paste0("c.", pos.provided, "..", pos.provided + del.size, "_del")
    }else
      del.string <- paste0("c.", del.end, "..", del.start, "_del")
  }
  
  dna.string <- paste0(dna.string, unlist(strsplit(x[1, 3], "/"))[1])
  
  return(c(dna.string, ""))  
  
}

#construct a list wich denotes the last coding base of every exon
getBorderBase <- function(exon.starts, exon.ends, translation.start, translation.end){
  border.base <- c()
  total.bases <- 0
  for(i in 1:(length(exon.starts)-1)) {
    if(translation.start > 0) {
      border.base <- c(abs(exon.starts[i] - exon.ends[i]) + 1 - abs(translation.start - exon.starts[i])) 
      translation.start <- 0
      total.bases <- border.base
    }else {
      tmp <- abs(exon.starts[i] - exon.ends[i]) + total.bases
      total.bases <- tmp
      border.base <- c(border.base, total.bases)
    }
    
  }
  return(border.base)
}

genRelCodingPos <- function(position, exon.starts, exon.ends, translation.start, translation.end, plus.strand = T){
  is.intron <- F
  is.exon <- F
  count <- 0
  rel.pos <- 0
  border.bases <- getBorderBase(exon.starts, exon.ends, translation.start)
  
  if(plus.strand){
    for(i in 1:length(exon.starts)){
      if(translation.start > position) {
        rel.pos <- translation.start - position + 1
        break
      } else if(exon.starts[i] <= position & exon.ends[i] >= position) {
        is.exon <- T
        count <- i
        if(count == 1){
          rel.pos <- position - tranlation.start + 1
        } else {
          rel.pos <- position - exon.starts[i] + 1 + border.bases[i]
        }
        break
      } else if(position < exon.ends[i] & position < exon.starts[i + 1]){
        is.intron <- T
        intron.length <- abs(exon.ends[i] - exon.starts[i + 1]) + 1
        if((exon.ends[i] + intron.length/2) < position){
          rel.pos <- paste0(border.bases[i], "+", (position - exon.ends[i] + 1))
        } else {
          rel.pos <- paste0(border.bases[i] + 1, "-", (position - exon.start[i + 1] - 1))
        }
        break
      }
    }
    
  }else{
    for(i in 1:length(exon.starts)){
      if(translation.start < position) {
        rel.pos <- abs(translation.start - position) + 1
        break
      } else if(exon.starts[i] >= position & exon.ends[i] <= position) {
        is.exon <- T
        count <- i
        print(i)
        if(count == 1){
          rel.pos <- abs(position - tranlation.start) + 1
        } else {
          rel.pos <- abs(position - exon.starts[i]) + 1 + border.bases[i - 1]
          
        }
        break
      } else if(position > exon.ends[i] & position > exon.starts[i + 1]){
        is.intron <- T
        intron.length <- abs(exon.ends[i] - exon.starts[i + 1]) + 1
        if((exon.ends[i] + intron.length/2) > position){
          rel.pos <- paste0(border.bases[i - 1]   + 1, "-", (abs(position - exon.starts[i]) + 1))
        } else {
          rel.pos <- paste0(border.bases[i - 1], "+", (abs(position - exon.ends[i - 1]) + 1))
        }
        break
      }
    }
  }
  #rel.pos <- paste0("c.", rel.pos)
  return(rel.pos)
}


mc1r <- read.csv("ALL-VariationTable-Homo_sapiens-Transcript-Variation_Transcript-75-ENST00000555427.csv", 
                 stringsAsFactors=F)

plus.one <- 89985667

#really required columns c(1:3, 9:11)
mc1r <- mc1r[380:605,]
mc1r <- subset(mc1r, mc1r[,3] != "HGMD_MUTATION")
mc1r <- mc1r[!grepl("somatic", mc1r[, 5]), ]

mc1r.mut.cleaned <- data.frame()

for(i in 1:nrow(mc1r)){
  type <- unlist(strsplit(mc1r[i, 3], "/"))
  if(type[1] == "-")
    mut <- gen.insertion(mc1r[i, ], plus.one)
  else if(grepl("^[ATCG]", type[1]) & type[2] == "-")
    mut <- gen.deletion(mc1r[i, ], plus.one)
  else if(grepl("^[ATCG]", type[1]) & grepl("^[ATCG]", type[2]))
    mut <- gen.subst(mc1r[i, ], plus.one)
  
  mc1r.mut.cleaned[i, 1] <- mut[1]
  mc1r.mut.cleaned[i, 2] <- mut[2]
}

write.xlsx(mc1r.mut.cleaned, "varianten.xlsx", sheetName="MC1R", col.names = T, row.names = F, append = T)




