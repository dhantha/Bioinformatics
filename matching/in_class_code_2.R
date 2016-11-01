source('https://dl.dropboxusercontent.com/u/907375/Bioinformatics_class/R_inclass_functions.R')

load_library(Biostrings)
load_library(ggplot2)
load_library(readr)
load_library(stringr)
load_library(reutils)
load_library(gridExtra)

cleanup_feat_table <- function(x){
   xx <- x$GBFeature_quals 
   xxx <- matrix(unlist(xx),ncol=2,byrow=TRUE)
   xxxx <- xxx[,2]
   names(xxxx) <- xxx[,1]
   xxxx
}

############################
############################
### FASTA AND BIOSTRINGS ###
############################
############################

fasta <- readDNAStringSet('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/f1fb586160d12c34f29532c731066fd8912a0e0c/example.fasta',format='fasta')
fasta

View(fasta)

length(fasta)
width(fasta)
names(fasta)

subseq(fasta,start=1,end=5)
subseq(fasta,start=5,end=25)
subseq(fasta,start=5,width=10)
subseq(fasta,start=5,end=-5)

alphabetFrequency(fasta)
alphabetFrequency(fasta,baseOnly=TRUE)
alphabetFrequency(fasta,baseOnly=TRUE,as.prob=TRUE)
round(alphabetFrequency(fasta,baseOnly=TRUE,as.prob=TRUE),3)

letterFrequency(fasta,'G',as.prob=TRUE)
letterFrequency(fasta,'GC',as.prob=TRUE)
letterFrequency(fasta,'ACGT',as.prob=TRUE)
letterFrequency(fasta,'ACGT',OR=0,as.prob=TRUE)

consensusMatrix(fasta,baseOnly=TRUE)
consensusMatrix(fasta,baseOnly=TRUE,width=100)
dinucleotideFrequency(fasta)
trinucleotideFrequency(fasta)
oligonucleotideFrequency(fasta,width=5)

RNAStringSet(fasta)
complement(fasta)
reverseComplement(fasta)
translate(reverseComplement(fasta))

############################
############################
########## FASTQ ###########
############################
############################

fastq <- readDNAStringSet('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/f1fb586160d12c34f29532c731066fd8912a0e0c/example.fastq',format='fastq')
fastq

############################
############################
## CREATING SEQUENCE SETS ##
############################
############################

problem_createsequencesets()
s1
s2
s3

S <- c(s1,s2,s3)
SS <- DNAStringSet(S)
SS

names(SS) <- c('sequence_1','sequence_2','sequence_3')


############################
############################
########### PASTE ##########
############################
############################


DOG <- c('dog1','dog2','dog3')
CAT <- c('cat1','cat2','cat3')

paste(DOG,CAT)
paste(DOG,CAT,sep='-')
paste(DOG,CAT,sep='_')
paste(DOG,CAT,sep='')
paste('dog','cat',1:3,sep='')
paste('dog','cat',1:3,sep='_')
paste('dog_','cat',1:3,sep='')

seq_names <- paste('sequence_',1:length(SS),' | User_12 | ',date(), sep='')
seq_names

names(SS) <- seq_names
SS
View(SS)

output_name <- 'seq_set_out.fasta'
writeXStringSet(SS,file=output_name,format="fasta")


############################
############################
######## METADATA ##########
############################
############################


FASTA <- readDNAStringSet('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/10bc2f50d1c739827ea2ba4edb146b36a6a4c14a/problems_metadata.fasta',format='fasta')
problem_metadata(FASTA)
FASTA
head(META)


which(META$Center == 'Philadelphia')
FASTA['Sequence_6333']

header_names <- META$ID[c(12,15,78)]
FASTA[header_names]

SS <- FASTA[META$Center == 'Houston' & META$Genus == 'Escherichia']
SS

output_name <- 'seq_subset_out.fasta'
writeXStringSet(SS,file=output_name,format="fasta")


############################
############################
####### GC FUNCTIONS #######
############################
############################


gc_calc <- function(x) (x['g']+x['c'])/sum(x)

gc <- function(s,pos){
  
  s <- stringr::str_to_lower(s)
  s <- unlist(strsplit(s,''))
  
  if (!missing(pos)) s <- s[seq(pos,length(s),3)]
  counts <- table(s)
  
  gc_calc(counts)
  
}


gc_skew_calc <- function(x) {counts <- table(x); (counts['g']-counts['c'])/(counts['g']+counts['c'])}

gc_skew <- function(s,win){
  
  s <- stringr::str_to_lower(s)
  s <- unlist(strsplit(s,''))
  
  if (missing(win)) {
    gc <- gc_skew_calc(s)
  }else{
    start <- seq(1,length(s),win)
    gc <- NULL
    for (i in start){
      gc <- c(gc, gc_skew_calc(s[(i):(i+win-1)]))
    }
  }
  
  gc
  
}


generate_random_dna_gc_s(len=1000,seed=1)

gc(s)
gc(s,1)
gc(s,2)
gc(s,3)

generate_random_dna_skew_s(len=1000,w=1,seed=5)

gc_skew(s)
gc_skew(s,100)
plot_skew(gc_skew(s,25))


############################
############################
######## NCBI SEARCH #######
############################
############################

ids1 <- esearch("CFTR AND human[Organism] AND complete",db='nucleotide',retmax=15,sort='relevance')
ids2 <- esearch("PKD1 AND human[Organism] AND complete",db='nucleotide',retmax=15,sort='relevance')
ids3 <- esearch("DMPK AND human[Organism] AND complete",db='nucleotide',retmax=15,sort='relevance')

ids_df <- reutils::content(esummary(ids1),'parsed')
View(ids_df)

efetch(ids1[1], rettype = "fasta", retmode = "text")
efetch(ids2[4], rettype = "fasta", retmode = "text")
efetch(ids3[5], rettype = "fasta", retmode = "text")

ids <- c(ids1[1],ids2[4],ids3[5])
ids

FASTA <- efetch(ids,db='nucleotide', rettype = "fasta", retmode = "xml")
LENS <- FASTA$xmlValue('//TSeq_length')
LENS
SEQS <- FASTA$xmlValue('//TSeq_sequence')
SEQS

tmp <- tempfile()
FASTA <- efetch(ids,db='nucleotide', rettype = "fasta", retmode = "text", outfile=tmp)
FASTA <- readDNAStringSet(tmp)

############################
############################
######### CDS SEARCH #######
############################
############################

ID <- esearch("Galdieria sulphuraria[Organism] AND whole genome",db='nucleotide',retmax=5,sort='relevance')
rec <- efetch(ID[1],db='nucleotide', rettype = "gb", retmode = "xml")
rec

prec <- reutils::content(rec,as='text')
prec <- xmlParse(prec)
prec <- xmlToList(prec)
prec

features <- prec$GBSeq$`GBSeq_feature-table`
cds_idx <- which(sapply(features,function(x) x[[1]]) == 'CDS')
features <- features[cds_idx]

features <- lapply(features,cleanup_feat_table)
na.omit(sapply(features,function(x) ifelse(grepl('ATPase',x['product']),x['protein_id'],NA)))

############################
############################
######## GC RESULTS ########
############################
############################


letterFrequency(FASTA,'GC',as.prob=TRUE)

skew <- gc_skew(FASTA[[1]],500)
p1 <- plot_skew(skew)

skew <- gc_skew(FASTA[[2]],500)
p2 <- plot_skew(skew)

skew <- gc_skew(FASTA[[3]],500)
p3 <- plot_skew(skew)

gridExtra::grid.arrange(p1,p2,p3,nrow=1)


skew <- lapply(seq_along(FASTA), function(i,w) {
  numer <- letterFrequencyInSlidingView(FASTA[[i]],'G',view.width=w) -
    letterFrequencyInSlidingView(FASTA[[i]],'C',view.width=w)
  denom <- letterFrequencyInSlidingView(FASTA[[i]],'GC',view.width=w)
  numer/denom
},w=500)

plot_skew(skew[[2]])


############################
############################
###### WHOLE GENOMES #######
############################
############################

load_library(BSgenome)

load_library(BSgenome.Athaliana.TAIR.04232008)
load_library(BSgenome.Osativa.MSU.MSU7)

Athaliana
Osativa

Athaliana$chr1


params <- new('BSParams',
              X=Athaliana,
              FUN = function(x) letterFrequency(x,'GC',as.prob=TRUE),
              exclude=c('M','C'))
unlist(bsapply(params))


params <- new('BSParams',
              X=Osativa,
              FUN = function(x) letterFrequency(x,'GC',as.prob=TRUE),
              exclude=c('ChrC','M','Un','Sy'))
unlist(bsapply(params))

