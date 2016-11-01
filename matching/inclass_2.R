source('https://dl.dropboxusercontent.com/u/907375/Bioinformatics_class/R_inclass_functions.R')

load_library(Biostrings)
load_library(ggplot2)
load_library(readr)
load_library(stringr)
load_library(reutils)
load_library(XML)

# fasta Vs fastQ

fasta <- readDNAStringSet('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/f1fb586160d12c34f29532c731066fd8912a0e0c/example.fasta',format='fasta')
fasta

problem_createsequencesets()
S <- c(s1,s2,s3)
SS <- DNAStringSet(S)

names(SS) <- c('sequence_1','sequence_2','sequence_3')

seq_names <- paste('sequence_',1:length(SS),' | User_12 | ',date(), sep='')
seq_names

names(SS) <- seq_names
output_name <- 'seq_set_out.fasta'
writeXStringSet(SS,file=output_name,format="fasta")


FASTA <- readDNAStringSet('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/10bc2f50d1c739827ea2ba4edb146b36a6a4c14a/problems_metadata.fasta',format='fasta')
problem_metadata(FASTA)

which(META$Center == 'Philadelphia')
FASTA['Sequence_6333']

# calculate GC function
gc_calc <- function(x)
{
  (x['g']+['c']/sum(x))
}

gc <- function(s,pos)
{
  s <- stringr::str_to_lower(s)
  s <- unlist(strsplit(s,''))
  
  if(!missing(pos)) s <- s[seq(pos,length(s),3)]
  
  counts <table(s)
  
  gc_calc(counts)
}

# calculate gc_skew

gc_skew_calc <- function(x)
{
  counts <- table(x)
  (counts['g']-counts['c'])/(counts['g']+counts['g']+counts['c'])
}

gc_skew <- function(s,win)
{
  s <- stringr::str_to_lower(s)
  s <- unlist(strsplit(s,''))
  
  if(missing(win))
  {
    gc <- gc_skew_calc(s)
  }
  else
  {
    start <- seq(1,length(s),win)
    gc <- NULL
    for(i in start)
    {
      gc <- c(gc,gc_skew_calc(s[(i):(i+win-1)]))
    }
  }
  gc
}

generate_random_dna_skew_s(len=1000,w=1,seed=5)
gc_skew(s)

ids1 <- esearch("CFTR AND human[Organism] AND complete",db='nucleotide',retmax=15) #,sort='relevance')
ids2 <- esearch("PKD1 AND human[Organism] AND complete",db='nucleotide',retmax=15) #,sort='relevance')
ids3 <- esearch("DMPK AND human[Organism] AND complete",db='nucleotide',retmax=15) #,sort='relevance')

ids_df <- reutils::content(esummary(ids1),'parsed')

efetch(ids1[1], rettype = "fasta", retmode = "text")
efetch(ids2[4], rettype = "fasta", retmode = "text")
efetch(ids3[5], rettype = "fasta", retmode = "text")

ids <- c(ids1[1],ids2[4],ids3[5])

FASTA <- efetch(ids,db='nucleotide', rettype = "fasta", retmode = "xml")
SEQS <- FASTA$xmlValue('//TSeq_sequence')

tmp <- tempfile()
FASTA <- efetch(ids,db='nucleotide', rettype = "fasta", retmode = "text") #, outfile=tmp)
FASTA <- readDNAStringSet(tmp)

letterFrequency(FASTA,'GC',as.prob=TRUE)
skew <- gc_skew(FASTA[[2]],500)


ID <- esearch("Galdieria sulphuraria[Organism] AND whole genome",db='nucleotide',retmax=5,sort='relevance')
rec <- efetch(ID[1],db='nucleotide', rettype = "gb", retmode = "xml")

prec <- reutils::content(rec,as='text')
prec <- xmlParse(prec)
prec <- xmlToList(prec)
features <- prec$GBSeq$`GBSeq_feature-table`
cds_idx <- which(sapply(features,function(x) x[[1]]) == 'CDS')
features <- features[cds_idx]
features <- lapply(features,cleanup_feat_table)
na.omit(sapply(features,function(x) ifelse(grepl('ATPase',x['product']),x['protein_id'],NA)))


load_library(BSgenome)
available.genomes()

load_library(BSgenome.Athaliana.TAIR.04232008)
load_library(BSgenome.Osativa.MSU.MSU7)

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
