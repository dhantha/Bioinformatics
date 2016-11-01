source('https://goo.gl/rVBfnT')

generate_random_dna_sequence_s(len=1000,seed=5)

ss <- unlist(strsplit(s,''))
table(ss) # split into indiv characters

table_ss <- table(ss)
total_ss <- sum(table_ss)

(table_ss['G'] + table_ss['C'])/total_ss

target_sequence <- c('CAT')
for (index in seq(from=1,to=length(ss),by=3))
{
  subsequence <- ss[index:(index+2)]
  subsequence <- paste(subsequence,collapse = '')
  
  if(subsequence == target_sequence) cat('Match at pos', index,'\n')
  
}


target_sequence <- c('CAT')
for (index in seq(from=1,to=length(ss)-1,by=3))
{
  subsequence <- ss[index:(index+2)]
  subsequence <- paste(subsequence,collapse='')
  if (subsequence == target_sequence) cat('Match at pos.',index,'\n')
}

target_sequence <- c('CAT')
total <- 0
for (index in seq(from=1,to=length(ss)-1,by=3))
{
  subsequence <- ss[index:(index+2)]
  subsequence <- paste(subsequence,collapse='')
  if (subsequence == target_sequence) {
    cat('Match at pos.',index,'\n')
    total <- total + 1
  }
}

