x = 'AGCAGACACGTGAT'
y = 'ATCACCGGTAT'

backtrack_lcs <- function(b,x,m,n)
{
  if(m==0 | n==0) return(NULL)
  if(b[m+1,n+1] == '\\')
  {
    return(c(x[m],backtrack_lcs(b,x,m-1,n-1)))
  }
  else if(b[m+1,n+1] == '|')
  {
    backtrack_lcs(b,x,m-1,n)
  }
  else
  {
    backtrack_lcs(b,x,m,n-1)
  }
}

#unlist(strsplit(x,''))

find_lcs <- function(x,y)
{
  options(expressions=10000)
  
  # unlist to each characters
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  # get the length of each seq
  m <- length(x)
  n <- length(y)
  
  # set up the backtrack key dict
  backtrack_key <- c('|','--','\\')
  
  # set up the score and backtrack matrix
  s <- matrix(0,m+1,n+1,dimnames=list(c('',x),c('',y))) 
  b <- matrix('',m+1,n+1,dimnames=list(c('',x),c('',y)))
  
  for(i in seq(2,m+1))
  {
    for (j in seq(2,n+1))
    {
      s_up <- s[i-1,j] # upper col
      s_left <- s[i,j-1] # left row
      s_diag <- s[i-1,j-1] + 1 # diag with +1 for a match
      
      if(x[i-1] == y[j-1])
      {
        scores <- c(s_up,s_left,s_diag)
      }
      else
      {
        scores <- c(s_up,s_left)
      }
      
      backtrack_update <- which.max(scores)
      score_update <- max(scores)
      
      s[i,j] <- score_update
      b[i,j] <- backtrack_key[backtrack_update]
                  
    }
  } 
  lcs <- backtrack_lcs(b,x,m,n)
  
  return(list(lcs=paste0(rev(lcs),collapse = ''),
              length=s[length(s)],
              scores=s,
              backtrack=b))
}

find_lcs('AGCAGACACGTGAT','ATCACCGGTAT')