pam250 <- as.matrix(read.table('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/0b783786e3797d1bb55172e750776e26021224b0/PAM250.dat'))
xy <- readr::read_lines('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/6617b2ceae8765ff7b91f4a600fc5460b335279d/local_sequences')

backtrack_global <- function(x,y,b,score_matrix,penalty,s,max_index){
  
  m <- max_index[1]
  n <- max_index[2]
  
  score <- 0
  align_x <- NULL
  align_y <- NULL
  
  while (m > 0 | n > 0){
    if (s[m,n] == 0){
      alignment <- c(paste0(rev(align_x),collapse=''),paste0(rev(align_y),collapse=''))
      
      return(list(score=score,alignment=alignment))
    }else{
      if (b[m+1,n+1] == '|'){
        align_x <- c(align_x,x[m])
        align_y <- c(align_y,'-')
        score <- score + penalty
        m <- m-1
      }else if(b[m+1,n+1] == '--'){
        align_x <- c(align_x,'-')
        align_y <- c(align_y,y[n])
        score <- score + penalty
        n <- n-1
      }else{
        align_x <- c(align_x,x[m]) 
        align_y <- c(align_y,y[n])
        score <- score + score_matrix[x[m],y[n]]
        n <- n-1
        m <- m-1
      }
      }
  }
  
  alignment <- c(paste0(rev(align_x),collapse=''),paste0(rev(align_y),collapse=''))
  
  return(list(score=score,alignment=alignment))
  
}

library(compiler)
local_alignment <- cmpfun(function(x,y,score_matrix,penalty){
  
  x <- unlist(strsplit(x,''))
  y <- unlist(strsplit(y,''))
  
  m <- length(x)
  n <- length(y)
  
  backtrack_key <- c('|','--','\\')
  
  s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  s[1,] <- rep(0,ncol(s))
  s[,1] <- rep(0,nrow(s))
  
  b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
  b[1,] <- '--'
  b[,1] <- '|'
  b[1,1] <- '\\'
  
  for (i in seq(2,m+1)){
    for (j in seq(2,n+1)){
      
      s_up <- s[i-1,j] + penalty
      s_left <- s[i,j-1] + penalty
      s_diag <- s[i-1,j-1] + score_matrix[x[i-1],y[j-1]]
      
      scores <- c(s_up,s_left,s_diag) 
      
      backtrack_update <- which.max(scores)
      score_matrix_update <- max(scores)
      if (score_matrix_update < 0){ 
        score_matrix_update <- 0
      }
      
      s[i,j] <- score_matrix_update
      b[i,j] <- backtrack_key[backtrack_update]
      
    }
  }
   max_index <- which(s == max(s),arr.ind = TRUE)
   
  return(backtrack_global(x,y,b,score_matrix,penalty,s,max_index))
  
})

x <- xy[1]
y <- xy[2]
score_matrix <- pam250
penalty <- -5
local_alignment(x,y,score_matrix,penalty)
