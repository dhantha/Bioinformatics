pam250 <- as.matrix(read.table('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/0b783786e3797d1bb55172e750776e26021224b0/PAM250.dat'))
xy <- readr::read_lines('https://gist.githubusercontent.com/sw1/8870a124624f31585d5c15be72fcfc21/raw/6617b2ceae8765ff7b91f4a600fc5460b335279d/local_sequences')

backtrack_global <- function(x,y,b,penalty,score_matrix,s,max_score_index)
{
	m <- max_score_index[1] #length(x)
	n <- max_score_index[2] #length(y)
	
	score <- 0
	align_x <- NULL
	align_y <- NULL
	
	while(m > 0 | n > 0)
	{
		if (s[m,n] == 0){
			# set the m,n to 0 to end the backtrack   
			n <- 0
			m <- 0
		}
		else if(b[m+1,n+1] == '|'){
			align_x <- c(align_x,x[m])
			align_y <- c(align_y,'-')
			score <- score + penalty
			m <- m-1
		}else if(b[m+1,n+1] == '--'){
			align_x <- c(align_x,'-')
			align_y <- c(align_y,y[n])
			score <- score + penalty
			n <- n-1
		}
		else{
			align_x <- c(align_x,x[m])
			align_y <- c(align_y,y[n])
			score <- score + score_matrix[x[m],y[n]]
			n <- n-1
			m <- m-1
		}
	}
	
	alignment <- c(paste0(rev(align_x),collapse=''),paste0(rev(align_y),collapse=''))
  
	return(list(score=score,alignment=alignment))
}

library(compiler)
# cmpfun
local_alignment <- cmpfun(function(x,y,penalty,score_matrix)
{
	# x,y is the input seq
	# penalty is the penalty
	# score_matrix_penalty is the scoring lookup table
	
	x <- unlist(strsplit(x,'')) # input seq
	y <- unlist(strsplit(y,''))
	
	m <- length(x)
	n <- length(y)
	
	backtrack_key <- c('|','--','\\')
	
	#print(m,n)
	# empty matrix s for scores 
	s <- matrix(0,length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
	s[1,] <- rep(0,ncol(s))
    s[,1] <- rep(0,nrow(s))
	
	# empty matrix for backtracking path 
	b <- matrix('',length(x)+1,length(y)+1,dimnames=list(c('',x),c('',y)))
	b[1,] <- '--'
	b[,1] <- '|'
	b[1,1] <- '\\'
	
	
	for (i in seq(2,m+1)){
		for (j in seq(2,n+1)){
      
		s_up <- s[i-1,j] + penalty
		s_left <- s[i,j-1] + penalty
		s_diag <- s[i-1,j-1] + score_matrix[x[i-1],y[j-1]]
	    
	    scores <- c(s_up,s_left,s_diag,0)  # and 0
	    
		# get the position of the max score
	    backtrack_update <- which.max(scores)
		# update the scores
	    score_matrix_update <- max(scores)
		
		if (score_matrix_update < 0){
			score_matrix_update <- 0
		}	
		
		s[i,j] <- score_matrix_update
		b[i,j] <- backtrack_key[backtrack_update]
	    
	  }
	}
	max_score_index <- which(s == max(s),arr.ind = TRUE)
	
	return(backtrack_global(x,y,b,penalty,score_matrix,s,max_score_index))
})


x <- xy[1]
y <- xy[2]
score_matrix <- pam250
penalty <- -5
local_alignment(x,y,score_matrix,penalty)