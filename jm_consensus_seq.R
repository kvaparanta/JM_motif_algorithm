jm_consensus_seq <-function(path_to_fasta_file, path_to_provided_PAM250_file){

  #download and install package 'seqinr'and 'nnet' 
  #as an input give the path to the .fasta file, where your sequences are in fasta format 
  #and the path to the provided PAM250.txt file
  #returns list of potential sites in the sequences and their position in the sequence
  PAM250<- read.table(path_to_provided_PAM250_file, sep=';', header=TRUE, row.names = 1)
  
  seqs <- read.fasta(file = path_to_fasta_file,      #reads the sequences from the fasta file
                     seqtype = "AA", as.string = TRUE)
  seqs2<- getSequence(seqs, as.string = TRUE) #sequences as strings
  s<-unlist(seqs2) #unlist to string array
  names<- getName(seqs, as.string = TRUE) # names of sequences
  
  #initialization
  n_seqs<-length(s) 
  max_length<- max(nchar(s))-5 
  scores <- array(rep(0, max_length*6*n_seqs), dim=c(max_length, 6, n_seqs))
  HNQRK <- c(2,3,6,9,12)
  EDQ <- c(4,7,6)
  HDE <- c(9,7,4)
  G <- c(5,8)
  
  
  for (y in 1:n_seqs){
    seq<-s[y]
    seq2<-strsplit(seq,NULL)
    seq3<-unlist(seq2)
    l <- length(seq3)-5
    for (i in 1:l){
      score_l <- 0
      score_s <- 0
      score_lr <- 0
      score_sr <- 0
      score_JMbL <- 0
      score_JMbs <- 0
      for (g in c(1,4,5,6)){
        if (g==1){
          k<-which(seq3[i+(g-1)]== colnames(PAM250))
          score_l <- score_l + PAM250[G,k]
          score_s <- score_s + PAM250[G,k]
          score_lr <- score_lr + PAM250[G,k]
          score_sr <- score_sr + PAM250[G,k]
          score_JMbL <- score_JMbL + PAM250[G,k]
          score_JMbs <- score_JMbs + PAM250[G,k]
          
        }
        if (g==4){
          k<-which(seq3[i+(g-1)]== colnames(PAM250))
          score_s <- score_s + mean(PAM250[HNQRK,k])
          score_sr <- score_sr + mean(PAM250[HDE,k])
          score_JMbs <- score_JMbs + mean(PAM250[EDQ,k])
        }
        if (g==5){
          k<-which(seq3[i+(g-1)]== colnames(PAM250))
          k_bef <-which(seq3[i+(g-2)]== colnames(PAM250))
          score_s <- score_s + mean(PAM250[HDE,k])
          score_l <- score_l + mean(PAM250[HNQRK,k])
          score_lr <- score_lr + mean(PAM250[HDE,k])
          score_sr <- score_sr + mean(PAM250[HNQRK,k])
          score_JMbL <- score_JMbL + mean(PAM250[EDQ,k])
          if (colnames(PAM250)[k] %in%  c('E','D')){
            score_JMbs <- score_JMbs + mean(PAM250[EDQ,k])
          }
          if (colnames(PAM250)[k] %in%  c('H','N','R','Q','K')){
            score_JMbs <- score_JMbs - 10
          }
          if (colnames(PAM250)[k_bef] %in% c('H','N','R','Q','K') & colnames(PAM250)[k] %in%  c('E','D')){
            score_JMbs <- score_JMbs - 10
          }
          if (colnames(PAM250)[k_bef] %in% c('H','N','R','Q','K') & colnames(PAM250)[k] %in%  c('S','T')){
            score_s <- score_s - 10
          }
          if (colnames(PAM250)[k_bef] %in% c('S','T') & colnames(PAM250)[k] %in% c('H','N','R','Q') ){
            score_sr <- score_sr - 10
          }
        }
        if (g==6){
          k<-which(seq3[i+(g-1)]== colnames(PAM250))
          k_bef <-which(seq3[i+(g-2)]== colnames(PAM250))
          score_l <- score_l + mean(PAM250[HDE,k])
          score_lr <- score_lr + mean(PAM250[HNQR,k])
          if (colnames(PAM250)[k] %in%  c('E','D') & colnames(PAM250)[k_bef] %in%  c('E','D')){
            score_JMbL <- score_JMbL + mean(PAM250[EDQ,k])
          }
          if (colnames(PAM250)[k] %in%  c('H','N','R','Q','K')){
            score_JMbL <- score_JMbL - 10
          }
          if (colnames(PAM250)[k_bef] %in% c('H','N','R','Q','K') & colnames(PAM250)[k] %in%  c('E','D')){
            score_JMbL <- score_JMbL - 10
          }
          if (colnames(PAM250)[k_bef] %in% c('H','N','R','Q','K') & colnames(PAM250)[k] %in%  c('S','T')){
            score_l <- score_l - 10
          }
          if (colnames(PAM250)[k_bef] %in% c('S','T') & colnames(PAM250)[k] %in% c('H','N','R','Q') ){
            score_lr <- score_lr - 10
          }
        }
        
      }
      scores[i,1,y-10]<- score_l
      scores[i,2,y-10]<-score_s
      scores[i,3,y-10]<- score_lr
      scores[i,4,y-10]<-score_sr
      scores[i,5,y-10]<-score_JMbL
      scores[i,6,y-10]<-score_JMbs
      
    }
  }
  
  q<-quantile(scores, probs = 0.98, na.rm = FALSE,
              +          names = TRUE)
  
  library('nnet')
  potential_sites <- matrix('',n_seqs,4)
  ind<-1
  
  for (j in 1:n_seqs){
    m[1]<-max(scores[,1,j])
    m[2]<-max(scores[,2,j])
    m[3]<-max(scores[,3,j])
    m[4]<-max(scores[,4,j])
    m[5]<-max(scores[,5,j])
    m[6]<-max(scores[,6,j])
    maximum <- which.is.max(m)
    if (maximum %in% c(1,2,3,4) & m[maximum]>=q){
      w1<-which(scores[,maximum,j]==m[maximum])
      seq<-s[j]
      seq2<-strsplit(seq,NULL)
      seq3<-unlist(seq2)
      potential_sites[ind,1]<- names[j+10]
      potential_sites[ind,2]<- paste(seq3[w1],seq3[w1+1],seq3[w1+2],seq3[w1+3],seq3[w1+4],seq3[w1+5], sep = "")
      potential_sites[ind,3]<- 'JMa'
      potential_sites[ind,4]<- w1
      ind <- ind+1
    }
    if (maximum %in% c(5,6) & m[maximum]>=q){
      w2<-which(scores[,maximum,j]==m[maximum])
      seq<-s[j]
      seq2<-strsplit(seq,NULL)
      seq3<-unlist(seq2)
      potential_sites[ind,1]<- names[j]
      potential_sites[ind,2]<- paste(seq3[w2[1]],seq3[w2[1]+1],seq3[w2[1]+2],seq3[w2[1]+3],seq3[w2[1]+4],seq3[w2[1]+5], sep="")
      potential_sites[ind,3]<- 'JMb'
      potential_sites[ind,4]<- w2[1]
      ind<-ind+1
    }
  }
  
  return(potential_sites[1:which(potential_sites[,1] == '')[1]-1,])
}