##
## survivalROC_NEW.R
##
## AUTHOR:  P. Heagerty
##
## DATE:  98/08/18
##
## (last modified:  98/11/17)
## Comments added by : P. SAHA
## DATE: 06/04/19
## Function name changed from roc.KM.calc() to survivalROC() as of 06/05/17
##---------------------------------------------------------------------------
survivalROC_NEW<-function( Stime, status, marker, entry=NULL,
                       predict.time, cut.values=NULL,
                       method="NNE", lambda=NULL, span=NULL,
                       window="symmetric" )
{
  ## DATE: July 12, 2006
  ## one of the option of method is changed from "smooth" to "NNE"
  ##
  ## PURPOSE:  calculations for Kaplan-Meier based ROC curve
  ##
  ##
  ###
  ### drop any missing
  ###
  ##
  times=Stime
  ## changed input from times to Stime, June 1, 2006, by Paramita Saha
  ## to make both the functions similar
  x <- marker
  
  if( is.null(entry) ) entry <- rep( 0, length(times) )
  ##
  bad <- is.na(times) | is.na(status) | is.na(x) | is.na(entry)
  entry <- entry[!bad]
  times <- times[!bad]
  status <- status[!bad]
  x <- x[!bad]
  if( sum(bad)>0 ) cat(paste("\n", sum(bad),
                             "records with missing values dropped. \n") )
  ## 
  if( is.null(cut.values) ) cut.values <- unique(x)
  cut.values <- cut.values[ order(cut.values) ]
  ncuts <- length(cut.values)
  ##
  ###
  ### sort the times
  ###
  ##
  ooo <- order( times )
  times <- times[ooo]
  status <- status[ooo]
  x <- x[ooo]
  ##
  ###
  ### overall survival probability
  ###
  ##
  s0 <- 1.0
  unique.t0 <- unique( times )
  unique.t0 <- unique.t0[ order(unique.t0) ]
  ## print( unique.t0 )
  n.times <- sum( unique.t0 <= predict.time )
  for( j in 1:n.times ){
    n<-sum( entry <= unique.t0[j] & times >= unique.t0[j] )
    d<-sum( (entry <= unique.t0[j])&(times==unique.t0[j])&(status==1) )
    if( n>0 ) s0<-s0*( 1 - d/n )	
  }
  s.pooled <- s0
  ## print( s.pooled )
  ##
  ###
  ### compute the sensitivity and specificity for each cut.value
  ###
  ##
  roc.matrix <- matrix( NA, ncuts, 2 )
  ##
  ## changed to strictly greater 98/11/17
  ##
  roc.matrix[ncuts,1] <- 0.0
  roc.matrix[ncuts,2] <- 1.0
  ##
  if( method=="KM" ){
    ##
    for( c in 1:(ncuts-1) ){
      ##
      s0 <- 1.0
      ##
      ## changed to strictly greater 98/11/17
      ##
      subset <- as.logical( x > cut.values[c] )
      ##
      e0 <- entry[subset]
      t0 <- times[subset]
      c0 <- status[subset]
      if( !is.null(t0) ){
        unique.t0 <- unique( t0 )
        unique.t0 <- unique.t0[ order(unique.t0) ]
        n.times <- sum( unique.t0 <= predict.time )
        if( n.times>0 ){
          for( j in 1:n.times ){
            n<-sum( e0 <= unique.t0[j] & t0 >= unique.t0[j] )
            d<-sum( (e0<=unique.t0[j])& (t0==unique.t0[j])&(c0==1) )
            if( n>0 ) s0<-s0*( 1 - d/n )	
          }
        }
      }
      p0 <- mean(subset) 
      roc.matrix[ c, 1 ] <-  (1-s0) * p0 / (1-s.pooled)
      roc.matrix[ c, 2 ] <-  1 - s0 * p0 / s.pooled
      ##
    }
  }## end-of method=="KM"
  ##
  ##
  if( method=="NNE" ){
    if( is.null(lambda) & is.null(span) ){
      ## cat("method = smooth requires either lambda or span! \n")
      cat("method = NNE requires either lambda or span! \n")
      stop(0)
    }
    x.unique <- unique(x)
    x.unique <- x.unique[ order(x.unique) ]
    S.t.x <- rep( 0, length(x.unique) )
    t.evaluate <- unique( times[status==1] )
    t.evaluate <- t.evaluate[ order(t.evaluate) ]
    t.evaluate <- t.evaluate[ t.evaluate <= predict.time ]
    for( j in 1:length(x.unique) ){
      if( !is.null(span) ){
        if( window=="symmetric" ){
          ##
          ###
          ### symmetric
          ###
          ##
          ddd <- ( x - x.unique[j] )
          n <- length(x)
          ddd <- ddd[ order(ddd) ]
          index0 <- sum( ddd< 0 ) + 1
          index1 <- index0 + trunc(n*span +.5)
          
          ## Date: May 5, 2008
          ## Someone sent us an email saying that the same span
          ## produced very different results in .R vs .C
          ## code. We found out that .R was using interval
          ## length = n*span/2 on either side while .C was using
          ## length = n*span+.5 on either side, so the interval
          ## length considered in .C was double the length for
          ## .R. Today, we fixed the .R code from the following
          ## line to the line above.
          ## index1 <- index0 + trunc(n*span +.5)
          
          if( index1>n ) index1 <- n
          lambda <- ddd[ index1 ]
          wt <- as.integer( ( (x-x.unique[j]) <= lambda ) & 
                              ( (x-x.unique[j]) >= 0 ) )
          index0 <- sum( ddd<= 0 ) 
          index2 <- index0 - trunc(n*span/2)
          if( index2<1 ) index2 <- 1
          lambda <- abs( ddd[ index1 ] )
          set.index<- ( (x-x.unique[j]) >= -lambda ) & 
            ( (x-x.unique[j]) <= 0 ) 
          wt[ set.index ] <- 1
          
        }
        if( window=="asymmetric" ){
          ##
          ###
          ### asymmetric
          ###
          ##
          ddd <- ( x - x.unique[j] )
          n <- length(x)
          ddd <- ddd[ order(ddd) ]
          index0 <- sum( ddd< 0 ) + 1
          index <- index0 + trunc(n*span)
          if( index>n ) index <- n
          lambda <- ddd[ index ]
          wt <- as.integer( ( (x-x.unique[j]) <= lambda ) & 
                              ( (x-x.unique[j]) >= 0 ) )
        }
      }else{
        ## wt <- as.integer( abs(x-x.unique[j]) <= lambda )
        wt <- exp( -(x-x.unique[j])^2 / lambda^2 )
      }
      s0 <- 1.0
      for( k in 1:length(t.evaluate) ){
        n <- sum( wt*(entry<=t.evaluate[k])&(times >= t.evaluate[k]) )
        d <- sum( wt*(entry<=t.evaluate[k])&(times==t.evaluate[k])*
                    (status==1) )
        if( n > 0 ) s0 <- s0 * ( 1 - d/n )
      }
      S.t.x[j] <- s0
    }
    S.all.x <- S.t.x[ match( x, x.unique ) ]
    n <- length( times )
    S.marginal <- sum( S.all.x )/n
    for( c in 1:(ncuts-1) ){
      ##
      ## changed to strictly greater 98/11/17
      ##
      p1 <- sum( x > cut.values[c] )/n
      Sx <- sum( S.all.x[ x > cut.values[c] ] )/n
      roc.matrix[ c, 1 ] <- (p1-Sx)/(1-S.marginal)
      roc.matrix[ c, 2 ] <-  1 - Sx/S.marginal
    }
  }## end-of method=="NNE"
  ##
  sensitivity=roc.matrix[,1]
  specificity=roc.matrix[,2]
  x <- 1 - c( 0.0, specificity )
  y <- c( 1.0, sensitivity )
  ##
  n <- length( x )
  ##
  dx <- x[-n] - x[-1]
  mid.y <- ( y[-n] + y[-1] )/2
  ##
  ## print( sum(dx) )
  area <- sum( dx*mid.y )
  
  list(cut.values=c(-Inf, cut.values),
       ## renamed x as cut.values and changed from (NA,unique.x) to this one. June 1, 2006, By Paramita Saha
       #Changed to TP and FP with 1.0 value added 06/05/15
       sensitivity=roc.matrix[,1],
       TP=y,
       specificity=roc.matrix[,2],
       FP=x,
       predict.time=predict.time, ## predict.time included in the return object, June 1, 2006
       Survival = s.pooled,
       AUC=area)
}
