library(coda)
library(parallel)


####################
#Plot raw frequency#
####################
freqplot1 <- function( DataFile, Sites, SaveFolder, M = 100 ){ #Three time points
  load( DataFile )
  n = ncol( rawdata ) / 3
  Colors = c("red","green","purple","blue")  
  Reads = c("A","C","G","T")
  topo.n = ceiling( max( colSums(rawdata) ) / M )
  Colors.M = topo.colors( topo.n, alpha=.8 )
  for (Site in Sites){
    sitedata = rawdata[ , (0 : 2) * n + Site ]
    freqlist = apply( sitedata, 2, function(x1){
      freq = x1/sum(x1); return( freq )
    } )   
    pdf( paste( SaveFolder, "/site", Site, ".pdf", sep = "" ) )
    #x: nucleotide
    #y: time
    readsums = colSums( sitedata )
    plot( 1:3, c( rep(0,2), 1 ), ylim = c(0,1.2), xlim =c(1,4), yaxt = "n", type = 'n', xlab = "passage", 
          ylab = "proportion", main = paste( "Site ", Site, sep = "" ), frame.plot = F )
    axis(2, at=seq(0,1.2,by=.2), labels = c( seq( 0, 1, by=.2 ), "total count" ) )
    read = sapply( 1:4, function(x){ points( freqlist[x,], col = Colors[x], cex =1, pch =x ) })
    read.sum = sapply( 1:3, function(z){ points( z, 1.2, pch = 15, cex = 3,
                                                 col = Colors.M[ceiling(readsums[z]/M)] )})
    legend( 3.5,.3, horiz = F, legend = Reads, col = Colors, cex =1, pch =1:4, bty = "n" )
    image.plot(legend.only = T, zlim = range(c(1,topo.n)*M), col=Colors.M, legend.shrink = .4,
               bigplot = c( .85,.72,1,1), midpoint = F, legend.line = 1, lwd= .5)
    dev.off()      
  }
}





#############################
#Random parameters (J=5) #
#############################
NucleoProb = function( J, J.main = 4 ){ 
  #J: number of possible reads; J.main: number of possible invariant reads 
  p = runif( J )
  main = rep( 0, J )
  r = runif(1)
  if( r > .5 ){ main[ sample.int( J.main, size = 1 ) ] = 100  
  }else if( r > .25 ){ main[ sample.int( J.main, size = 2 ) ] = c(25, 75)  
  }else{ main[ sample.int( J.main, size = 2 ) ] = c(50, 50) } 
  newp = rdirichlet( 1, p + main ) 
  return( newp )  
}

NewNucleoProb = function( oldp, read.main = 4 ){
  oldp.main = oldp[ 1:read.main ]
  oldp.2max = which( oldp.main == sort( oldp.main, decreasing = T )[2] )
  newp = oldp; newp[ oldp.2max ] = 10; newp = newp / sum( newp )
  return( newp )
}


##################################
#Consolidate homogenous positions#
##################################
consoData<-function( Data )
{
  J = nrow( Data )
  N0 = ncol( Data ) 
  C0 = rep( 0, N0 )
  for ( j in 1:J ){
    for ( i in 1:N0 ){
      if ( sum( Data[ -j, i ] == 0 ) == ( J - 1 ) ){ C0[ i ] = j }
    }
  }  
  g0 = list( NULL ) #Original indices for the consolidated sites
  na.list = list( NULL )
  newD = matrix( NA, nrow = J, ncol = J )
  for ( j in 1:J ){
    g0[[ j ]] = which( C0 == j )
    l = length( which( C0 == j ) )
    if ( l == 0 ){ na.list
    }else if ( l == 1 ){ newD[ , j ] = Data[ , g0[[ j ]] ]
    }else{ newD[ , j ] = rowSums( Data[ , g0[[ j ]] ] ) 
    }
  }
  g0[[ J+1 ]] = which( C0==0 )
  g0.length = lapply( g0, length )
  emptyg0 = which( g0.length == 0 )
  if ( length( emptyg0 ) > 0 ){ newD = newD[ , -emptyg0 ] }  
  newData = cbind( newD, Data[ , g0[[ J + 1 ]] ] )
  return( list( newData, g0 ) )
}

#########################################
#Create test data                       #          
#K=15; Substitutions=5; n=300; N=1200;  #
#########################################
createtestdata<-function(copy, K0 = 15, newK = 5, d = 5, n = 300, 
                         sitect = 5000, tests = 4, postD = 4 )
{
  PP = sapply( 1:K0, function(x){ return( NucleoProb( d ) ) } )
  newPP = apply( PP[ , 1:newK ], 2, NewNucleoProb )
  b = n / K0
  rawdata = matrix( NA, nrow = d, ncol = tests * n )
  for( i in 1:K0 ){
    for( t in 1:tests ){
      rawdata[ , ( b * ( i - 1 ) + 1 ) : ( b * i ) + n * ( t - 1 ) ] = 
        rmultinom( b, sitect, PP[ , i ] )  
    }
  }
  for( i in 1:newK ){
    for( t in postD ){
      rawdata[ , ( b * ( i - 1 ) + 1 ) + n * ( t - 1 ) ] = 
        rmultinom( 1, sitect, newPP[ , i ] )  
    }
  }
  save( rawdata, PP, newPP, tests, sitect, postD, 
        file = paste("Tests/Data/TestData_",copy,".Rdata",sep=""))
  return(paste("Test", copy," created!",sep=""))
}

#########################################################
#Thresholding the reads by percentage of the total count#
#########################################################
filterdata<-function(Data,threshold)
{
  threshold.list = lapply( 1 : ncol( Data ), function( x ){
    Col = Data[ , x ]; kills = which( Col < (sum( Col ) * threshold) ); Col[ kills ] = 0; return( Col ) } )
  newData = do.call( cbind, threshold.list)
  return( newData )
}


##################################
#Log posterior without normalizer#
##################################
LP<-function(Y,C,K,a) #a is the even weight in the Dirichlet prior, whose concentration=5a.
{  
  J=nrow(Y)
  lgempty=J*lgamma(a)-lgamma(J*a)
  if (is.null(dim(Y))==T) 
  {	
    lp=sum(lgamma(Y+a))-lgamma(sum(Y)+a*J)+(K-1)*lgempty
  }
  else 
  {	
    M=colSums(Y)
    kk=rep(0,K)
    for (k in 1:K)
    {
      ind=which(C==k)	
      if (length(ind)==0) {kk[k]=lgempty}
      else if (length(ind)==1) {kk[k]=sum(lgamma(Y[,ind]+a))-lgamma(M[ind]+J*a)}
      else {kk[k]=sum(lgamma(rowSums(Y[,ind])+a))-lgamma(sum(M[ind])+J*a)}
    }
    lp=sum(kk)
  }
  return(lp)
}

###################
#Geweke diagnostic#
###################
geweke.diag2 <- function (x, frac1 = 0.1, frac2 = 0.5) 
{
  if (frac1 < 0 || frac1 > 1) {
    stop("frac1 invalid")
  }
  if (frac2 < 0 || frac2 > 1) {
    stop("frac2 invalid")
  }
  if (frac1 + frac2 > 1) {
    stop("start and end sequences are overlapping")
  }
  if (is.mcmc.list(x)) {
    return(lapply(x, geweke.diag2, frac1, frac2))
  }
  x <- as.mcmc(x)
  xstart <- c(start(x), floor(end(x) - frac2 * (end(x) - start(x))))
  xend <- c(ceiling(start(x) + frac1 * (end(x) - start(x))), 
            end(x))
  y.variance <- y.mean <- vector("list", 2)
  for (i in 1:2) {
    y <- window(x, start = xstart[i], end = xend[i])
    y.matrix <- as.matrix(y)
    y.mean[[i]] <- apply(y.matrix, 2, mean)
    if (identical(sum(y.matrix==y.matrix[1])==length(y.matrix), TRUE)){
      y.variance[[i]] <- 0 
    }else{
      y.variance[[i]] <- spectrum0.ar(y)$spec/niter(y)
    }
  }
  z <- (y.mean[[1]] - y.mean[[2]])/sqrt(y.variance[[1]] + y.variance[[2]])
  out <- list(z = z, frac = c(frac1, frac2))
  class(out) <- "geweke.diag"
  return(out)
}



################################
#2-component Metropolis-Hasting#
################################
split<-function(data,G,s,a)
{
  Y=data[,G]
  N=length(G)
  if (N==1) {MC=list(left=G,right=NULL)}
  else if (sum(0==(Y-Y[,1]))==N*(nrow(Y))) {{MC=list(left=G,right=NULL)}}
  else{
    if (N==2){Co=c(1,2)
    }else{Co=kmeans(t(Y),2)[[1]]}
    lp=rep(NA,s)
    CC=matrix(NA,nrow=s,ncol=N)
    po=LP(Y,Co,2,a)
    t=0
    absZ=10
    while (absZ>2) #Metroplis checked with Geweke diagnostic
    {
      C=Co
      i=ceiling(N*runif(1))
      C[i]=C[i]%%2+1
      p=LP(Y,C,2,a)
      d=p-po
      if (runif(1)<exp(d)){ #Accept the proposal
        Co=C
        po=p
      }
      lp[t%%s+1]=po
      CC[t%%s+1,]=Co			
      if (t%%s+1==s) 
      {
        MClp=mcmc(lp)
        geweke=geweke.diag2(MClp)
        z=geweke$z
        if (is.na(z)){absZ=0 
        }else if (z=="Inf"){z=10
        }else{absZ=abs(z)}
        c=colSums(CC==1)/s
        MC=list(RunningC=CC, LP=lp, left=G[which(c>=.5)],right=G[which(c<.5)])
      }
      t=t+1
    }
    c=colSums(CC==1)/s
    MC=list(RunningC=CC, LP=lp, left=G[which(c>=.5)],right=G[which(c<.5)])
  }
  return(MC)
}

#############################
#Log ratio of two posteriors#
#############################
LR<-function(Y,C1,C2,i,a)
{
  A=which(C1==C1[i])
  B=which(C2==C2[i])
  M=colSums(Y)
  mi=M[i]
  J=nrow(Y)
  lr=rep(0,J)
  for (j in 1:J) 
  {
    y=Y[j,i]
    if (y>0) {lr[j]=-lbeta(sum(Y[j,B])-y+a,y)+lbeta(sum(Y[j,A])-y+a,y)}
  }
  LR=sum(lr)+lbeta(sum(M[B])-mi+a*J,mi)-lbeta(sum(M[A])-mi+a*J,mi)
  return(LR)
}

#########################################
#Bundle data according to tree end nodes#
#########################################
grpData<-function(Y,Gp)
{
  K=length(Gp)
  Z=matrix(NA,nrow=nrow(Y),ncol=K)
  for (k in 1:K) 
  {
    if (length(Gp[[k]])==1) {Z[,k]=Y[,Gp[[k]]]}
    else {Z[,k]=rowSums(Y[,Gp[[k]]])}	
  }
  return(Z)
}

###################
#Combine groups#
###################
blockcomb<-function(Z,Gp,s,a,bk) #bk is the burnin number in terms of multiples of K.
{
  K=length(Gp)
  if (K==1) {print("one group")
  }else {
    lp=rep(NA,s)
    BCC=matrix(NA,nrow=s,ncol=K)
    BC0=seq(1,K,1)
    p0=LP(Z,BC0,K,a)
    t=0
    absZ=10
    burn=bk*K
    maxtotal=burn+100000*K #Set an upper bound for number of iterations.
    while (absZ>=2) #Check stationary using Geweke diagnostic.
    {
      BC=BC0
      i=ceiling(K*runif(1))
      Bci=(ceiling((K-1)*runif(1))+BC[i]) %% K
      if (Bci==0) {BC[i]=K}
      else{BC[i]=Bci}
      d=exp(LR(Z,BC0,BC,i,a))
      u=runif(1)
      if (t<burn){
        if (u<d) {BC0=BC}
        else {BC0=BC0}	
      }else{
        tb=t-burn
        if (u<d)	
        {
          BC0=BC
          p0=LP(Z,BC,K,a)
        }
        else 
        {
          BC0=BC0
          p0=p0
        }
        lp[tb%%s+1]=p0
        BCC[tb%%s+1,]=BC0			
        
        if (tb%%s+1==s) 
        {
          MClp=mcmc(lp)
          geweke=geweke.diag2(MClp)
          z=geweke$z
          if (is.na(z)){absZ=0 
          }else if (z=="Inf"){absZ=10
          }else{absZ=abs(z)}
        }
      }
      t=t+1
      if (t==maxtotal)
      {
        absZ=0
        print(paste("Stopped after ",maxtotal," iterations",sep=""))
      }
    }   
    MC=list(RunningBC=BCC, LP=lp,t=t,zscore=z)
    return(MC)
  }
}

##############################
#treeGp&BC->DatafinalGp&DataC#
##############################
newGp<-function(treeGp,BC,consoData)
{
  n=dim(consoData)[2]
  C=rep(NA,n)
  K=length(BC)
  nK=length(unique(BC))
  if (nK==K) 
  {
    for (k in 1:K){C[ treeGp[[k]] ]=k} 
    return(list(group=treeGp, index=C))
  }
  else
  {
    G=list(rep(NA,nK))
    for (k in 1:K) {C[treeGp[[k]]]=BC[k]}
    for (j in 1:nK) {G[[j]]=which(C==unique(BC)[j])}
    for (i in 1:nK){C[G[[i]]]=i} 
    return(list(group=G, index=C))
  }	
}

###############################
#runningBC=>GibbsinitialDataCs#
###############################
GibbsC0<-function(runningBC,treeGp,consoData)
{
  if (is.null(dim(runningBC))){return(newGp(treeGp,runningBC,consoData)$index)
  }else{
    GibbC<-function(BC){return(newGp(treeGp,BC,consoData)$index)}
    C0s=apply(runningBC,1,GibbC)
    return(t(C0s))
  }
}

#############################
#rename C
#############################
renameC<-function(C)
{
  newC = rep( NA,length( C ) )
  uniqueC = unique( C )
  for ( i in 1:length(uniqueC) ){ 
    newC[ which( C == uniqueC[i] ) ] = i }
  return( newC )	
}

###############################
#sums of y & m with indicators#
###############################
myfcn<-function(data,C,K)
{
  J=nrow(data)
  csY=colSums(data)
  yInd=matrix(NA,nrow=J,ncol=K)	
  for (k in 1:K)
  {
    Ind=which(C==k)
    l=length(Ind)
    if (l==0){yInd[,k]=rep(0,J)
    }else{ for (j in 1:J) {yInd[j,k]=sum(data[j,Ind])} }
  }
  mInd=colSums(yInd)
  my=list(sMind=mInd,sYind=yInd)	
  return(my)
}

##########################
#Choose a C from runningC#
##########################
convC<-function(CC,lp)
{
  aclp=abs(lp-mean(lp))
  state<-which(aclp==min(aclp))
  if (length(state)==1) {return(CC[state,])}
  else 
  {
    s=state[1]
    return(CC[s,])
  }
}

convC.median<-function(CC,lp)
{
  state<-which(lp==median(lp))
  if (length(state)==1) {return(CC[state,])}
  else 
  {
    s=state[1]
    return(CC[s,])
  }
}
###############
#Gibbs
##############
LR.Gibbs<-function(Y,C,i,ki,a) #i is the sequence index, ki is the new group label for position i. 
{
  A=which(C==C[i])
  tempC=C
  tempC[i]=ki
  B=which(tempC==ki)
  M=colSums(Y)
  mi=M[i]
  J=nrow(Y)
  lr=rep(0,J)
  for (j in 1:J) 
  {
    y=Y[j,i]
    if (y>0) {lr[j]=-lbeta(sum(Y[j,B])-y+a,y)+lbeta(sum(Y[j,A])-y+a,y)}
  }
  LR=sum(lr)+lbeta(sum(M[B])-mi+a*J,mi)-lbeta(sum(M[A])-mi+a*J,mi)
  return(LR)
}

GibbsGeweke<-function(Data,initialC,s,a,bk)
{
  n=length(initialC)
  K=length(unique(initialC))
  Kseq=seq(1,K,1)
  t=0
  C0=initialC  
  LPs=rep(NA,s)
  Cs=matrix(NA,nrow=s,ncol=n)
  absZ=10
  burn=K*bk
  z=0
  maxtotal=burn+1000000*K
  while(absZ>=2)
  {
    i=ceiling(n*runif(1))
    lrs=sapply(Kseq, LR.Gibbs, Y=Data, C=C0, i=i, a)
    rs=exp(lrs-max(lrs))
    pp=rs/sum(rs)
    cumPi=cumsum(pp)
    u=runif(1)
    C0[i]=sum(u>cumPi)+1
    if(t>=burn)
    {
      tb=t-burn
      LPs[tb%%s+1]=LP(Data,C0,K,a)
      Cs[tb%%s+1,]=C0
      if(tb%%s+1==s)  
      {
        MClp=mcmc(LPs)
        geweke=geweke.diag2(MClp)
        z=geweke$z
        if (is.na(z)){absZ=0 
        }else if (z=="Inf"){absZ=10
        }else{absZ=abs(z)}   
        C=convC(Cs,LPs)
        MC=list(RunningC=Cs, LP=LPs, t=t, zscore=z, C=C)
      }
    }
    t=t+1
    if (t==maxtotal)
    {
      absZ=0
      print(paste("Stopped after ",maxtotal," iterations",sep=""))
    }
  }
  C=convC(Cs,LPs)
  MC=list(RunningC=Cs, LP=LPs, t=t, zscore=z, C=C)
  return(MC)
}

GibbsRandomScan<-function(Data,initialC,Scan,a,bk)
{
  n=length(initialC)
  K=length(unique(initialC))
  Kseq=seq(1,K,1)
  t=0
  C0=initialC  
  LPs=rep(NA,Scan)
  Cs=matrix(NA,nrow=Scan,ncol=n)
  burn=K*bk
  maxtotal=burn+Scan
  while(t<maxtotal)
  {
    i=ceiling(n*runif(1))
    lrs=sapply(Kseq, LR.Gibbs, Y=Data, C=C0, i=i, a)
    rs=exp(lrs-max(lrs))
    pp=rs/sum(rs)
    cumPi=cumsum(pp)
    u=runif(1)
    C0[i]=sum(u>cumPi)+1
    if(t>=burn)
    {
      tb=t-burn
      LPs[tb+1]=LP(Data,C0,K,a)
      Cs[tb+1,]=C0
    }
    t=t+1
  }
  C=convC(Cs,LPs)
  MC=list(RunningC=Cs, LP=LPs, t=t, C=C, TotalScan=t)
  return(MC)
}

GibbsFixScan<-function(Data,initialC,Scan,a, Record = F, CORE = 2, TempFile)
{
  N=length(initialC)
  K=length(unique(initialC))
  Kseq=seq(1,K,1)
  C0=initialC
  LPs = Cs = NULL
  if( Record == T ){
    LPs=rep(NA,Scan)
    Cs=matrix(NA,nrow=Scan,ncol=N)  
  }
  Scan0 = 0
  while(Scan0 < Scan){
    for( ii in 1:N){
      lrs=simplify2array(mclapply(Kseq, LR.Gibbs, Y=Data, C=C0, i=ii, a = a, mc.cores = CORE))
      rs=exp(lrs-max(lrs))
      pp=rs/sum(rs)
      cumPi=cumsum(pp)
      u=runif(1)
      C0[ii]=sum(u>cumPi)+1  
      if( ii %% 100 == 0 ){
        save( ii, Scan0, C0, K, file = TempFile)
      }
    }
    Scan0 = Scan0 + 1
    if( Record == T ){
      LPs[Scan0] = LP(Data,C0,K,a)
      Cs[Scan0, ] = C0  
    }
  }
  MC=list(ScanC=Cs, ScanLP=LPs, Scan=Scan, C = C0)
  return(MC)
}

GibbsFixScan.cont<-function(TempFile, Data, Scan, a, CORE = 2)
{
  load(TempFile)
  N=length(C0)
  Kseq=seq(1,K,1)
  LPs = Cs = NULL
  ii0 = ii + 1
  while(Scan0 < Scan){
    for( newii in ii0:N){
      lrs = simplify2array(mclapply(Kseq, LR.Gibbs, Y = Data, C = C0, i = newii, 
                                    a = a, mc.cores = CORE))
      rs=exp(lrs-max(lrs))
      pp=rs/sum(rs)
      cumPi=cumsum(pp)
      u=runif(1)
      C0[newii]=sum(u>cumPi)+1  
      if( newii %% 100 == 0 ){
        ii = newii
        save( ii, Scan0, C0, K, file = TempFile)
      }
    }
    Scan0 = Scan0 + 1
  }
  MC=list(ScanC=Cs, ScanLP=LPs, Scan=Scan, C = C0)
  return(MC)	
}


###########################
#consoDataC=>originalDataC#
###########################
CtoC<-function(consoC,g0,originalData)
{
  J=nrow(originalData)
  C=rep(NA,ncol(originalData))
  gl = lapply( g0, length )
  gl = do.call( rbind, gl )
  glJ = which( gl[ 1:J ] > 0 )
  jj = length( glJ )
  if( jj > 0 ){ 
    for (i in 1 : jj ) {
      j = glJ[ i ]
      C[g0[[j]]]=consoC[i]    
    }
  }
  s=seq(1,length(g0[[J+1]]))
  C[g0[[J+1]][s]]=consoC[s+jj]  
  return(C)
}

#######################
#Gather Gibbs results#
#######################
Gibbsresult2<-function(C0s,name,rawdata,GibbsMC=T)
{
  if( GibbsMC == T ){
    cputs = maxt = NULL
    GibbsFiles = Sys.glob(paste(name,"/FixGibbsDone/FixGibbsMC*.Rdata",sep="")) 
    finalCs = matrix(NA,nrow=length(GibbsFiles),ncol=ncol(rawdata))
    i = 1
    for (filename in GibbsFiles)
    {
      load(filename)
      finalCs[i,]=CtoC(renameC(GibbsMC$C), g0=G0, originalData=rawdata)
      cputs[i]=cput
      i = i + 1
    }
    maxt=max(cputs)
  }else{
    finalCs = t( apply( C0s, 1, CtoC, g0=G0, originalData=rawdata ) )
    cputs = maxt = NULL
  }  
  Gibbs=list(FinalCs=finalCs, Runningtimes=cputs, maxCPUtime=maxt)
  return(Gibbs)
}

################################
#Transformed Hellinger distance#
#ln(1-ln(1-H^2))               #
################################
HtMatrix<-function(data,C,a)
{
  K=length(unique(C))
  if (K==1)
  {
    return(paste("Single group, Ht matrix is zero!",sep=""))
  }else{
    Ht=diag(K,x=0)
    indsum=myfcn(data,C,K)
    J=nrow(data)
    sMind=indsum[[1]]
    sYind=indsum[[2]]
    for (i in 1:(K-1))
    {
      for (j in (i+1):K)
      {
        L1j=sYind[,i]+a
        L2j=sYind[,j]+a
        L1=sMind[i]+a*J
        L2=sMind[j]+a*J
        logfrac=sum(lgamma((L1j+L2j)/2))-lgamma((L1+L2)/2)-
          .5*(sum(lgamma(L1j)+lgamma(L2j))-lgamma(L1)-lgamma(L2))
        fH=log(1-logfrac)
        Ht[i,j]=fH
        Ht[j,i]=fH
      }
    }
  } 
  return(Ht)
}

HtDiv<-function(data,runningC,a) #Three time points pairwise Ht distances
{ 
  L=nrow(runningC)
  n=ncol(data)/3
  Ht12=matrix(NA,nrow=L,ncol=n)
  Ht13=matrix(NA,nrow=L,ncol=n)
  Ht23=matrix(NA,nrow=L,ncol=n)
  for (l in 1:L)
  {
    Cl=runningC[l,]
    Htl=HtMatrix(data,Cl,a)
    for (i in 1:n)
    {
      Ht12[l,i]=Htl[Cl[i],Cl[i+n]]
      Ht13[l,i]=Htl[Cl[i],Cl[i+n*2]]
      Ht23[l,i]=Htl[Cl[i+n],Cl[i+n*2]]
    }
  }
  return(list(Ht12,Ht13,Ht23))
}

HtList<-function(data,runningC,a){ 
  L=nrow(runningC)
  Htls=lapply(1:L, function(x){return(HtMatrix(data,runningC[x,],a))} )
  return(Htls)
}

pairHt<-function(Htls,t1,t2,runningC,n)
{
  L=nrow(runningC)
  if (t1==t2){Ht=matrix(0,nrow=L,ncol=n)}
  else{
    Ht=matrix(NA,nrow=L,ncol=n)
    subCs=runningC[,c(((t1-1)*n+1):(t1*n),((t2-1)*n+1):(t2*n))]  
    for (l in 1:L)
    {
      Htl=Htls[[l]]
      subC=subCs[l,]
      Ht[l,]=sapply(1:n,function(x){return(Htl[subC[x],subC[x+n]])})
    }
  }
  return(Ht)
}

HtMed <- function( GibbsCs, Rawdata, a, PostD, n, FolderName, Control = T, mHtsFile = "mHts.Rdata" ){
  Finalcs = t( apply( GibbsCs, 1, renameC ) )
  Htlist = HtList( Rawdata, Finalcs, a )
  Hts.med = list( t1t2 = NULL, N = NULL, D = NULL ) 
  ht1 = pairHt( Htlist, 1, 2, Finalcs, n )
  Hts.med[[1]] = apply( ht1, 2, quantile, .5 )
  D.length = length( PostD )
  
  if( Control == T){
    for( i in 2:3 ){
      for( j in 1:2 ){
        Hts.med[[i]][[j]] = matrix( NA, nrow = D.length, ncol = n )
        for ( t in 1:D.length ){
          k = (i-2) * D.length + 2 + t
          htijt = pairHt( Htlist, j, k, Finalcs, n )
          Hts.med[[ i ]][[ j ]][ t, ] = apply( htijt, 2, quantile, .5 )
        }
      }
    }
    mHts.pool = c( Hts.med[[1]], as.vector( t( do.call( rbind, Hts.med[[2]] ) ) ),
                   as.vector( t( do.call( rbind, Hts.med[[3]] ) ) ) ) 
  }else{
    Hts.med[[2]] = NA
    for( j in 1:2 ){
      Hts.med[[3]][[j]] = matrix( NA, nrow = D.length, ncol = n )
      for ( t in 1:D.length ){
        k = 2 + t
        htijt = pairHt( Htlist, j, k, Finalcs, n )
        Hts.med[[ 3 ]][[ j ]][ t, ] = apply( htijt, 2, quantile, .5 )
      }
    }
    mHts.pool = c( Hts.med[[1]], as.vector( t( do.call( rbind, Hts.med[[3]] ) ) ) ) 
  }
  mHts.pos = mHts.pool[ mHts.pool > 0 ]
  
  save( Finalcs, Hts.med, mHts.pool, mHts.pos, 
        file = paste( FolderName, mHtsFile, sep = "" ) )
} 

freaksite0 <- function( InfFile, Pairs, PostD, Cut, Control = T ){ 
  load( InfFile )
  datasets = length( Pairs ); n = length( mHts.pool )/datasets; D.length = length( PostD )
  if( Control == T ){
    Hts.NoD = rbind( Hts.med[[1]], do.call( rbind, Hts.med[[2]] ) )  
    maxHts.NoD = apply( Hts.NoD, 2, max )
  }else{
    maxHts.NoD = Hts.NoD = Hts.med[[1]]
  }
  Hts.Dend = rbind( Hts.med[[3]][[1]][D.length, ], Hts.med[[3]][[2]][D.length, ] )
  minHts.Dend = apply( Hts.Dend, 2, min ) #; maxHts.Dend = apply( Hts.Dend, 2, max )
  D.Site = which( minHts.Dend > Cut ); DN.Site = intersect( which( minHts.Dend < maxHts.NoD ), which( maxHts.NoD > Cut ) ) 
  #     DN.Site = which( maxHts.NoD > Cut )
  S = TailT = TailT.noise = NULL#; Cut2 = "Site Dependent"
  if( length( D.Site ) > 0 ){
    S = setdiff( D.Site, DN.Site )
    if( length( S ) > 0 ){
      Ht.S = t(sapply(S, function(x){ return( mHts.pool[x + (0:(datasets - 1)) * n] ) }))
      colnames( Ht.S ) = Pairs
      TailT = cbind( data.frame( Site = S ), Ht.S )
    } 
  }
  if( length( DN.Site ) > 0 ){
    Ht.noise = t( sapply( DN.Site, function(x){
      return( mHts.pool[x + (0:(datasets - 1)) * n] ) 
    }))
    colnames( Ht.noise ) = Pairs
    TailT.noise = cbind( data.frame( Site = DN.Site ), Ht.noise )
  }
  S.list = list( Substitution = S, table.S = TailT, potential = D.Site, noise = DN.Site, 
                 table.noise = TailT.noise, cutoff = Cut )  
  return( S.list )
} 

freak.result0 <- function( InfFile, ResultFile, SaveFile, delta = 3, Control = T, buffer = 100 ){
  load( InfFile ); load( ResultFile ) #Provide mHts.pool, Pairs, PostD, noise.test
  lwr = which(Cuts > quantile(mHts.pool, .9))[1]; noise = noise.test[ 3, lwr:ncol(noise.test) ]
  noise.unique.length = length(unique(noise))
  if( noise.unique.length < delta * 2 ){ 
    delta = 1 
  }
  if( noise.unique.length == 1 ){
    cutoff.s = cutoff.s0 = 1; noi.l = length(noise); noi.ldiff = 0
  }else{
    noi.l = sapply( 0:max(noise), function(x){return(sum(noise == x))} )
    delta2 = delta*2-1
    noi.ldiff = sapply( 1:( length(noi.l) - delta2 ), function(x){
      del = sum(noi.l[x + delta:delta2]) - sum(noi.l[x + 0:(delta - 1)])
      return(del)
    })
    cutoff.s0 = max(which(noi.ldiff == min(noi.ldiff)))
    cutoff.s = cutoff.s0 + delta - 1
  }    
  if(cutoff.s0 == 1){
    delta0 = delta
    while( noi.l[delta0] < buffer){
      cutoff.s = cutoff.s - 1 
      delta0 = delta0 - 1
      if(delta0 == 0){print("Too much noise!")}
    }     
  }
  d = Cuts[min(which(noise.test[3,] == cutoff.s - 1 ) ) + buffer ]     
  result = freaksite0( InfFile, Pairs, PostD, d, Control )
  save(delta, noi.l, noi.ldiff, cutoff.s, d, result, buffer, file = SaveFile)
  return(result)
}

freak.result1 <- function( InfFile, ResultFile, SaveFile, delta = 3, Control = T, alpha = .5 ){
  load( InfFile ); load( ResultFile ) #Provide mHts.pool, Pairs, PostD, noise.test
  lwr = which(Cuts > quantile(mHts.pool, .9))[1]; noise = noise.test[ 3, lwr:ncol(noise.test) ]
  noise.unique.length = length(unique(noise))
  if( noise.unique.length < delta * 2 ){ 
    delta = 1 
  }
  if( noise.unique.length == 1 ){
    cutoff.s = cutoff.s0 = 1; noi.l = length(noise); noi.ldiff = 0
  }else{
    noi.l = sapply( 0:max(noise), function(x){return(sum(noise == x))} )
    delta2 = delta*2-1
    noi.ldiff = sapply( 1:( length(noi.l) - delta2 ), function(x){
      del = sum(noi.l[x + delta:delta2]) - sum(noi.l[x + 0:(delta - 1)])
      return(del)
    })
    cutoff.s0 = max(which(noi.ldiff == min(noi.ldiff)))
    cutoff.s = cutoff.s0 + delta - 1
  }    
  if(cutoff.s0 == 1){
    delta0 = delta
    buffer = mean( noi.l[1:delta0] ) * alpha
    while( noi.l[delta0] < buffer){
      cutoff.s = cutoff.s - 1 
      delta0 = delta0 - 1
      buffer = mean( noi.l[1:cutoff.s] ) * alpha
    }     
  }else{ buffer = mean( noi.l[(cutoff.s - delta + 1):cutoff.s] ) * alpha }                    
  d = Cuts[min(which(noise.test[3,] == cutoff.s - 1 ) ) + buffer ]     
  result = freaksite0( InfFile, Pairs, PostD, d, Control )
  save(delta, buffer, noi.l, noi.ldiff, cutoff.s, d, result, alpha, file = SaveFile)
  return(result)
}


plotNiDi <- function( InfFile, ResultFile, SaveFile, SavePlot = F, PlotName = NULL,
                      Main = "Summary Stats Plot", d0.value = T, Legend = T, 
                      Control = T ){
  load(InfFile); load(ResultFile); load(SaveFile)
  datasets = length( Pairs ); n = length( mHts.pool )/datasets; D.length = length( PostD )
  if( Control == T ){
    Hts.NoD = rbind( Hts.med[[1]], do.call( rbind, Hts.med[[2]] ) )  
    maxHts.NoD = apply( Hts.NoD, 2, max )
  }else{ maxHts.NoD = Hts.NoD = Hts.med[[1]] }
  Hts.Dend = rbind( Hts.med[[3]][[1]][D.length, ], Hts.med[[3]][[2]][D.length, ] )
  minHts.Dend = apply( Hts.Dend, 2, min ) 
  signals = result$Substitution
  Ylim = c(0,max(mHts.pos)+2)
  if( SavePlot == F ){
    plot(minHts.Dend, cex= .5, col = "green", ylim = Ylim,
         font.lab = 2, main = Main, xlab = "Genome site index", ylab = "Ht" )
    points(signals, minHts.Dend[signals], cex= 1.5, col = "red")
    points(maxHts.NoD, pch = 3, cex= .5, col = "blue")
    abline(h=d, lty = 2)
    if(d0.value == T){text(n*.05, d+.5, substitute('d'[0] == dd, list(dd = d)), cex =.8 )
    }else{ text(0, d+.5, expression(bold('d'[0] )), cex =.8 )}
    if(Legend == T){
      legend(n*.1, max(Ylim)*1.1, legend = c("Ht(Di) ", " Ht(Ni)", "Signal"), pch = c(1,3,1), 
             pt.cex = c( .5, .5, 1.5), cex = .8, 
             col = c("green", "blue", "red"), bty = "n", horiz = T)  
    } 
  }else{
    pdf( PlotName )
    plot(minHts.Dend, cex= .5, col = "green", ylim = Ylim,
         font.lab = 2, main = Main, xlab = "Genome site index", ylab = "Summary statistics Ht" )
    points(signals, minHts.Dend[signals], cex= 1.5, col = "red")
    points(maxHts.NoD, pch = 3, cex= .5, col = "blue")
    abline(h=d, lty = 2)
    if(d0.value == T){text(n*.05, d+.5, substitute('d'[0] == dd, list(dd = d)), cex =.8 )
    }else{ text(0, d+.5, substitute('d'[0] ), cex =.8 )}
    if(Legend == T){
      legend(n*.1, max(Ylim)*1.1, legend = c("Ht(Di) ", " Ht(Ni)", "Signal"), pch = c(1,3,1), 
             pt.cex = c( .5, .5, 1.5), cex = .8, 
             col = c("green", "blue", "red"), bty = "n", horiz = T)
    }
    dev.off()
  }  
}





#################################
#Extract inference result
##################################
TestResult <- function( TestFiles, mut.true = c(1, 21, 41, 61, 81) ){ 
  mut.results = sapply( TestFiles, function(x){ 
    load( x )
    freaks = as.numeric( result$Substitution )
    if( identical(freaks, mut.true) ){ return( c(1, 0, 0, 0) )                  
    }else{
      diff.length = length( setdiff(freaks, mut.true) )
      if (diff.length == 0){ return(c(0, 1, 0, 0))
      }else{
        diff.length2 = length( setdiff(mut.true, freaks) )
        if (diff.length2 == 0){ return(c(0, 0, 1, 0))
        }else{ return(c(0, 0, 0, 1)) } 
      }
    }
  } #c(PR, FN, FP, FPN)
  )  
  return( rowSums(mut.results) )
}
TestTime <- function( TimeFiles ){
  testtime = sapply( TimeFiles, function(x){
    load(x); return(CPUtime)
  })
  return(testtime)
}

############################

