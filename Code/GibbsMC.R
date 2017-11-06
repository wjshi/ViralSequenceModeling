source("Allfcns.R")
source("library_cmdline.R")




OnlyGibbs<-cmdline.numeric("isGibbs")
Ci<-cmdline.numeric("Ci")
name<-cmdline.strings("name")
copy<-cmdline.numeric("copy")


load(paste(name,"/Data.Rdata",sep=""))



if (OnlyGibbs==1)
{
  rep<-cmdline.numeric("replicate")
  set.seed(copy*20000+Ci*10+rep)
#   load(paste("Tests/Data/TestData_",copy,".Rdata",sep=""))
  t0=proc.time()
  initialC=kmeans(t(Data),Ci)$cluster
  GibbsMC=GibbsGeweke(Data,initialC,s=10000,a=1/25,bk=1000)
  t1=proc.time()-t0
  cput=t1[1]
  save(GibbsMC,cput,file=paste(name,"/GibbsOnly",Ci,"_",rep,".Rdata",sep=""))
  
  
#   save(initialC,GibbsMC,cput,file=paste(name,"/GibbsOnly",Ci,"_",rep,".Rdata",sep=""))  
}else
{
  a=1/(nrow(Data))^2
#  a=1/(nrow(Data)*2)
  
  set.seed(copy*20000+Ci*10+5+1/a)
  load(paste(name,"/GibbsResult.Rdata",sep=""))
  t0=proc.time()
  initialC=newC0s[Ci,]
  GibbsMC=GibbsGeweke(Data,initialC,10000,a,0)
  t1=proc.time()-t0
  cput=t1[1]
#  save(initialC,GibbsMC,cput,file="test6.Rdata")
  save(initialC,GibbsMC,cput,file=paste(name,"/GibbsMC",Ci,".Rdata",sep=""))
#   save(initialC,GibbsMC,file=paste(name,"/GibbsMC",Ci,".Rdata",sep=""))
}

