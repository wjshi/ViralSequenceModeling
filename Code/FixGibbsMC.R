source("Allfcns.R")
source("library_cmdline.R")




Ci<-cmdline.numeric("Ci")
name<-cmdline.strings("name")
copy<-cmdline.numeric("copy")
CORE<-cmdline.numeric("core")


load(paste(name,"/Data.Rdata",sep=""))
a=1/(nrow(Data))^2

set.seed(copy*20000+Ci*10+5+1/a)
load(paste(name,"/FixGibbsResult.Rdata",sep=""))
tempfile = paste(name, "/FixGibbsMC", Ci, "temp.Rdata", sep="")

t0=proc.time()
initialC=newC0s[Ci,]
GibbsMC=GibbsFixScan(Data, initialC, 1, a, Record=F, CORE=CORE, TempFile=tempfile)
t1=proc.time()-t0
cput=t1[1]
save(initialC,GibbsMC,cput,file=paste(name,"/FixGibbsMC",Ci,".Rdata",sep=""))

