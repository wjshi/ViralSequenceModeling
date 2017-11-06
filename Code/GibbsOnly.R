source("Allfcns.R")
source("library_cmdline.R")
source("source.R")

# shortname<-cmdline.strings("shortname")
# copy<-cmdline.numeric("copy") #For possible multiple tests 
# load(paste("Tests/Data/TestData_",copy,".Rdata",sep=""))
shortname = "Test"

# a=1/25; n =300
# Pairs = c( "t1t2", "t1t3", "t2t3", "t1t3_D", "t2t3_D" ); PostD = 4
# 
# load(paste(name,"/Data.Rdata",sep=""))
for (copy in 11:100){
  name=paste("Tests/", shortname,copy,"_a25",sep="")
  for ( Ci in c(20, 40, 60, 80) )
  {
    sub=paste("bsub -q week -J Gibbs", shortname, copy,"_",Ci,"[1] -o ~/LSFouts/Gibbs.out R CMD BATCH --vanilla --args --name=",
              name," --isGibbs=",1," --replicate=\\\u0024LSB_JOBINDEX --copy=",copy," --Ci=",Ci," GibbsMC.R ", name,"/Gibbs",Ci,"\\\u0024LSB_JOBINDEX.out", sep="")
    system(sub) 
  }  
}


# Gibbsresult2<-function(name,rawdata,Rep=100)
# {
#   finalCs=matrix(NA,nrow=Rep,ncol=ncol(rawdata))
#   cputs=rep(NA,Rep)
#   for (rep in 1:Rep)
#   {
#     load(paste(name,"/GibbsOnly",Ci,"_",rep,".Rdata",sep=""))
#     finalCs[rep,]=CtoC(GibbsMC$C, g0=G0, originalData=rawdata)
#     cputs[rep]=cput
#   }
#   maxt=max(cputs)
#   Gibbs=list(FinalCs=finalCs, Runningtimes=cputs, maxCPUtime=maxt)
#   return(Gibbs)
# }
# 
# # Sys.sleep(60) 
# 
# 
# for (Ci in c(50,100) )
{
#   system(paste("bjobs -J Gibbs", shortname, copy,"_",Ci, " > JOB_",shortname, copy,sep=""))
#   jobsize=file.info(paste("JOB_",shortname, copy,sep=""))$size
#   
#   while(jobsize > 0){
#     Sys.sleep(60) 
#     system(paste("bjobs -J Gibbs", shortname, copy,"_",Ci, "* > JOB_",shortname, copy,sep=""))
#     jobsize=file.info(paste("JOB_",shortname, copy,sep=""))$size
#     print(jobsize)
#   }
#   system(paste("rm JOB_",shortname, copy, sep="" ))
  
  
#   Gibbs=Gibbsresult2(name,rawdata)
#   CPUtime=Gibbs$maxCPUtime
#   
#   folder = paste(name, "/Gibbs", Ci, "/", sep = ""); dir.create( folder )
#   HtMed( Gibbs$FinalCs, rawdata, a, PostD, n, folder, Control = T )
#   inffile = paste( folder, "mHts.Rdata", sep = "" )
#   resultfile =  paste( folder, "Local.Rdata", sep = "" ) 
#   load(inffile); Cuts = seq(0, max(mHts.pos), by = .001)
#   noise.test = sapply(Cuts, function(y){
#     result = freaksite0(inffile, Pairs, PostD, y, Control = T )
#     noise.l = c( length(result$Substitution), length(result$potential), 
#                  length(result$noise))
#     return( noise.l )
#   })
#   save( Pairs, PostD, Cuts, noise.test, file = resultfile ) 
#   
#   savefile =  paste( folder, "Final.Rdata", sep = "" )     
#   result = freak.result0( inffile, resultfile, savefile, delta = 3, Control = T, buffer = 100 )
#   
#   save(Gibbs, result, CPUtime, file = paste( folder, "Inference.Rdata",sep=""))
#   
# }

