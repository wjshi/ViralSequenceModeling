source("library_cmdline.R")
source("source.R")
source("Allfcns.R")

shortname<-cmdline.strings("shortname") #For simulated data, shortname = Test 
copy<-cmdline.numeric("copy") #For possible multiple tests 
CORE<-cmdline.numeric("core") #Number of cores needed for mclapply

set.seed(copy*10)



# Create a folder "Tests" and a subfolder "Tests/Data" for the simulated data sets
PostD = 4; K.control = 15; K.mut = 5; J = 5; n = 300; m = 5000
createtestdata(copy, K0 = K.control, newK = K.mut, d = J, n = n, sitect = m, tests = max(PostD), 
               postD = PostD )

load(paste("Tests/Data/TestData_",copy,".Rdata",sep=""))
Pairs = c( "t1t2", "t1t3", "t2t3", "t1t3_D", "t2t3_D" )
a=1/(nrow(rawdata)^2)
Jname = paste(shortname, copy, "_a",1/a,sep="") #e.g. Jname = Test1_25
name = paste( "Tests/", Jname, sep = "" )
dir.create(name)

source("SHJ_algor.R")






