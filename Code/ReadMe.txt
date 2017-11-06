
Main script: 

SHJTest.R illustrates how to create a test file and run through the algorithm. To run the code, go to the directory where all the supporting scripts are located, then enter 
bsub -o ~/LSFouts/Test.out R CMD BATCH --vanilla --args --shortname=Test --core=2 --copy=1 SHJTest.R test1.out



Supporting scripts:

Allfcns.R provides all the necessary functions.

library_cmdline.R provides all the source code needed for LSF command line input.

source.R enables automatic hierarchical divisive tree.

SHJ_algor.R, Tree.R, GibbsMC.R, FixGibbsMC.R are the scripts for the full algorithm. The default is designed for a week queue.

