#######################
#   John Henderson    #
#    UC-Berkeley      #            
#       &             #
#   Sara Chatfield    #
#    UC-Berkeley      #       
#                     # 
#    Who Matches:     #
#   Replication Code  #
#                     #
#   Released          #
#     Version 1.0     #
#                     #
#   Nov 17, 2010      #
#######################

Included are ten files:

- ReadMe.txt
- data_all_waves.dta
- WhoMatches.Rdata
- indicator.Rdata
- funcs.R
- objects.R
- indicator.R
- GenMatchCode.R
- ReplicationCode.R
- itersPscoreExample.R

The files GenMatchCode.R and ReplicationCode.R contain the replication code for the main portion of the paper. Each of these files sources the objects.R and funcs.R files, and explains how the code recovers the findings in the paper.  

The objects.R file builds the relevant data objects from the Stata file, data_all_waves.dta, and the funcs.R loads the particular functions that are used in addition to the 'MASS', 'Matching', 'foreign' packages which are also required. 

(The file itersPscoreExample.R contains code that illustrates the pscore iteration analysis.  The main functionality here is contained in the funcs. R file however.  Because of their large size, the files containing the propensity score simulation data, e.g. the data in Figure 2, are included separately.)


END