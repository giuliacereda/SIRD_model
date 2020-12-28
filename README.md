# SI2R2Dmodel
Analysis of COVID spread in Tuscany and Italy


 The repository contains the following files:
 
 Run.R
 Functions.R
 Plot.R
 
 The repository contains the following folders:
 WS
 Figures
 
## Run.R
1.The main code is in the file Run.R It calls the functions of the file “Functions.R” 
Put all files in the same folder and set it as the working directory through the command setwd(pathname) in Run.R
2. Define the time window for the analysis appropriately setting the variables data.ini and data.exit
3. If you want to use your cores in parallel, set: vuoi.parallelizzare<-1, otherwise 0
4. If you want to modify the value of the infection fatality rate, change the value assigned to the variable p_seq (now p_seq=0.0114)
5. Set the number of boostrap, for example n.boot<-100
 

N.B. The code models R0(t) through a natural cubic regression spline with 4 internal knots placed  at the quantiles. One can increase the number of knots by adding new parameters to the R0(t) curve in the code (i.e. by adding new variables r6, r7…). This could be needed when extending the study period.  

N.N.B with  the code takes around tot hours with 16 virtual chores to run.

When the codes has ended running, a copy of the workingspace is created in the folder WS, the name of the file contains the final data 


## Plot.R

This code produces the plots that are saved in the folder Results
In order to produce your own plots, you have to set the following variables:

1. schermo=0 if you want the code to produce a pdf in the folder Result, schermo=1 if you want the figures to be plotted on screen
2. load the working space you want via the command load. Remember the working space are in the subfolder WS created by the code Run.R

