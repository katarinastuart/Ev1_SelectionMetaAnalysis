# 	 This file is used to plot figures for the software Bayescan in R.

#    This program, BayeScan, aims at detecting genetics markers under selection,
#	 based on allele frequency differences between population. 
#    Copyright (C) 2010  Matthieu Foll
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Arguments:
# - file is the name of your file ex: "output_fst.txt"
# - there are three ways to choose a threshold value to detect loci under selection
#   a. directly using a Posterior Odds threshold (PO)
#   b. controling for a maximum False Discovery Rate (in the list of outliers, FDR is 
#      the expected proportion of false positives ,i.e. misclassified). PO is 
#      estimated from the  maximum FDR you choose (the estimated FDR can be lower)
#   c. choosing a weight w (0-1), the function will find the PO which minimizes 
#   S=w*FNDR+(1-w)*FDR. In this case a second plot is drawn showing the value of the
#	score S depending of the PO threshold chosen.
#   One of PO or FDR parameter must be present (not both)
# - n is the number of points used to find the optimal PO using case c. above
# - size is the size of the points and text labels for outliers
# - pos is the distance between the points and the labels 

# Output:
# This function returns different paremeters in a list
# - PO: posterior odds used as a threshold
# - FDR: the estimated False Discovery Rate 
# - FNDR: the estimated False Non Discovery Rate (in the list of non-outliers, 
#   FNDR is the expected proportion of false negative, i.e. misclassified)
# - p: the posterior probability used as a threshold (PO=p/(1-p))
# - outliers: the list of outliers
# - nb_outliers: the number of outliers

# Typical usage: 
# - load this file into R (file/source R code)
# - in R, go to the directory where "output_fst.txt" is (file/change current dir)
# - at the R prompt, type 
# > plot_bayescan("output_fst.txt",PO=10) 
# or
# > plot_bayescan("output_fst.txt",FDR=0.05) (certainly the prefered option)
# or
# > plot_bayescan("output_fst.txt",w=0.9) 
# if you save the output in a variable, you can recall the different results:
# results<-plot_bayescan("output_fst.txt",FDR=0.05)
# results$PO
# results$FDR
# results$FNDR
# results$p
# results$outliers
# results$nb_outliers


#
# plotting posterior distribution is very easy in R with the output of BayeScan:
# first load the output file *.sel produced by BayeScan
# > mydata=read.table("bi.sel",colClasses="numeric")
# choose the parameter you want to plot by setting for example:
# parameter="Fst1"
# then this line will make the plot for:
# > plot(density(mydata[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))
# you can plot population specific Fst coefficient by setting
# parameter="Fst1"
# if you have non-codominant data you can plot posterior for Fis coefficients in each population:
# parameter="Fis1"
# if you test for selection, you can plot the posterior for alpha coefficient for selection:
# parameter="alpha1"
# you also have access to the likelihood with:
# parameter="logL"
# if you have the package "boa" installed, you can very easily obtain Highest Probability 
# Density Interval (HPDI) for your parameter of interest (example for the 95% interval):
# > boa.hpd(mydata[[parameter]],0.05)


plot_bayescan<-function(file,PO=-1,FDR=-1,w=-1,n=1000,size=1,pos=0.35)
{
res=read.table(file)

PO2FDR<-function(res,PO)
{
q=PO/(1+PO)
significant=res[res$prob>=q,]
FDR=sum(1-significant$prob)/nrow(significant)
FDR
}

FDR2PO<-function(res,FDR)
{
q_vals=sort(res$prob)
for (q in q_vals) {
  significant=res[res$prob>=q,]
  if(sum(1-significant$prob)/nrow(significant)<=FDR) {
      break
      }
  } 
q/(1-q) 
}


BE<-function(res,PO,w)
{
q=PO/(1+PO)
non_significant=res[res$prob<q,]
FNDR=sum(non_significant$prob)/nrow(non_significant)

significant=res[res$prob>=q,]
FDR=sum(1-significant$prob)/nrow(significant)

(1-w)*FDR+w*FNDR
}

if (FDR==-1 && PO>0 && w==-1) {
  FDR=PO2FDR(res,PO)  
}
else if (PO==-1 && FDR>0 && w==-1) {
  PO=FDR2PO(res,FDR)
  FDR=PO2FDR(res,PO)
}
else if (PO==-1 && FDR==-1 && w>=0) {
  logPO=seq(0,4,length=n)
  errors=0
  for (i in 1:n)
  {
    errors[i]=BE(res,10^(logPO[i]),w)
  }
  PO=10^(logPO[max(which(errors==min(errors)))])
  FDR=PO2FDR(res,PO)
}
else {  
  stop("One of PO, FDR or w is required as an argument of the function")
}
if (PO==Inf)
	p=1
else
	p=PO/(1+PO)
non_significant=res[res$prob<p,]
FNDR=sum(non_significant$prob)/nrow(non_significant)

# when p=1 replaces by p=0.9999 to plot with logit scale
if(nrow(res[res$prob>=0.9999,])>=1)
res[res$prob>=0.9999,]$prob=0.9999
if(nrow(res[res$prob<=0.001,])>=1)
res[res$prob<=0.001,]$prob=0.001

outliers=as.integer(row.names(res[res$prob>=p,]))

ok_outliers=TRUE
if (PO==Inf | sum(res$prob>=p)==0)
	ok_outliers=FALSE;

if (w>=0) {
  #par(mfrow=c(2,1))
  par(mar = c(5, 4, 4, 4) + 0.3)
}
# plot
plot(log10(res$prob/(1-res$prob)),res$fst,xlab="log(PO)",ylab="Fst",pch=19,cex=size)
# add names of loci over p and vertical line
if (ok_outliers) {
	text(log10(res[res$prob>=p,]$prob/(1-res[res$prob>=p,]$prob))+pos*(round(runif(nrow(res[res$prob>=p,]),1,2))*2-3),res[res$prob>=p,]$fst,row.names(res[res$prob>=p,]),cex=size)
	lines(c(log10(PO),log10(PO)),c(-1,1),lwd=2)
}

if (w>=0) {
  par(new = TRUE)
  plot(logPO,errors,xlim=range(log10(res$prob/(1-res$prob))), type = "l", axes = FALSE, lty=3, bty = "n", xlab = "", ylab = "")
#  plot(logPO,errors,type="l",xlab="log(PO)",ylab="Posterior Odds")
  axis(side=4, at = pretty(range(errors)))
  mtext("Posterior Odds", side=4, line=3)
}
#par(mfrow=c(1,1))
return(list("PO"=PO,"FDR"=FDR,"FNDR"=FNDR,"p"=p,"outliers"=outliers,"nb_outliers"=length(outliers)))
}



