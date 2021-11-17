library("Rtsne")
require(VariantAnnotation)
require(reshape)
library(ggplot2)
library(gridExtra)
library(plotly)

preplotVariantCrossSamples=function(variant_vector=NA, samples,col,new_plot=TRUE){
  if (is.na(variant_vector)){
    if (new_plot){
      plot(correctAF[strictSomatic,samples][1,],type='l',col='white',ylim=c(0,1),ylab='AF')
    }
    for(i in 1:sum(strictSomatic)){
      lines(correctAF[strictSomatic,samples][i,], col=col)
    }
  }else{
    if (new_plot){
      plot(correctAF[strictSomatic,][variant_vector,samples][1,],type='l',col='white',ylim=c(0,1),ylab='AF')
    }
    for(i in 1:sum(variant_vector)){
      lines(correctAF[strictSomatic,][variant_vector,samples][i,], col=col)
    }
  }
}

patient=""

vcf=readVcf(paste(patient, ".somatic.nocnv.vcf", sep=""), "hg38")

AO = as.data.frame(geno(vcf)$AO); # Number of alt allele oberservations
RO = geno(vcf)$RO; # Number of ref allele overservations
DP = geno(vcf)$DP; # Read depth
SAF = as.vector(info(vcf)$SAF); # Alt obs on forward strand
SAR = as.vector(info(vcf)$SAR); # Alt obs on reverse strand

selector=complete.cases(DP);
selector = selector & SAF>1 & SAR>1 

AF=matrix(rep(0, NROW(AO)*NCOL(AO)), ncol=NCOL(AO), byrow=T);
for(i in 1:NROW(AO)) {
  for(j in 1:NCOL(AO)) {
    AF[i, j] = AO[[i,j]] / (AO[[i,j]] + RO[[i, j]]);
  }
}
colnames(AF)=colnames(AO)
rownames(AF)=rownames(AO)

for(i in 1:NROW(DP)){
  for(j in 1:NCOL(DP)){
    if (DP[[i,j]]<50){
      selector[[i]]=FALSE
    }
  }
}
topo.colors(2)
print("depth >=50")
print(sum(selector))

binary=matrix(rep((0), NROW(AO)*NCOL(AO)), ncol=ncol(AO), byrow=T)
binary[,1]=(AO[,1] >5 | AF[,1]>0.1)
strictSomatic=selector & !(binary[,1])
print("strictSomatic")
print(sum(strictSomatic))
correctAF=AF

preplotVariantCrossSamples(samples=c(1:length(AO)),col="#00000090")

strictAF=AF[strictSomatic,]
data = as.data.frame(strictAF)

a <- ggplot(data, aes(data$`17566X`,data$`17566X`, alpha=.5,label=rownames(data))) +
  geom_point() +
  xlab("Baseline") +
  ylab("Y0.5") +
  xlim(0,1) +
  ylim(0,1) +
  theme_classic() +
  theme(legend.position = "none")

b <- ggplot(data, aes(data$`17566X`,data$`17566X`, alpha=.5,label=rownames(data))) +
  geom_point() +
  xlab("Y0.5") +
  ylab("Y1") +
  xlim(0,1) +
  ylim(0,1) +
  theme_classic() +
  theme(legend.position = "none")

c <- ggplot(data, aes(data$`17566X`,data$`17566X`, alpha=.5,label=rownames(data))) +
  geom_point() +
  xlab("Y1") +
  ylab("Y2") +
  xlim(0,1) +
  ylim(0,1) +
  theme_classic() +
  theme(legend.position = "none")

d <- ggplot(data, aes(data$`17566X`,data$`17566X`, alpha=.5,label=rownames(data))) +
  geom_point() +
  xlab("Y2") +
  ylab("Y3") +
  xlim(0,1) +
  ylim(0,1) +
  theme_classic() +
  theme(legend.position = "none")

e <- ggplot(data, aes(data$`17566X`,data$`17566X`, alpha=.5,label=rownames(data))) +
  geom_point() +
  xlab("Y3") +
  ylab("Y3.5") +
  xlim(0,1) +
  ylim(0,1) +
  theme_classic() +
  theme(legend.position = "none")

########################################################################################################################################
##                                           Removal of any variant with AF > 0.009 at germline                                       ##
########################################################################################################################################

selector=complete.cases(DP);
selector = selector & SAF>1 & SAR>1 

AF=matrix(rep(0, NROW(AO)*NCOL(AO)), ncol=NCOL(AO), byrow=T);
for(i in 1:NROW(AO)) {
  for(j in 1:NCOL(AO)) {
    AF[i, j] = AO[[i,j]] / (AO[[i,j]] + RO[[i, j]]);
  }
}
colnames(AF)=colnames(AO)
rownames(AF)=rownames(AO)

for(i in 1:NROW(DP)){
  for(j in 1:NCOL(DP)){
    if (DP[[i,j]]<50){
      selector[[i]]=FALSE
    }
  }
}
topo.colors(2)
print("depth >=50")
print(sum(selector))
unfiltered = selector

for(i in 1:NROW(AF)){
  if (AF[[i,1]]>0.009){
    selector[[i]]=FALSE
  }
}
print("Remove Germline")
print(sum(selector))

binary=matrix(rep((0), NROW(AO)*NCOL(AO)), ncol=ncol(AO), byrow=T)
binary[,1]=(AO[,1] >5 | AF[,1]>0.1)
strictSomatic=selector & !(binary[,1])
print("strictSomatic")
print(sum(strictSomatic))
correctAF=AF

preplotVariantCrossSamples(samples=c(1:length(AO)),col="#00000090")

print("With Germline")
print(sum(unfiltered))
print("strictSomatic")
print(sum(strictSomatic))

#############################################################
## Add color to the line plot based on apparent subclones. ##
## Values will change for each patient.                    ##
#############################################################

plotVariantCrossSamples=function(variant_vector=NA, samples,col,new_plot=TRUE){
  if (new_plot){
    plot(correctAF[strictSomatic,samples][1,],type='l',col='white',ylim=c(0,1),ylab='AF')
  }
  for(i in 1:sum(strictSomatic)){
    if (correctAF[strictSomatic,samples][i,2]>0.3){
      lines(correctAF[strictSomatic,samples][i,], col="blue")
    } else if (correctAF[strictSomatic,samples][i,2]<0.3 && correctAF[strictSomatic,samples][i,2]>0.1 ){
      if(correctAF[strictSomatic,samples][i,4]<0.1){
        lines(correctAF[strictSomatic,samples][i,], col="green")
      } else{
        lines(correctAF[strictSomatic,samples][i,], col="red")
      }
    } else{
      lines(correctAF[strictSomatic,samples][i,], col="black")
    }}}

plotVariantCrossSamples(samples=c(1:length(AO)),col="#00000090")