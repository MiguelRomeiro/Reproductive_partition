# Description of this script

# This script aims to:

# Calculate heterozygosity measures for samples specimens
# Calculate land cover around sampled nests
# Generate processed data files as outputs (already available at ~/data/processed/)

## Developed by Miguel P. Pereira-Romeiro
# Authors: Miguel P. Pereira-Romeiro, Marianne Azevedo-Silva, Henrique Silva Florindo, Gustavo Maruyama Mori, Paulo Silva Oliveira, Anette

# Install required packages

packages <- c("raster", "rgeos", "piecewiseSEM", "ggplot2", "dplyr", "car", "DHARMa", "performance", "gtools", "corrplot", "effectsize")

# Loop through the vector and install each package
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

install.packages("raster")
install.packages("rgeos")
install.packages("piecewiseSEM")
install.packages("terra")
# installed.packages("ggplot2")
# install.packages("dplyr")
# install.packages("car")
# install.packages("DHARMa")
# install.packages("performance")
# install.packages("gtools")
# install.packages("corrplot")
# install.packages("effectsize")
# install.packages("renv")

#Load required packages

library(raster)
library(rgeos)
library(piecewiseSEM)
library(ggplot2)
library(dplyr)
library(car)
library(DHARMa)
library(performance)
library(gtools)
library(corrplot)
library(renv)
library(terra)

#### Loading Data ####

landuse <- raster("~/Pesquisa/Reproductive_partition/data/raw/land_1m.asc")
landuse1 <- rast("~/Pesquisa/Reproductive_partition/data/raw/land_1m.asc")

data_rufi <- read.csv("~/Pesquisa/Reproductive_partition/data/raw/data_rufi.csv",sep =",", head=T)
data_rufi[,c(3,4, 12:length(data_rufi))] <- as.numeric(unlist(data_rufi[,c(3,4,12:length(data_rufi))]))
str(data_rufi)

data_reng <- read.csv("~/Pesquisa/Reproductive_partition/data/raw/data_reng.csv",sep =",", head=T)
data_reng[,c(3,4, 12:length(data_reng))] <- as.numeric(unlist(data_reng[,c(3,4,12:length(data_reng))]))
str(data_reng)


#### Calculating Heterozigosity ####

## Define GENHET function (Coulon, 2009)
"GENHET"<-
  function(dat,estimfreq,locname,alfuser){
    
    nbloc=(ncol(dat)-1)/2
    nbind=nrow(dat)
    
    #estimation of allele frequencies (only if estimfreq=T)
    
    if(estimfreq=="T")
      
    {
      
      #creation of the list of alleles
      datv=vector(length=nbind*nbloc*2)
      for (i in 2:ncol(dat)) datv[(nrow(dat)*(i-2)+1):(nrow(dat)*(i-1))]=dat[,i]
      al=sort(na.omit(unique(datv)))
      
      #count of the number of times each allele appears + nb of missing data
      alcount=matrix(nrow=(length(al)+1),ncol=(nbloc+1))
      alcount[,1]=c(al,NA)
      for(j in 1:(nrow(alcount)-1))
        for(k in 1:(ncol(alcount)-1))
          alcount[j,(k+1)]=sum(dat[,(k*2):(k*2+1)]==alcount[j,1],na.rm=T)
      for(l in 2:ncol(alcount))
        alcount[nrow(alcount),l]=(2*nbind-sum(alcount[1:(nrow(alcount)-1),l]))
      
      
      #creation of the table of allele frequencies
      alfreq=matrix(nrow=length(al),ncol=(nbloc+1))
      colnames(alfreq)=c("Allele",locname)
      alfreq[,1]=al
      for(m in (1:nrow(alfreq)))
        for (n in 2:ncol(alfreq)) alfreq[m,n]=alcount[m,n]/(nbind*2-alcount[nrow(alcount),n])
      
    }
    
    else alfreq=alfuser
    
    
    dat=as.data.frame(dat)
    library(gtools)
    res=matrix(nrow=nrow(dat),ncol=6)
    colnames(res)=c("sampleid","PHt","Hs_obs","Hs_exp","IR","HL")
    res[,1]=as.character(dat[,1])
    
    
    
    #estimation of E per locus (for HL and Hs_exp)
    E=vector(length=nbloc)
    alfreq2=alfreq[,2:ncol(alfreq)]*alfreq[,2:ncol(alfreq)]
    for(k in 1:ncol(alfreq2)) E[k]=1-sum(alfreq2[,k],na.rm=T)
    
    
    #estimation of the mean heterozygosity per locus
    mHtl=vector(length=nbloc)
    ctNAl=0
    ctHtl=0
    for(l in 1:ncol(dat))
    {if (even(l)==T)
    {
      for (m in 1:nrow(dat))
      { if (is.na(dat[m,l])==T) ctNAl=(ctNAl+1)
      else if (is.na(dat[m,(l+1)])==T) ctNAl=(ctNAl+1)
      else if (dat[m,l]!=dat[m,(l+1)]) ctHtl=(ctHtl+1)
      }
      mHtl[l/2]=ctHtl/(nrow(dat)-ctNAl)
      ctNAl=0
      ctHtl=0
    }
    }
    
    #the program in itself
    
    ctHt=0
    ctNA=0
    
    ctHm=0
    smHtl=0
    mmHtl=0
    sE=0
    mE=0
    
    sfl=0
    
    sEh=0
    sEj=0
    
    for(i in 1:nrow(dat))
    { for (j in 2:(nbloc*2))
    { if (even(j)==T)
    {
      if (is.na(dat[i,j])==T) ctNA=(ctNA+1)
      else if (is.na(dat[i,(j+1)])==T) ctNA=(ctNA+1)
      else {
        if (dat[i,j]!=dat[i,(j+1)])
        {
          ctHt=(ctHt+1)
          sEj=sEj+E[j/2]
        }
        else sEh=sEh+E[j/2]
        smHtl=smHtl+mHtl[j/2]
        sE=sE+E[j/2]
        sfl=sfl+alfreq[alfreq[,1]==as.numeric(dat[i,j]),(j/2+1)]+alfreq[alfreq[,1]==as.numeric(dat[i,j+1]),(j/2+1)]
      }
    }
    }
      res[i,2]=ctHt/(nbloc-ctNA)
      mmHtl=smHtl/(nbloc-ctNA)
      res[i,3]=(ctHt/(nbloc-ctNA))/mmHtl
      mE=sE/(nbloc-ctNA)
      res[i,4]=(ctHt/(nbloc-ctNA))/mE
      ctHm=nbloc-ctHt-ctNA
      res[i,5]=(2*ctHm-sfl)/(2*(nbloc-ctNA)-sfl)
      res[i,6]=sEh/(sEh+sEj)
      ctHt=0
      ctNA=0
      ctHm=0
      smHtl=0
      mmHtl=0
      sE=0
      mE=0
      sfl=0
      sEh=0
      sEj=0
    }
    return(res)
  }

## Preparing data ##

locname_ru <- colnames(data_rufi[,seq(15,47, by=2)])
locname_ru

locname_re <- colnames(data_reng[,seq(15,47, by=2)])
locname_re

## Calculating heterozigosity ##

het_rufi <- GENHET(dat = data_rufi[,c(2,15:48)], estimfreq = "T", locname = locname_ru)
het_rufi <- as.data.frame(het_rufi)
het_rufi$sampleid <- as.factor(het_rufi$sampleid)

data_rufi <- merge(data_rufi, het_rufi, by = "sampleid")
data_rufi[,c(49:53)] <- as.numeric(unlist(data_rufi[,c(49:53)]))


het_reng <- GENHET(dat = data_reng[,c(2,15:48)], estimfreq = "T", locname = locname_re)
het_reng <- as.data.frame(het_reng)
het_reng$sampleid <- as.factor(het_reng$sampleid)

data_reng <- merge(data_reng, het_reng, by = "sampleid")
data_reng[,c(49:53)] <- as.numeric(unlist(data_reng[,c(49:53)]))


#### Calculating Buffers ####

coord_rufi <- SpatialPoints(unique(data_rufi[,c(4,3)])) #longitude has to come before latitude
coord_reng <- SpatialPoints(unique(data_reng[,c(4,3)])) #longitude has to come before latitude

# plot(landuse)
# points(coord_rufi, col = as.factor(data_reng$nq), pch = 16)
# points(coord_reng, col = as.factor(data_reng$nq), pch = 17)
# table(values(landuse))

## C. rufipes ##

result_ru=matrix(NA, nrow = length(unique(data_rufi$nest)), ncol = 7)
colnames(result_ru)= c(1,2,3,4,5,6,7)
rownames(result_ru)<-unique(data_rufi$nest)
head(result_ru)

names_ru <- unique(data_rufi$nest)

for(i in 1:nrow(coord_rufi@coords)){
  buffer_rufipes <- buffer(coord_rufi[i], width=4) 
  crop<-crop(landuse,extent(buffer_rufipes)) 
  mask<-mask(crop,buffer_rufipes) 
  plot(mask)
  y<-table(values(mask)) 
  y<-t(as.matrix(y)) 
  x<-matrix(0,1,7) 
  colnames(x)<-c("1","2","3","4","5","6","7") 
  z<-merge(x,y,all=T)
  z<-z[-1,]
  z<-z[,c("1","2","3","4","5","6","7")]
  z<-as.matrix(z)
  row = nrow(z)
  if (row==0){
    result_ru[(i),]<-NA
  } else {
    result_ru[(i),]<-z
  }
  rownames(result_ru)[i]<-names_ru[i]
}

result_ru

buffer_ru <- as.data.frame(result_ru)/rowSums(result_ru,na.rm=T)
colnames(buffer_ru)= c("agr","cer","cers","cult","wtr","urb","road")
buffer_ru$'cert' <- rowSums(buffer_ru[,c(2,3)],na.rm = T)
buffer_ru[is.na(buffer_ru)] <- 0
buffer_ru["nest"]<- row.names(buffer_ru)
buffer_ru<-buffer_ru %>% relocate(nest)
buffer_ru$nest<-as.factor(buffer_ru$nest)

buffer_ru

data_rufi<-merge(data_rufi, buffer_ru, by = "nest")

write.csv(data_rufi, file = "~/Pesquisa/Reproductive_partition/data/processed/rufi_processed.csv", row.names = FALSE)

## C. rengerii ##

result_re=matrix(NA, nrow = length(unique(data_reng$nest)), ncol = 7)
colnames(result_re)= c(1,2,3,4,5,6,7)
rownames(result_re)<-unique(data_reng$nest)
head(result_re)

names_re <- unique(data_reng$nest)

for(i in 1:nrow(coord_reng@coords)){
  buffer_renggeri <- buffer(coord_reng[i], width=4) 
  crop<-crop(landuse,extent(buffer_renggeri))
  mask<-mask(crop,buffer_renggeri) 
  plot(mask)
  y<-table(values(mask))
  y<-t(as.matrix(y)) 
  x<-matrix(0,1,7)
  colnames(x)<-c("1","2","3","4","5","6","7")
  z<-merge(x,y,all=T)
  z<-z[-1,]
  z<-z[,c("1","2","3","4","5","6","7")]
  z<-as.matrix(z)
  row = nrow(z)
  if (row==0){
    result_re[(i),]<-NA
  } else {
    result_re[(i),]<-z
  }
  rownames(result_re)[i]<-names_re[i]
}

result_re

buffer_re <- as.data.frame(result_re)/rowSums(result_re,na.rm=T)
colnames(buffer_re)= c("agr","cer","cers","cult","wtr","urb","road")
buffer_re$'cert' <- rowSums(buffer_re[,c(2,3)],na.rm = T)
buffer_re[is.na(buffer_re)] <- 0
buffer_re["nest"]<- row.names(buffer_re)
buffer_re<-buffer_re %>% relocate(nest)
buffer_re$nest<-as.factor(buffer_re$nest)

buffer_re

data_reng<-merge(data_reng, buffer_re, by = "nest")

write.csv(data_reng, file = "~/Pesquisa/Reproductive_partition/data/processed/reng_processed.csv", row.names = FALSE)

## Correlation of Variables ##

cor_rufi <- cor(data_rufi[,c(6,8,14,49,61)])
test_rufi = cor.mtest(data_rufi[,c(6,8,14,49,61)], conf.level = 0.95)
corrplot(cor_rufi, p.mat = test_rufi$p, method = 'color', type = 'lower', addCoef.col ='black', number.cex = 1.5, cl.cex = 1, order = 'alphabet', diag=FALSE, main="Camponotus rufipes")

cor_reng <- cor(data_reng[,c(6,8,14,49,61)])
test_reng = cor.mtest(data_reng[,c(6,8,14,49,61)], conf.level = 0.95)
corrplot(cor_reng, p.mat = test_reng$p, method = 'color', type = 'lower', addCoef.col ='black', number.cex = 1.5, cl.cex = 1, order = 'alphabet', diag=FALSE, main="Camponotus renggeri")

#CRIAR CSV --> DADOS PROCESSADOS