#general functions
packages=c("preprocessCore","pROC","pheatmap","ggplot2","gridExtra","grid","ggpubr","scales","stringr","qvalue","vegan","stringr","tidyverse","pROC","survival","readr","gdata","DescTools","moments","fastcluster","Rfast","dynamicTreeCut","data.table","plyr","dplyr")
sapply(packages,require,character=TRUE)
#####################################################################################
###General Modules
#####################################################################################
#scale the max of a vector to 1 and min to 0
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
#filter by one list
filter <- function (x,filter,inlist) {
  if (inlist){
    x <- x[x %in% filter]  
  }
  else{
    x <- x[!x %in% filter]  
  }
  return (x)
}
#filter by two list
filter.double <- function (x,include,exclude) {
  x <- x[(x %in% include) & (!x %in% exclude)]  
  return (x)
}
#turn genotype list into a matrix fast
genolist2matrix<-function(genolist){
un.genolist <- unlist(genolist)
cell<-sort(unique(un.genolist))
tmp.bin <- matrix(0, nrow = length(genolist), ncol = length(cell))
dimnames(tmp.bin) <- list(names(genolist), cell)
ij <- cbind(rep(1:length(genolist), lengths(genolist)), match(un.genolist,cell))
tmp.bin[ij] <- 1
return(tmp.bin)
}
#turn matrix into genolist file
matrix2genofile<-function(matrix,outfile){
  for (i in 1:nrow(matrix)){
    tmp.out=c(unlist(strsplit(rownames(matrix)[i], ":(?=[^:]+$)", perl=TRUE)),names(which(matrix[i,]==1)))#corrected from: strsplit(rownames(matrix)[i],":")
    write(tmp.out,file=outfile,append=TRUE,sep="\t",ncolumns = length(tmp.out))
  }
}
#take nth root
nthroot<-function(x,n){
  abs(x)^(1/n)*sign(x)
}
#read gmt files into a list
read_gmt <- function(fileName, min = 5) {
  Lines <- fread(fileName, sep = "")[[1]] #read lines. This seems to remove the annotation line start with #
  read.list <- lapply (Lines, function(x) {
    genotype=strsplit(x, "\t")[[1]]
    gname=paste(genotype[1],genotype[2],sep=":")
    return(list(gname,genotype[-c(1:2)]))
    })
  genotype.list=lapply(read.list, `[[`, 2) 
  names(genotype.list)= lapply(read.list, `[[`, 1)
  genotype.list=genotype.list[lengths(genotype.list)>=min]
  return(genotype.list)
}
#write a list to gmt file
write_gmt <- function(genotype.list, file.name) {
  if (file.exists(file.name)) {
    message(paste(file.name, "exists, deleting..."))
    file.remove(file.name)
  }
  write("#GMT file written from a genotypelist", file=file.name, sep="\t",append=TRUE, ncolumns=1)
  for (i in seq_along(genotype.list)) {
    genotype <- genotype.list[[i]]
    gname=names(genotype.list)[i]
    ncolumns <- 2 + length(genotype)
    write(c(unlist(strsplit(gname, ":(?=[^:]+$)", perl=TRUE)), genotype), file=file.name, sep="\t",append=TRUE, ncolumns=ncolumns)
  }
}
#read preCal gmt file into a list object
read_preCal <- function(fileName) {
  Lines <- fread(fileName, sep = "")[[1]] #read lines. This seems to remove the annotation line start with #
  preCal.list=list()
  for (Line in Lines){
    tmp.line<-unlist(strsplit(Line,"\t"))
    tmp.split=strsplit(tmp.line[3:length(tmp.line)],"@")
    feature.redun=as.numeric(sapply(tmp.split, `[[`, 2))
    names(feature.redun)= sapply(tmp.split, `[[`, 1)
    preCal.list[[tmp.line[1]]]=feature.redun
  }
return(preCal.list)
}
#find concepts that contains the cosmicID in a conceptList
findConcepts<-function(cosmicID,conceptList){
  concepts<-c()
  for(i in 1:length(conceptList)){
    if(cosmicID %in% conceptList[[i]]){
      concepts<-append(concepts,names(conceptList)[i])
    }
  }
  return(concepts)
}
#####################################################################################
###General Modules for GDSC and CCLE cell line datasets
#####################################################################################
#define sensitive and resistant cell lines based on AUC
define.response.AUC<-function(AUC){ #AUC should be a vector of AUC
  iAUC.sort<-sort(1-AUC)
  mySlope=max(iAUC.sort)/which.max(iAUC.sort)
  cutoff.sen<-1-iAUC.sort[which.max((mySlope*1:length(iAUC.sort)-iAUC.sort)/sqrt(1+mySlope^2))]
  cutoff.res<-median(AUC)+mad(AUC,constant=1)
  return(c(cutoff.sen,cutoff.res))
}
#define sensitive and resistant cell lines based on ActArea
define.response.ActArea<-function(ActArea){#ActArea should be a vector of ActArea
  iActArea.sort<-sort(ActArea)
  mySlope=max(iActArea.sort)/which.max(iActArea.sort)
  cutoff.sen<-iActArea.sort[which.max((mySlope*1:length(iActArea.sort)-iActArea.sort)/sqrt(1+mySlope^2))]
  cutoff.res<-median(ActArea)-mad(ActArea,constant=1)
  return(c(cutoff.sen,cutoff.res))
}
#generate permutated GDSC testing sets based on n folds
# select 20 testing cell line lists for all drugs. Testing cell lines: 1/10, training cell lines: 9/10
permuTestingSetV2<-function(drug.selection,drugData,resultdir,chartdir,testperc=0.2,ks.p=0.1,min.sensitive=5){
  #drugID.vec<-drug.selection$drugID[drug.selection$SensitiveCellLine.n.. >= 20 & drug.selection$Skewness < 0]
  drugID.vec<-drug.selection$drugID
  for(drugID.call in drugID.vec){
    subData.call<-drugData[drugData$DRUG_ID==drugID.call,]
    density.all<-density(subData.call$AUC)
    pdf(file=paste(chartdir,"/density_allvspermuTestingSet_drugID",drugID.call,".pdf",sep=""))
    plot(density.all,col="black",xlab="AUC",ylab="Density",main=paste("Density plots of drugID ",drugID.call,sep=""))
    cutoff=drug.selection$Sensitive.cutoff[drug.selection$drugID==drugID.call]
    senData.call=subData.call[subData.call$AUC<cutoff,]
    restData.call=subData.call[subData.call$AUC>=cutoff,]
    senCount=length(senData.call$COSMIC_ID)
    restCount=length(restData.call$COSMIC_ID)
    i=1
    skip=0
    while(i <=20 & skip<=100){
      senTest<-sample(1:senCount,senCount*testperc)
      restTest<-sample(1:restCount,restCount*testperc)
      senTrain<-which(!1:senCount %in% senTest)
      restTrain<-which(!1:restCount %in% restTest)
      if (length(senTest)<3){
        skip=skip+1
        next
      }
      suppressWarnings(senKS<-ks.test(senData.call$AUC[senTrain],senData.call$AUC[senTest])) #there will be a warning if two or more cell lines have same auc value
      suppressWarnings(restKS<-ks.test(restData.call$AUC[restTrain],restData.call$AUC[restTest]))
      if(senKS$p.value<=ks.p|restKS$p.value<=ks.p){
        print("Skip one permutation attempt because of significant deviation from original AUC distribution")
        skip=skip+1
        next
      }
      subData.test <- rbind(senData.call[senTest,],restData.call[restTest,])
      density.sub <- density(subData.test$AUC)
      lines(density.sub$y~density.sub$x,type="l",col="grey")
      test.senstive.n=length(subData.test$COSMIC_ID[subData.test$COSMIC_ID %in% senData.call$COSMIC_ID])
      if(test.senstive.n<min.sensitive){
        print(paste0("Skip one permutation attempt because sensitive subject n=",test.senstive.n,",less than ",min.sensitive))
        skip=skip+1
        next        
      }
      write.table(subData.test[,c(1,3,7)],file=paste(resultdir,"/testingSet_drugID",drugID.call,"_permu",i,".txt",sep=""),
                  row.names=FALSE,col.names=T,sep="\t")
      i<-i+1
    }
    lines(density.all$y~density.all$x,col="black")
    print(paste("20 permutated testing lists generated for drugID:#",drugID.call,sep=""))
    dev.off()
  }
}
#####################################################################################
###Freeze version (Trim Corrected): Generate 6 Overlapping Expression Bins Generated 
###Based On Cutoffs created using trimmed mean and sd of fold change using unlog TPM 
###expression data
#####################################################################################
#sdlevels means the levels of standard deviations must be symmetric and must include zero (the default is -3 to 3)
binarize.expfeature<-function(expData,outfilepath,sdpercfilter=0.2,sdlevels=-6:6,genomewide=TRUE){
  suppressWarnings(file.remove(outfilepath))
  if (length(grep("_1$", colnames(expData)))>0){
    expData<-expData[, -grep("_1$", colnames(expData))]#remove duplcate columns (Duplicated column names deduplicated: '1503362' => '1503362_1' [845])  
  }
  rownames(expData)<-gsub("\\.\\d+$", "", rownames(expData)) #remove .\d+ in ensemble_gene
  expMatrix<-as.matrix(expData)
  if (genomewide==TRUE){
    expMatrix=normalize.quantiles(expMatrix)
    dimnames(expMatrix) <- list(rownames(expData), colnames(expData))
    exp.sd <- apply(expMatrix, 1, sd)
    expMatrix<-expMatrix[which(exp.sd > quantile(exp.sd, sdpercfilter)), ] 
  }
  samples=colnames(expMatrix)
  writeLines(paste(c(paste("#Overlapping Expression Bins Generated Based On Cutoffs created using Trimmed (0.1) Mean + ",paste(sdlevels[sdlevels!=0],collapse="/"),"*SD based on log2 transformed fold change for quantile normalized TPM",sep=""),samples),collapse='\t'),outfilepath)
  for (i in 1:nrow(expMatrix)){
    expi=as.numeric(expMatrix[i,])
    names(expi)=names(expMatrix[i,])
    #expi=expi+min(expi[expi!=0])/2
    expi=expi+1
    expi=log2(expi/mean(expi,trim=0.1))
    expi.trim=Trim(expi,trim=0.1)
    if (sd(expi.trim)==0){
      next
    }
    cutoffs=mean(expi.trim)+sd(expi.trim)*sort(sdlevels[sdlevels!=0])
    expbins<-cut(expi,breaks=c(-Inf,cutoffs,Inf), labels=sort(sdlevels))
    names(expbins)=names(expi)
    binlevels=as.numeric(levels(expbins))
    for (bin in binlevels){
        if (bin>0){
          tmp.out=c(rownames(expMatrix)[i],paste("UP","_Level",abs(bin),sep=""),names(expbins[which(expbins %in% binlevels[binlevels>=bin])]))
        }else if(bin<0){
          tmp.out=c(rownames(expMatrix)[i],paste("DN","_Level",abs(bin),sep=""),names(expbins[which(expbins %in% binlevels[binlevels<=bin])]))
        }else{
          tmp.out=c()
        }
        if (length(tmp.out)>6){
          write(tmp.out,file=outfilepath,append=TRUE,sep="\t",ncolumns = length(tmp.out))  
      }
    }
  } 
}
#################################################################################   
###########Compute Feature Redundancy New Version v20210715 Clustering based#####    
#################################################################################   
#fast calculate TCGA jaccard index matrix for genotypes of a subject and generate precalculation files
calSimilarity_subject_tcga<-function(subject.call,tcga.genotype.list,subject.genotype.list,method="ochiai",cutoff=0.1,by.cluster=TRUE,cutTree.method=c("quantile","dynamic"),dynamic.method="tree",minClusterSize,split.depth=2,visualize.cluster=FALSE,log.dir){ #method=="ochiai" or "jaccard"
  match.arg(cutTree.method)
  tmp.genoCell<-findConcepts(subject.call,subject.genotype.list)
  tmp.genoCell<-tmp.genoCell[tmp.genoCell %in% names(tcga.genotype.list)]
  tmp.tcga.genoList<-tcga.genotype.list[tmp.genoCell]
  if (length(tmp.tcga.genoList)==0){
    tmp.out<-c(0)
  }else{
    tmp.genolist.bin<-genolist2matrix(tmp.tcga.genoList)
    if (method=="ochiai"){
      dist.matrix<-as.matrix(designdist(tmp.genolist.bin, method = "1-J/sqrt(A*B)"))
    }else if (method=="phi"){
      dist.matrix<-as.matrix(designdist(tmp.genolist.bin, method = "((A-J)*(B-J)-J*(P-A-B+J))/sqrt(A*B*(P-A)*(P-B))"))
    }else if(method=="jaccard"){
      dist.matrix<-as.matrix(designdist(tmp.genolist.bin, method = "(A+B-2*J)/(A+B-J)"))
    }else if(method=="raup-crick"){
      dist.matrix<-as.matrix(designdist(tmp.genolist.bin, method = "1-phyper(J-1, A, P-A, B)"))
    }else if (method=="kulczynski"){
      dist.matrix<-as.matrix(designdist(tmp.genolist.bin, method = "1-(J/2)*(1/A+1/B)"))
    }else{
      print ("method must be set as one of the following (pleae match case): ochiai,phi,jaccard,raup-crick, or kulczynski")
    }
    if (method=="phi"){
      sim.matrix=-dist.matrix  
    }else{
      sim.matrix=1-dist.matrix
    }
    sim.matrix <- ifelse(sim.matrix<cutoff,0,sim.matrix)
    if (by.cluster==TRUE & nrow(dist.matrix)>20){
      cluster.feature=fastcluster::hclust(as.dist(dist.matrix), method = "ward.D2")
      #groups=cutree(cluster.feature,h=Rfast::nth(unique(cluster.feature$height), 5, descending = T))
      if (cutTree.method=="quantile"){
        quantile.tree <- quantile(unique(cluster.feature$height), probs = seq(0, 1, 0.0001))
        groups=cutree(cluster.feature,h=quantile.tree["99.95%"]) 
      }else if (cutTree.method=="dynamic"){
        log <- capture.output({
          groups=cutreeDynamic(cluster.feature,minClusterSize = minClusterSize,method = dynamic.method,distM = dist.matrix,deepSplit = split.depth)
          names(groups)=cluster.feature$labels
        }) 
      }
      write(c(paste(length(unique(groups)),"feature groups detected for subject",subject.call)),file=paste0(log.dir,"/preCal.cluster.number.log"),append=TRUE)
      if (visualize.cluster==TRUE){
        dist.matrix.subject<-as.matrix(designdist(t(tmp.genolist.bin), method = "1-J/sqrt(A*B)"))
        cluster.subject=fastcluster::hclust(as.dist(dist.matrix.subject), method = "ward.D2")
        jpeg(paste0(log.dir,"/",gsub(" ",".",gsub(":","-",Sys.time())),".SUBJECT",subject.call,".cluster.heatmap.jpg"), width = 5000, height =5000)
        pheatmap(tmp.genolist.bin, cluster_rows=cluster.feature, cluster_cols=cluster.subject, cutree_rows=length(unique(groups)), show_rownames = FALSE, show_colnames = FALSE)
        graphics.off()
        jpeg(paste0(log.dir,"/",gsub(" ",".",gsub(":","-",Sys.time())),".SUBJECT",subject.call,".cluster.cormap.jpg"), width = 5000, height =5000)
        quantile.range <- quantile(sim.matrix, probs = seq(0, 1, 0.01))
        myBreaks <- seq(quantile.range["10%"], quantile.range["90%"], 0.01)
        myColor  <- colorRampPalette(c("skyblue", "white", "red"))(length(myBreaks) - 1)
        pheatmap(sim.matrix, cluster_rows=cluster.feature, cluster_cols=cluster.feature, breaks=myBreaks,color=myColor,cutree_rows=length(unique(groups)), cutree_cols=length(unique(groups)), show_rownames = FALSE, show_colnames = FALSE)
        graphics.off()
      }
      sum.similarity=c()
      for(i in 1:ncol(sim.matrix)){
        featurei=colnames(sim.matrix)[i]
        groupi=groups[featurei]
        groupfeatures=names(groups)[groups==groupi]
        sumi=sum(sim.matrix[groupfeatures,i])
        if (sumi==0){
          stop(paste("subject:",subject.call,"col:",i,"feature:",featurei,"group:",groupi,"sumi=0"))
        }
        sum.similarity=setNames(c(sum.similarity, sumi), c(names(sum.similarity), featurei))
      }
    }else{
      sum.similarity<-colSums(sim.matrix)
    }
    tmp.out=sum.similarity
  }
  return(tmp.out)
}
#batch precalculate dataset genotype redundacy based on TCGA and genotype lists
batch_calSimilarity_genolist<-function(subject.genotype.list,REF.genotype.list,subject.preCalfile,minTCGA=20,method="ochiai",cutoff=0.1,by.cluster=TRUE,cutTree.method="dynamic",dynamic.method="hybrid",minClusterSize=40,split.depth=2,log.dir,visualize.cluter=FALSE){
  subject.id<-unique(unlist(subject.genotype.list))
  suppressWarnings(file.remove(subject.preCalfile))
  lostperc<-length(names(subject.genotype.list)[!names(subject.genotype.list) %in% names(REF.genotype.list)])/length(names(subject.genotype.list))
  header=paste0("#Precalculations of genotype redundancies: ","minTCGA=",minTCGA,", method=",method,", cutoff=",cutoff,", by.cluster=",by.cluster,", cutTree.method=",cutTree.method,", dynamic.method=",dynamic.method,", minClusterSize=",minClusterSize,", split.depth=",split.depth,". Loss of genotype from TCGA: ",100*lostperc,"%","\n")
  feature.redund=list()
  for (i in 1:length(subject.id)){
      if (i %% 100==0){
        print (paste("Processed",i,"subjects"))
        feature.redund[[subject.id[i]]]<-calSimilarity_subject_tcga(subject.call=subject.id[i],tcga.genotype.list=REF.genotype.list,subject.genotype.list=subject.genotype.list,method=method,cutoff=cutoff,by.cluster=by.cluster,cutTree.method=cutTree.method,dynamic.method=dynamic.method,minClusterSize=minClusterSize,split.depth=split.depth,visualize.cluster=visualize.cluter,log.dir = log.dir)
      }else{
        feature.redund[[subject.id[i]]]<-calSimilarity_subject_tcga(subject.call=subject.id[i],tcga.genotype.list=REF.genotype.list,subject.genotype.list=subject.genotype.list,method=method,cutoff=cutoff,by.cluster=by.cluster,cutTree.method=cutTree.method,dynamic.method=dynamic.method,minClusterSize=minClusterSize,split.depth=split.depth,visualize.cluster=F,log.dir = log.dir)
      }
    }
  result=list(feature.redundancy=feature.redund,feature.redundancy.parameters=header)
  save(result,file=subject.preCalfile)
}
#####################################################################################
### Modules for weighted random walk k-s test to calculate NES
#####################################################################################
# calculate ES based on weights (sorted in decreasing order)
cal_ES<-function(weight,posList){
  if (sum(weight<0)>0) {
    stop (paste("negative values found in weight, aborting calculation of ES",weight,sep=" "))#there cannot be nagative values in weight
  }
  weight<-sort(weight,decreasing=TRUE)
  onIndice<-which(names(weight) %in% posList)
  stepDown<-1/(length(weight)-length(onIndice))
  norWeight<-rep(-stepDown,length(weight))
  names(norWeight)<-names(weight)
  norWeight[onIndice]<-weight[onIndice]/sum(weight[onIndice]) #corrected: should not be abs(weight[onIndice]/sum(weight[onIndice])) and there cannot be negative values
  cumNorWeight<-cumsum(norWeight)
  minES<-min(cumNorWeight)
  maxES<-max(cumNorWeight)
  if (abs(minES)>=abs(maxES)){
    return(0)    
  }else{
    return(maxES)    
  }
}
#calculate permutated ES
perm_weightedKS<-function(weights,numOnList,nPermu=2000){
  permu=sapply(1:nPermu,function(x) cal_ES(weights,sample(names(weights),numOnList)))
  permu<-as.numeric(permu)
  return(permu)
}
#perform weighted KS tests
weightedKSV2<-function(weights,genotype,myPermu,p.cut=0.05,min.numOnList=5,correct.train=TRUE){ 
  cat("correct.train=", correct.train, "\n")
  mytest = lapply(genotype, function (x) {
    ES = cal_ES(weights, x)
    onIndice = which(names(weights) %in% x)
    numOnList = length(onIndice)
    if(numOnList < min.numOnList|ES==0) {
      return(c(0, 0, "", rep("",length(weights))))
    } else {
      permu = myPermu[[numOnList]]
      NES = ES/mean(permu)
      pValue = length(permu[which(permu >= ES)])/length(permu)
      NES.cell = rep("",length(weights))
      if (pValue<p.cut & correct.train==T){
        for (i in onIndice){
          genotype.i=x[x != names(weights[i])]
          ESi<-cal_ES(weights,genotype.i)
          NESi<-ESi/mean(myPermu[[numOnList-1]])
          NES.cell=replace(NES.cell, i, NESi)      
        } 
      }
      return(c(NES,pValue,"NA",NES.cell))
    }
  })
  mytest.out = bind_rows(mytest) %>% t()
  colnames(mytest.out) = c("NES","pValue","qValue",names(weights))
  mytest.out = cbind(Genotype = rownames(mytest.out), mytest.out)
  rownames(mytest.out) = NULL
  mytest.out = mytest.out[mytest.out[, 'qValue'] == "NA", ]
  mytest.out[, 'qValue']=qvalue(as.numeric(mytest.out[, 'pValue']), pi0 = 1)$qvalues
  mytest.out=mytest.out[as.numeric(mytest.out[, 'pValue'])<p.cut,]
  mytest.out = as.data.frame(mytest.out,stringsAsFactors=F)
  return(mytest.out)
}
#calculate weighted KS using weightedKSV2 based on AUC
weightedKS.AUCV2<-function(weight,genotype,preCal=TRUE,sensitive=TRUE,root=2,p.cut=0.05,permuN=2000,correct.outlier=FALSE,correct.train=TRUE){
   if(sensitive){
    print("Using 1-AUC as weight to rank sensitive cell lines")
    tmp.weight<-sort(1-weight) ##use 1-AUC to favor the sensitive cell lines
    if (correct.outlier==TRUE){
      mySlope=max(tmp.weight)/which.max(tmp.weight)
      cut=tmp.weight[which.max((mySlope*1:length(tmp.weight)-tmp.weight)/sqrt(1+mySlope^2))]
      tmp.weight=ifelse(tmp.weight>cut,nthroot(tmp.weight,root), tmp.weight) #weight transform nthroot = root of 1/n, here we used square root 
    }
  }else{
    print("Using AUC as weight to favor resistant cell lines")
    tmp.weight=weight
  }
  myPermu<-list()
  matchn=sapply(genotype,function(x){length(intersect(names(tmp.weight),x))})
  matchn=sort(unique(c(matchn,matchn-1)))
  k=0
  print(paste("Pre-calculating",permuN, "ES permutations",sep=" ")) #corrected: precalculate all possible positive numbers less than the max number
  for(j in matchn){
    myPermu[[j]]<-perm_weightedKS(tmp.weight,j,permuN)
    names(myPermu)[j]<-j
  }
  print("Performing weighted KS tests")
  ks.result<-weightedKSV2(weights=tmp.weight,genotype=genotype,myPermu=myPermu,p.cut=p.cut,min.numOnList=5,correct.train=correct.train)
  return(ks.result)
}
#run weighted KS file a permutated testing set file using weightedKSV2 module
## the "minsize" argument is to filter out genotypes that define a small number of cell lines. 
## For example, "minsize=5" will remove genotypes that define less than 5 cell lines.
run.weightedKS.fileV2<-function(test.file,drugData,genotype.list,outdir,minsize=5,correct.outlier=FALSE,correct.train=T){
  drugID.call<-as.numeric(unlist(str_extract_all(test.file, "(?<=drugID)\\d+"))[1])
  permuID<-as.numeric(unlist(str_extract_all(test.file, "(?<=permu)\\d+"))[1])
  subData.call<-drugData[drugData$DRUG_ID==drugID.call,]
  print(paste("Calculating Weighted KS for drugID:",drugID.call," permutation #",permuID,sep=""))
  #myTesting<-read.table(file=test.file,sep="\t",header=TRUE)
  myTesting<-read.delim(test.file,stringsAsFactors = F,check.names = F,header = T,sep="\t")
  subData.select<-subData.call[!subData.call$COSMIC_ID %in% myTesting$COSMIC_ID,]
  genotype.select<-lapply(genotype.list,filter.double,include=subData.call$COSMIC_ID,exclude=myTesting$COSMIC_ID)
  genotype.select<-genotype.select[which(sapply(genotype.select,function(x) length(x)>=minsize))]
  auc.weight<-subData.select[,"AUC"]
  names(auc.weight)<-as.character(subData.select[,"COSMIC_ID"])
  ks.sensitive<-weightedKS.AUCV2(weight=auc.weight,genotype=genotype.select,preCal = TRUE,sensitive = TRUE,correct.outlier=correct.outlier,correct.train=correct.train)
  ks.resistant<-weightedKS.AUCV2(weight=auc.weight,genotype=genotype.select,preCal = TRUE,sensitive = FALSE,correct.outlier=correct.outlier,correct.train=correct.train)            
  print("Outputing results")
  write.table(ks.sensitive,file=paste(outdir,"/weightTransf_NES_drugID_",drugID.call,"_permuID_",permuID,".sensitive.xls",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(ks.resistant,file=paste(outdir,"/weightTransf_NES_drugID_",drugID.call,"_permuID_",permuID,".resistant.xls",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
  print(paste("drugID_",drugID.call," finished",sep=""))
}
#run weighted KS for all cell lines of a drug ID using weightedKSV2 module
run.weightedKS.drugV2<-function(drugID.call,drugData,genotype.list,outdir,transform.root=2,p.cut=1,minsize=5,correct.outlier=FALSE,correct.train=T){
  subData.select<-drugData[drugData$DRUG_ID==drugID.call,]
  print(paste("Calculating Weighted KS for drugID:",drugID.call,sep=""))
  genotype.select<-lapply(genotype.list,filter,filter=subData.select$COSMIC_ID,inlist=TRUE) 
  genotype.select<-genotype.select[which(sapply(genotype.select,function(x) length(x)>=minsize))]
  auc.weight<-subData.select[,"AUC"]
  names(auc.weight)<-as.character(subData.select[,"COSMIC_ID"])
  ks.sensitive<-weightedKS.AUCV2(auc.weight,genotype.select,preCal = TRUE,sensitive = TRUE,root=transform.root,p.cut=p.cut,correct.outlier=correct.outlier,correct.train=correct.train)
  ks.resistant<-weightedKS.AUCV2(auc.weight,genotype.select,preCal = TRUE,sensitive = FALSE,root=transform.root,p.cut=p.cut,correct.outlier=correct.outlier,correct.train=correct.train)            
  print("Outputing results")
  write.table(ks.sensitive,file=paste(outdir,"/weightTransf_NES_drugID_",drugID.call,".sensitive.xls",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(ks.resistant,file=paste(outdir,"/weightTransf_NES_drugID_",drugID.call,".resistant.xls",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
  print(paste("drugID_",drugID.call," finished",sep=""))
}
#####################################################################################
### Modules for building iGenSig module for cell line datasets
#####################################################################################
#filter features and remove equivocal features based on feature weight table
#sig.feature.weight must contain Feature.List column
#feature.list.train must only contain the subject used in training sets
filter.feature<-function(feature.weight,feature.list.train,minsize,weight.cut,low.weight.cut,p.cut,low.p.cut,q.cut,low.q.cut,feature.col=1,highlevel.minsize=10,filter.criteria=c(1,2)){ #filter 1: Genotypes identified as low-level OE/UE but not high-level OE/UE, 2:Genotypes identified as both increased and decreased activity of the same gene
  feature.weight$selected=ifelse(as.numeric(feature.weight$NES)>weight.cut&as.numeric(feature.weight$pValue)<p.cut&as.numeric(feature.weight$qValue)<q.cut,1,0)
  feature.weight$selected.lowcut=ifelse(as.numeric(feature.weight$NES)>low.weight.cut&as.numeric(feature.weight$pValue)<low.p.cut&as.numeric(feature.weight$qValue)<low.q.cut,1,0)
  colnames(feature.weight)[feature.col]="Feature.List"
  feature.weight$ENSG=sub("\\:.*$","",feature.weight$Feature.List)
  feature.weight$activity=ifelse(grepl("UP_Level1",feature.weight$Feature.List),1,
                                 ifelse(grepl("UP\\_",feature.weight$Feature.List),2,
                                        ifelse(grepl("DN\\_Level1",feature.weight$Feature.List),-1,
                                               ifelse(grepl("DN\\_",feature.weight$Feature.List),-2,0))))
  up1only=unique(feature.weight$ENSG[feature.weight$activity == 1 & feature.weight$selected==1][which(!feature.weight$ENSG[feature.weight$activity == 1 & feature.weight$selected==1] %in% feature.weight$ENSG[feature.weight$activity == 2 & feature.weight$selected.lowcut==1])])
  dn1only=unique(feature.weight$ENSG[feature.weight$activity == -1 & feature.weight$selected==1][which(!feature.weight$ENSG[feature.weight$activity == -1 & feature.weight$selected==1] %in% feature.weight$ENSG[feature.weight$activity == -2 & feature.weight$selected.lowcut==1])])
  feature.length=stack(lengths(feature.list.train))
  feature.length=feature.length[str_count(feature.length$ind,":")==1,]
  feature.length=separate(data = feature.length, col = ind, into = c("ENSG", "ALT"), sep = ":")
  highlevel.up=feature.length[grepl("UP_Level([2-9]|10)",feature.length$ALT) & feature.length$values>=highlevel.minsize,]
  highlevel.dn=feature.length[grepl("DN_Level([2-9]|10)",feature.length$ALT) & feature.length$values>=highlevel.minsize,]
  ENSG.contraversal=unique(c(intersect(feature.weight$ENSG[feature.weight$selected==1 & feature.weight$activity>0],feature.weight$ENSG[feature.weight$selected.lowcut==1 & feature.weight$activity<0]),
                             intersect(feature.weight$ENSG[feature.weight$selected.lowcut==1 & feature.weight$activity>0],feature.weight$ENSG[feature.weight$selected==1 & feature.weight$activity<0])))
  feature.weight.selected=feature.weight[feature.weight$selected==1,] #this step is critical for the following filtering determination. do not delete.
  feature.weight.selected$filter=ifelse(feature.weight.selected$activity==1 & (feature.weight.selected$ENSG %in% intersect(up1only,unique(highlevel.up$ENSG))),1,
                                        ifelse(feature.weight.selected$activity == -1 & (feature.weight.selected$ENSG %in% intersect(dn1only, unique(highlevel.dn$ENSG))),1,
                                               ifelse(feature.weight.selected$ENSG %in% ENSG.contraversal,2,0)))
  print (paste0(sum(feature.weight.selected$filter==1)," Genotypes identified as low-level OE/UE but not high-level OE/UE are found to be significant"))
  print (paste0(sum(feature.weight.selected$filter==2)," Genotypes identified as both increased and decreased activity of the same gene are found to be significant"))
  print (paste0(sum(feature.weight.selected$filter %in% filter.criteria)," features will be filered from ",nrow(feature.weight.selected)," significant features"))
  feature.weight.selected=feature.weight.selected[!feature.weight.selected$filter %in% filter.criteria,]
  colnames(feature.weight.selected)[feature.col]="Genotype"
  return(feature.weight.selected[,!(colnames(feature.weight.selected) %in% c("ENSG","activity","filter","selected","selected.lowcut"))])
}
#calculate sensitive/resistant igenSig scores based on NES of sensitive/resistant KS tests and precalculation of genotype overlaps using square fomular
cal.genSig<-function(gensig.model=NULL,file.sen=NULL,file.res=NULL,preCalmatrix,genotype.list,trainset,outfile,minsize=10,weight.cut=0,low.weight.cut=0,p.cut=0.01,low.p.cut=0.05,q.cut=0.25,low.q.cut=1,highlevel.minsize=10,power=1,root=1,ECNpenalty=1,rm.equivocal=TRUE,correct.train=T){ #please provide filtered genotype.list based on screened cell lines for each drug if rm.equivocal=TRUE
  if (typeof(gensig.model)=="NULL"){
    nes.sen.input<-read.delim(file.sen,stringsAsFactors = F,header = T,sep="\t")
    nes.sen.input[is.na(nes.sen.input)] <- 0
    nes.res.input<-read.delim(file.res,stringsAsFactors = F,header = T,sep="\t")
    nes.res.input[is.na(nes.res.input)] <- 0
    if (rm.equivocal){#remove genes with both positive activity (i.e. upregulation) and negative activitty (i.e. downregulation) in significant genotypes and genes with high level (2|3) over or under expression genotypes with more than ksMinSize samples, but only level 1 genotypes are significant
      nes.sen.sub=filter.feature(feature.weight=nes.sen.input,feature.list.train=lapply(genotype.list,function(x) x[x %in% trainset]),weight.cut=weight.cut,low.weight.cut=low.weight.cut,p.cut=p.cut,low.p.cut=low.p.cut,q.cut=q.cut,low.q.cut=low.q.cut,highlevel.minsize=highlevel.minsize)
      nes.res.sub=filter.feature(feature.weight=nes.res.input,feature.list.train=lapply(genotype.list,function(x) x[x %in% trainset]),weight.cut=weight.cut,low.weight.cut=low.weight.cut,p.cut=p.cut,low.p.cut=low.p.cut,q.cut=q.cut,low.q.cut=low.q.cut,highlevel.minsize=highlevel.minsize)
    }else{
      nes.sen.sub<-nes.sen.input[nes.sen.input$NES>weight.cut & nes.sen.input$pValue<p.cut & nes.sen.input$qValue<q.cut,-c(3,4)]
      nes.res.sub<-nes.res.input[nes.res.input$NES>weight.cut & nes.res.input$pValue<p.cut & nes.res.input$qValue<q.cut,-c(3,4)]
    }
    nes.sen <- as.matrix(nes.sen.sub[-1])
    row.names(nes.sen) <- nes.sen.sub$Genotype
    nes.res <- as.matrix(nes.res.sub[-1])
    row.names(nes.res) <- nes.res.sub$Genotype
  }else{
    list
    nes.sen=gensig.model$nes.sen
    nes.res=gensig.model$nes.res
    print (paste("loss of significant sensitive genotypes: ",length(which(!row.names(nes.sen) %in% names(genotype.list)))/length(row.names(nes.sen))*100,"%",sep=" "))
    print (paste("loss of significant resistant genotypes: ",length(which(!row.names(nes.res) %in% names(genotype.list)))/length(row.names(nes.res))*100,"%",sep=" "))
  }
  sen.list=list()
  res.list=list()
  result<-data.frame(matrix(ncol = 3, nrow = 0))
  subject.all=unique(unlist(genotype.list))
  for(i in 1:length(preCalmatrix)){
    cellid=names(preCalmatrix)[i]
    if (length(preCalmatrix[[i]])==0|!cellid %in% subject.all) next
    tmp.epsilon<-data.frame(Genotypes=names(preCalmatrix[[i]]),Epsilon=preCalmatrix[[i]],stringsAsFactors = F)
    tmp.epsilon[,"Epsilon"]=(as.numeric(tmp.epsilon[,"Epsilon"]))^root
    ECN=nrow(tmp.epsilon)/mean(as.numeric(tmp.epsilon[,"Epsilon"]),trim=0.3) #Previous Version: ECN=sum(1/as.numeric(tmp.epsilon[,"Epsilon"])) 
    sen.epsilon<-tmp.epsilon[tmp.epsilon$Genotypes %in% row.names(nes.sen),]
    res.epsilon<-tmp.epsilon[tmp.epsilon$Genotypes %in% row.names(nes.res),]
    if (correct.train==T & paste("X",cellid,sep="") %in% colnames(nes.sen)){
      sen.weight<- nes.sen[which(row.names(nes.sen) %in% sen.epsilon$Genotypes),c("NES",paste("X",cellid,sep="")),drop=F]
      sen.weight <- as.data.frame(replace(sen.weight[,"NES"], which(sen.weight[,paste("X",cellid,sep="")]!=0),sen.weight[which(sen.weight[,paste("X",cellid,sep="")]!=0),paste("X",cellid,sep="")]),stringsAsFactors=FALSE)
    }else{
      sen.weight<- as.data.frame(nes.sen[which(row.names(nes.sen) %in% sen.epsilon$Genotypes),"NES"],stringsAsFactors=FALSE)
    }
    sen.weight$Genotypes <- as.character(row.names(sen.weight))
    colnames(sen.weight)=c("NES","Genotypes")
    sen.data<-as.matrix(merge(sen.epsilon,sen.weight,by.x="Genotypes",by.y="Genotypes"))
    sen.data=cbind(sen.data,normWeight=(as.numeric(sen.data[,"NES"])^power)/as.numeric(sen.data[,"Epsilon"]))
    score.sensitive<-sum(as.numeric(sen.data[,"normWeight"]))/(ECN^ECNpenalty)
    sen.data=as.data.frame(sen.data,stringsAsFactors=FALSE)
    colnames(sen.data)[-1]=c(paste(cellid,colnames(sen.data[-1]),sep=":"))
    sen.list[[paste(cellid,sep="")]]<-sen.data
    if (correct.train==T & paste("X",cellid,sep="") %in% colnames(nes.res)){
      res.weight<- nes.res[which(row.names(nes.res) %in% res.epsilon$Genotypes),c("NES",paste("X",cellid,sep="")),drop=F]
      res.weight <- as.data.frame(replace(res.weight[,"NES"], which(res.weight[,paste("X",cellid,sep="")]!=0),res.weight[which(res.weight[,paste("X",cellid,sep="")]!=0),paste("X",cellid,sep="")]),stringsAsFactors=FALSE)
    }else{
      res.weight<- as.data.frame(nes.res[which(row.names(nes.res) %in% res.epsilon$Genotypes),"NES"],stringsAsFactors=FALSE)
    }
    res.weight$Genotypes <- as.character(row.names(res.weight))
    colnames(res.weight)=c("NES","Genotypes")
    res.data<-as.matrix(merge(res.epsilon,res.weight,by.x="Genotypes",by.y="Genotypes"))
    res.data=cbind(res.data,normWeight=(as.numeric(res.data[,"NES"])^power)/as.numeric(res.data[,"Epsilon"]))
    score.resistant<-sum(as.numeric(res.data[,"normWeight"]))/(ECN^ECNpenalty)  
    res.data=as.data.frame(res.data,stringsAsFactors=FALSE)
    colnames(res.data)[-1]=c(paste(cellid,colnames(res.data[-1]),sep=":"))
    res.list[[paste(cellid,sep="")]]<-res.data
    result[nrow(result) + 1,]<-c(cellid,score.sensitive,score.resistant)
    if(i %% 100==0){
      print(paste("Processed ",i," subjects",sep=""))
    }
  }
  colnames(result)<-c("CosmicID","GeneSig.sensitive","GeneSig.resistant")
  gensig.model=list(nes.sen=nes.sen, nes.res=nes.res, parameters=list(minsize=minsize,weight.cut=weight.cut,low.weight.cut=low.weight.cut,p.cut=p.cut,low.p.cut=low.p.cut,q.cut=q.cut,low.q.cut=low.q.cut,highlevel.minsize=highlevel.minsize,power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal),sen.list=sen.list, res.list=res.list)
  save(gensig.model, file = paste(tools::file_path_sans_ext(outfile),".rda",sep=""))
  write.table(result,file=outfile,sep="\t",row.names=FALSE)
}
#batch build iGenSig models and calculate iGenSig scores for GDSC dataset
batchCalGenSig.GDSC<-function(GDSC.genotype.list,GDSC.preCalfile,drug.vec,GDSC.drugData,GDSC.testsetdir,GDSC.weightdir,GDSC.gensigdir,p.cut=1,low.p.cut=1,q.cut=0.1,low.q.cut=0.3,power=1,root=1,ECNpenalty=1,rm.equivocal=T,correct.train=T,allSubject=T,use.genSig.model=T,skip.completed=T){
  if (grepl(".gmt",GDSC.preCalfile)){
    preCalmatrix=read_preCal(GDSC.preCalfile) 
  }else if (grepl(".RData",GDSC.preCalfile,ignore.case=TRUE)){
    preCalmatrix=mget(load(GDSC.preCalfile, envir=(tmp<- new.env())), envir=tmp)$result$feature.redundancy
  }else{
    print ("precalfile name wrong")
  }
  files<-list.files(path = GDSC.weightdir, pattern = "sensitive.xls")
  if (length(files)==0) stop(paste("weight files not detected in",GDSC.weightdir,sep=" "))
  runs<-c()
  for (file in files){
    drugID.call=as.numeric(unlist(str_extract_all(file, "(?<=drugID\\_)\\d+"))[1])
    if (drugID.call %in% drug.vec){
      prefix<-sub(".sensitive.xls", "", file)
      runs=c(prefix,runs)  
    }
  }
  runs=unique(runs)
  for (i in 1:length(runs)){
    prefix=runs[i]
    drugID.call=as.numeric(unlist(str_extract_all(prefix, "(?<=drugID\\_)\\d+"))[1])
    permuID.call=as.numeric(unlist(str_extract_all(prefix, "(?<=permuID\\_)\\d+"))[1])
    file.sen=paste(GDSC.weightdir,"/",prefix,".sensitive.xls",sep="")
    file.res=paste(GDSC.weightdir,"/",prefix,".resistant.xls",sep="")
    outfile<-paste(GDSC.gensigdir,"/",prefix,"_genSig.xls",sep="")
    if(file.exists(outfile)& skip.completed==T){
      print(paste(outfile," found, calculation skipped",sep=""))
      next
    }
    print(paste("Building iGeneSig model for",prefix,sep=" "))
    if (allSubject==T){
      drug.genotype.list=GDSC.genotype.list
    }else{
      drug.genotype.list=lapply(GDSC.genotype.list,filter,filter=GDSC.drugData$COSMIC_ID[GDSC.drugData$DRUG_ID==drugID.call],inlist=TRUE)
    }
    drug.subjects=GDSC.drugData$COSMIC_ID[GDSC.drugData$DRUG_ID==drugID.call]
    if (is.na(permuID.call)){
      print (paste0("permuID not found for file prefix ",prefix,", using all subjects as training set"))
      trainset=drug.subjects
    }else{
      testset=read.delim(paste0(GDSC.testsetdir,"/testingSet_drugID",drugID.call,"_permu",permuID.call,".txt"),stringsAsFactors = F,check.names = F,header = T,sep="\t")
      trainset=drug.subjects[!drug.subjects %in% testset$COSMIC_ID]
    }
    if (use.genSig.model==T & file.exists(paste(tools::file_path_sans_ext(outfile),".rda",sep=""))){
      file.model=paste(tools::file_path_sans_ext(outfile),".rda",sep="")
      gensig.model=mget(load(file.model, envir=(tmp<- new.env())), envir=tmp)$gensig.model
    }else{
      gensig.model=NULL 
    }
    cal.genSig(gensig.model=gensig.model,file.sen=file.sen,file.res=file.res,preCalmatrix=preCalmatrix,genotype.list=drug.genotype.list,trainset=trainset,outfile=outfile,p.cut=p.cut,low.p.cut=low.p.cut,q.cut=q.cut,low.q.cut=low.q.cut,power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal,correct.train=correct.train)
  }
}
#####################################################################################
### Modules for apply iGenSig module to validation datasets
#####################################################################################
#Batch Calculate GenSig scores for validation sets based on GDSC KS result models
batchCalGenSig.validationset<-function(genotype.list,preCalfile,drug.vec,GDSC.gensigdir,result.gensigdir,minSubject=1){
  if (grepl(".gmt",preCalfile)){
    preCalmatrix=read_preCal(preCalfile) 
  }else if (grepl(".RData",preCalfile,ignore.case=TRUE)){
    preCalmatrix=mget(load(preCalfile, envir=(tmp<- new.env())), envir=tmp)$result$feature.redundancy
  }else{
    print ("precalfile name wrong")
  }
  files<-list.files(path = GDSC.gensigdir, pattern = "genSig.rda")
  if (length(files)==0){
    stop(paste("weight files not detected in",GDSC.gensigdir,sep=" "))
  }
  runs<-list()
  count=0
  for (i in 1:length(files)){
    drugID.call<-as.numeric(unlist(str_extract_all(files[i], "(?<=drugID_)\\d+"))[1])
    prefix<-sub("_genSig.rda", "", files[i])
    if (drugID.call %in% drug.vec){
      count=count+1
      runs[[count]]=c(drugID.call,prefix) 
    }
  } 
  runs=unique(runs)
  if (length(runs)==0){
    stop(paste("the iGenSig model for the drug.vec is not found in the GDSC genSig folder:",GDSC.gensigdir))
  }
  for (i in 1:length(runs)){
    drugID.call=runs[[i]][1]
    prefix=runs[[i]][2]
    file.model=paste(GDSC.gensigdir,"/",prefix,"_genSig.rda",sep="")
    gensig.model=mget(load(file.model, envir=(tmp<- new.env())), envir=tmp)$gensig.model
    list2env(gensig.model$parameters, env = environment())
    outfile<-paste(result.gensigdir,"/",prefix,"_genSig.xls",sep="")
    print(paste("Computing GeneSig for",prefix,sep=" "))
    if(!file.exists(outfile)){
      cal.genSig(gensig.model=gensig.model,preCalmatrix=preCalmatrix,genotype.list=genotype.list,trainset=NA,outfile=outfile,p.cut=p.cut,low.p.cut=low.p.cut,q.cut=q.cut,low.q.cut=low.q.cut,power=power,root=root,ECNpenalty=ECNpenalty,rm.equivocal=rm.equivocal)
    }else{
      print(paste(outfile," found, calculation skipped",sep=""))
    }
  } 
}
model.GDSC2trial<-function (GDSC.gensigdir,Trial.gensigdir,Trial.drug,Trial.genotype.list,TCGA.catype.genotype.list,Trial.preCalfile,
                            sample.col,response.col,event.col,time.col,calGenSig=T,
                            by.cluster=TRUE,cutTree.method="dynamic",dynamic.method="hybrid",cutTree.depth=2,minClusterSize=40){
  log.dir=paste0(Trial.gensigdir,"/",cutTree.method,dynamic.method,"MinSize",minClusterSize,"Depth",as.integer(cutTree.depth))
  dir.create(log.dir)
  if (!file.exists(Trial.preCalfile)){
    batch_calSimilarity_genolist(subject.genotype.list=Trial.genotype.list,REF.genotype.list=TCGA.catype.genotype.list,subject.preCalfile=Trial.preCalfile,by.cluster=by.cluster,cutTree.method=cutTree.method,dynamic.method=dynamic.method,minClusterSize=minClusterSize,split.depth=cutTree.depth,log.dir=log.dir,visualize.cluter=F)
  }
  if (calGenSig==T){
    batchCalGenSig.validationset (
      genotype.list=Trial.genotype.list,
      preCalfile=Trial.preCalfile,
      drug.vec=Trial.drug,
      GDSC.gensigdir=GDSC.gensigdir,
      result.gensigdir=Trial.gensigdir
    ) 
  }
}
#####################################################################################
### Benchmark cell line iGenSig models based on sensitive GenSig scores
#####################################################################################
benchmark.genSig<-function(iGenSig.resultfile,testset.dir,drugData,score="GeneSig.sensitive",response=c("AUC","ActArea")){  
  match.arg(response)
  drugID.call=as.numeric(unlist(str_extract_all(iGenSig.resultfile, "(?<=drugID_)\\d+"))[1])
  permuID.call=as.numeric(unlist(str_extract_all(iGenSig.resultfile, "(?<=permuID_)\\d+"))[1])
  result.plot<-read.delim(iGenSig.resultfile,stringsAsFactors = F,check.names = F,header = T,sep="\t")
  result.plot$GeneSig.sensitive=normalize(result.plot$GeneSig.sensitive)
  result.plot$GeneSig.resistant=normalize(result.plot$GeneSig.resistant)
  subData<-drugData[!is.na(drugData$DRUG_ID) & drugData$DRUG_ID==drugID.call,]
  if (response=="AUC"){
    cutoffs=define.response.AUC(subData[,response])
    subData$label<-ifelse(subData[,response]<cutoffs[1],"sen",ifelse(subData[,response]>cutoffs[2],"res","other"))
    if (is.na(permuID.call)){
      result.plot$testSet<-"Y"
    }else{
      myTest<-read.delim(paste(testset.dir,"/testingSet_drugID",drugID.call,"_permu",permuID.call,".txt",sep=""),stringsAsFactors = F,check.names = F,header = T,sep="\t")
      result.plot$testSet=ifelse(result.plot$CosmicID %in% myTest$COSMIC_ID,"Y","N")
    }
    result.plot=merge(result.plot,subData,by.x="CosmicID",by.y="COSMIC_ID")
    drug.name<-subData$DRUG_NAME[subData$DRUG_ID==drugID.call][1]
  }else if (response=="ActArea"){
    cutoffs=define.response.ActArea(subData[,response])
    subData$label<-ifelse(subData[,response]>cutoffs[1],"sen",ifelse(subData[,response]<cutoffs[2],"res","other"))
    result.plot$testSet="Y"
    result.plot=merge(result.plot,subData,by.x="CosmicID",by.y="ACHID")
    drug.name<-subData$Compound[subData$DRUG_ID==drugID.call][1]
  }
  if (is.na(permuID.call)){
    plot.gensig<-ggplot(result.plot[order(result.plot$label),],aes(x=GeneSig.resistant,y=GeneSig.sensitive)) + 
      geom_point(shape=16,size=1,aes(color=label)) + 
      theme(text = element_text(size=6),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black",size=0.3), axis.ticks = element_line(colour = "black",size=0.3), legend.position="bottom",legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), legend.key=element_blank()) + 
      scale_color_manual(values = c("sen" = "red", "res" = "deepskyblue", "mid" = "grey"))
  }else{
    plot.gensig<-ggplot(result.plot[order(result.plot$label,result.plot$testSet),],aes(x=GeneSig.resistant,y=GeneSig.sensitive)) + 
      geom_point(shape=21,stroke=0.2,size=1,aes(color=testSet,fill=label)) + 
      theme(text = element_text(size=6),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black",size=0.3), axis.ticks = element_line(colour = "black",size=0.3), legend.position="bottom",legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), legend.key=element_blank()) + scale_color_manual(breaks = c("N", "Y"), values=c(rgb(255, 255, 255, alpha=0, names = NULL, maxColorValue = 255), "black")) + 
      scale_fill_manual(values = c("sen" = "red", "res" = "deepskyblue", "mid" = "grey"))
  }
  result.test<-result.plot[result.plot$testSet=="Y",]
  result.test$label<-ifelse(result.test$label=="sen","sen","other")
  result.plot=result.plot[order(result.plot$testSet),]
  if (sum(result.test$label=="sen")>=3 & sum(result.test$label=="other")>=3){
    test.roc<-roc(result.test$label,result.test[,score],levels=c("other","sen"),quiet=FALSE,direction = "<")
    response.sum=table(result.test$label)
    opt.cut<-coords(test.roc, x="best", input="threshold", best.method="youden",best.weights = c(2,0.2),transpose=TRUE)#response.sum["sen"]/(response.sum["sen"]+response.sum["other"]))
    plot.score<-ggplot(result.plot,aes(x=result.plot[,score],y=result.plot[,response])) + 
      geom_point(shape=21,fill="grey",stroke=0.2, size=1,aes(color=result.plot[,"testSet"])) + 
      geom_hline(yintercept=as.numeric(cutoffs[1]), linetype="dotted", colour = "red", size=0.3, alpha=0.8) + 
      geom_hline(yintercept=cutoffs[2], linetype="dotted", colour = "deepskyblue", size=0.3, alpha=0.8) + 
      geom_vline(aes(xintercept=opt.cut[1]), color="darkgrey", linetype="dotted", size=0.3)+
      xlab(score)+
      theme(text = element_text(size=6),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black",size=0.3), axis.ticks = element_line(colour = "black",size=0.3), legend.position="bottom",legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), legend.key=element_blank()) + 
      scale_color_manual(breaks = c("N", "Y"), values=c(rgb(255, 255, 255, alpha=0, names = NULL, maxColorValue = 255), "black"))
    plot.roc<-ggroc(test.roc,color="red")+ 
      geom_point(aes(x=opt.cut[2], y=opt.cut[3]), colour="black", size=0.3)+
      geom_text(aes(x=opt.cut[2], y=opt.cut[3],label=paste("Cutoff:",round(opt.cut[1],digits=2),"(",percent(opt.cut[2]),",",percent(opt.cut[3]),")",sep="")),hjust=0, vjust=-0.5,colour="black",size=1.5)+
      geom_text(aes(x=0.75, y=0.25,label=paste("AUROC=",round(test.roc$auc,digits = 3))),hjust=0, vjust=-0.5,colour="black",size=1.5)+
      theme(text = element_text(size=6),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black",size=0.3), axis.ticks = element_line(colour = "black",size=0.3), legend.position="bottom",legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), legend.key=element_blank())
  }else{
    plot.score<-ggplot(result.plot,aes(x=result.plot[,score],y=result.plot[,response])) + 
      geom_point(shape=21,fill="grey",stroke=0.2, size=1,aes(color=result.plot[,"testSet"])) + 
      geom_hline(yintercept=as.numeric(cutoffs[1]), linetype="dotted", colour = "red", size=0.3, alpha=0.8) + 
      geom_hline(yintercept=cutoffs[2], linetype="dotted", colour = "deepskyblue", size=0.3, alpha=0.8) + 
      xlab(score)+
      theme(text = element_text(size=6),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black",size=0.3), axis.ticks = element_line(colour = "black",size=0.3), legend.position="bottom",legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), legend.key=element_blank()) + 
      scale_color_manual(breaks = c("N", "Y"), values=c(rgb(255, 255, 255, alpha=0, names = NULL, maxColorValue = 255), "black"))
    plot.roc=plot.new()
  }
  plots.list<-list(plot.gensig,plot.score,plot.roc)
  names(plots.list)<-c(paste("DrugID_",drugID.call,"_PermuID_",permuID.call,"_GenSig_Plot",sep=""),paste("DrugID_",drugID.call,"_PermuID_",permuID.call,"_",score,"_Plot",sep=""),paste("DrugID_",drugID.call,"_PermuID_",permuID.call,"_ROC_Plot",sep=""))
  if (exists("test.roc")){
    tmp.benchmark <- c(drugID.call,as.character(drug.name),permuID.call,test.roc$auc) 
  }else{
    tmp.benchmark <- c(drugID.call,as.character(drug.name),permuID.call,NA)
  }
  return(list(plots.list,tmp.benchmark))
}
batch.benchmarkGenSig <- function(gensig.dir,testset.dir=NULL,drugData,dataset=c("GDSC","CCLE")){
  match.arg(dataset)
  response=ifelse(dataset=="GDSC","AUC","ActArea")
  files<-list.files(path = gensig.dir, pattern = "permuID_\\d+_genSig.xls")
  benchmark.df<-data.frame(matrix(ncol = 4, nrow = 0))
  plots<-list()
  for (file in files){
    drugID.call=as.numeric(unlist(str_extract_all(file, "(?<=drugID_)\\d+"))[1])
    permuID.call=as.numeric(unlist(str_extract_all(file, "(?<=permuID_)\\d+"))[1])
    if (drugID.call %in% unique(drugData$DRUG_ID)){
      tmp.results<-benchmark.genSig(iGenSig.resultfile=paste(gensig.dir,file,sep="/"),testset.dir=testset.dir,drugData=drugData,score="GeneSig.sensitive",response=response)
    }
    benchmark.df[nrow(benchmark.df) + 1,]<-tmp.results[[2]]
    plots<-c(plots,tmp.results[[1]])
  }
  colnames(benchmark.df)<-c("drugID","drugName","permuID","AUC")
  write.table(benchmark.df,file=paste(gensig.dir,"/benchmark.genSig.result.xls",sep=""),sep="\t", row.names=FALSE, col.names = TRUE)
  if (length(plots)>0){
    drugID.prev=0
    for (i in 1:length(plots)){
      drugID.call<-as.numeric(unlist(str_extract_all(names(plots)[i], "(?<=DrugID_)\\d+"))[1])
      permuID.call<-as.numeric(unlist(str_extract_all(names(plots)[i], "(?<=PermuID_)\\d+"))[1])
      if (drugID.call!=drugID.prev){
        drug.name<-drugData$DRUG_NAME[drugData$DRUG_ID==drugID.call][1]
        if (drugID.prev!=0){
          dev.off()
        }
        pdf(file=paste(gensig.dir,"/DrugID_",drugID.call,"_benchmark.genSig.plot.pdf",sep=""),width=4, height=6)
        print (paste("Printing DrugID_",drugID.call,"_benchmark.genSig.plot.pdf",sep=""))
        grid.newpage() # Open a new page on grid device
        pushViewport(viewport(height=1,layout = grid.layout(round(length(plots)/3), 3))) # Assign to device viewport with plotN by 3 grid layout
      }
      i.file=(i-1) %% 18 +1
      row=(i.file-1) %/%3 +1
      col=(i-1) %% 3 +1
      print(plots[[i]]+theme(legend.position="none"), vp = viewport(layout.pos.row = row, layout.pos.col = col))
      drugID.prev=drugID.call
    }
    dev.off()   
  }
}