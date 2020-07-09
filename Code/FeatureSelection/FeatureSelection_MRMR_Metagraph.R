library(RANKS) 
library("e1071")
library(caret)
library(ggplot2)

library(clusterSim)
library(plyr)
library('mRMRe')

datasets <- c("IntAct", "NCBI", "STRING")

diseases <- c( "prostateCancer","breastCancer","diabetesMellitus","obesity","lungCancer","alzheimer", "colorectalCancer")
memory.limit(size=50000)
sized<-length(datasets)
k=5
dis=1
z=1
j=1

path1<-"Data_for_Feature_Selection\\"
pathWrite<-""
for(dis in 1:7){
  disease=diseases[dis]
  labelfile<-paste0(path1,"2018OMIMDiseaseLabels/")
  
  for(z in 1:2){
    dataset=datasets[z]
    print(paste0("Currently handling: ",disease," ",dataset))

      for (j in 1:k){
        dataset.name <-paste0("..\\..\\","FinalRepresentations\\",disease, dataset,j,"MatrixforSVMvector.txt")
        keyw = read.table(dataset.name,header =FALSE,sep="\t")
        K <- as.matrix(keyw)
     
        n <- nrow(K)
        p <- numeric(n)

        testIndices <- paste0(path1,"TestIndices/", dataset,disease, j,"indicesofTestSet.txt")  

        testIndices= read.table(testIndices,header =FALSE,sep="\n")
        testIndices<-as.array(testIndices$V1)

        xtest<-as.data.frame(K[testIndices,])#test indices start with 1 so no need any increase or decrease by 1
        
        xtrain <- K[-(testIndices),]
      
        allLabels <-paste0(labelfile , dataset , "SVMclassfor",disease,".txt")
        allLabel= read.table(allLabels,header =FALSE,sep="\n")
        
        L <- as.matrix(allLabel)
        
        ytrain <-as.matrix( L[-(testIndices),])
        ytest<-L[testIndices,]##test indices start with 1 so no need any increase or decrease by 1  
        ytrain<-as.numeric(ytrain)
        ytrain<-as.data.frame(ytrain)

        trainMean <- apply(xtrain,2,mean)
        trainSd <- apply(xtrain,2,sd)
        xtrain<- sweep(sweep(xtrain, 2L, trainMean), 2, trainSd, "/") # using the default "-" to subtract mean column-wise   
        # ## centered AND scaled
        xtest<-sweep(sweep(xtest, 2L, trainMean), 2, trainSd, "/")
        trainsize<-nrow(xtrain)
        merged<-rbind(xtrain,xtest)
        merged<-merged[,!(colSums(is.na(merged)))]
        xtrain<-merged[1:trainsize,]
        xtest<-merged[(trainsize+1):nrow(merged),]
        # calculate correlation matrix
        # which(is.na(keywords))
        #>>>>>>>>>>REMOVE ZERO VARIANCE VARIABLES
        zv <- apply(xtrain, 2, function(x) length(unique(x)) == 1)
        dfr <- xtrain[, !zv]
        #>>>>>>>>>>>>>>>>>
        correlationMatrix <- cor(dfr[,1:ncol(dfr)],use="complete.obs")
        # summarize the correlation matrix
        #print(correlationMatrix)
        # find attributes that are highly corrected (ideally >0.75)
        highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
        # print indexes of highly correlated attributes
        xtrain<-xtrain[,-highlyCorrelated ]

        if (nrow(xtrain) != nrow(ytrain)) 
          stop("train data and label matrices do not agree")
        colnames(ytrain) <- "Labels"
        DATA=cbind(xtrain,ytrain)
        DATA<-as.data.frame(DATA)
        names(DATA)
  
        showClass("mRMRe.Filter")
        ## build an mRMRe.Data object
        myData <- mRMR.data(data = DATA)

        filter <- mRMR.classic("mRMRe.Filter", data = myData, target_indices = ncol(DATA),feature_count=ncol(DATA)-1)
        scores<-unlist(scores(filter))
        causality<- unlist(causality(filter))
        featureNames<-featureNames(filter)
        filter <- mRMR.classic("mRMRe.Filter", data = myData, target_indices = ncol(DATA),feature_count=ncol(DATA)-1)
        scores<-unlist(scores(filter))
        causality<- unlist(causality(filter))
        featureNames<-featureNames(filter)
        wrtscores<-cbind(featureNames[1:length(scores)], scores)
        wrtcausality<-cbind(featureNames[1:length(causality)], causality)
        
        ScoresSorted<-wrtscores[order(wrtscores[,2],decreasing = TRUE),]
        causalitySorted<-wrtcausality[order(wrtcausality[,2],decreasing = TRUE),]
        
        path<- paste0(pathWrite,"mRMR/")
        dir.create(path,recursive = TRUE,showWarnings = FALSE)
        write.table(ScoresSorted, file = paste0(path,disease,dataset,j,"KeywordsEmbeddingsUsefulSubset.txt"),append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE, sep = "\t")
        write.table(causalitySorted, file = paste0(path,disease,dataset,j,"causalitySorted.txt"),append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE, sep = "\t")
        
        topvector<-seq.int(10, 500, 10) #from , to, by =10
        for(to in 1:length(topvector)){
          path<- paste0(pathWrite,"mRMRTOP",topvector[to],"/")
          dir.create(path,recursive = TRUE,showWarnings = FALSE)
          
          write.table(ScoresSorted[1:topvector[to],], file = paste0(path,disease,dataset,j,"KeywordsEmbeddingsUsefulSubset.txt"),append = FALSE,row.names = FALSE,col.names = FALSE,quote = FALSE, sep = "\t")
        }        
      }
      }
    }
    
    
    
  
