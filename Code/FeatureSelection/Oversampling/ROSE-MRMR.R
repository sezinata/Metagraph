library(RANKS) 
library("e1071")
library(caret)
library(ggplot2)
library('stringr')
library(clusterSim)
library('ROSE')
library('AUC')

networks<-c("PPI","PPIK")

datasets <- c("IntAct", "NCBI", "STRING")

diseases <- c( "prostateCancer","breastCancer","diabetesMellitus","obesity","lungCancer","alzheimer", "colorectalCancer")
features<-array()
featuresFolders<-array()
ifeat<-1
for(to in seq(10, 500, by = 10) ){
  features[ifeat]<- paste0("mRMRTOP",to)
  ifeat<-ifeat+1
}


OversamplingMethods<-c("ROSE","SMOTE")
OversamplingMethod<-OversamplingMethods[1]
parameters<-as.numeric(c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'))
parametersFolder<-c('01','02','03','04','05','06','07','08','09')
sized<-length(datasets)
k=5
pv=1
q=1
net=1
feat=1
dis=1
z=1
j=1
net=1
empty=0
para=1

empties<-data.frame()

path1<-"..\\Data_for_Feature_Selection"
pathRead<-"..\\"
pathWrite<-""
for(dis in 1:7){
  disease=diseases[dis]
  labelfile<-paste0(path1,"\\2018OMIMDiseaseLabels/")
  
  for(z in 1:2){
    dataset=datasets[z]
    
    for(feat in 1:length(features)){
      feature<-features[feat]
      featuresFolder<-featuresFolders[feat]
      print(paste0("Currently handling: ",disease," ",dataset," ",feature))
      
      result<-array()  

      for(para in 1:9){
        parameter<-parameters[para]
        Oversampling<-paste0(OversamplingMethod,parametersFolder[para])
        
        
        sum=0
        
        for (j in 1:k){
          dataset.name <-paste0("..\\..\\..\\","FinalRepresentations\\",disease, dataset,j,"MatrixforSVMvector.txt")
          
          metagraph = read.table(dataset.name,header =FALSE,sep="\t")
          K <- metagraph
          path<- paste0(pathRead,feature,"/")
          
          selectedFeatures= read.table(file = paste0(path,disease,dataset,j,"KeywordsEmbeddingsUsefulSubset.txt"),header =FALSE,sep="\t")
          Selectedvar_names<-as.vector(selectedFeatures$V1)
          if(length(Selectedvar_names)>0){

            filtered=as.data.frame( K[,Selectedvar_names])  #each keyword individually a list element with its values for all proteins
            n <- nrow(filtered)
            testIndices <- paste0(path1,"\\TestIndices/", dataset,disease, j,"indicesofTestSet.txt")  

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
            
            if (nrow(xtrain) != nrow(ytrain)) 
              stop("train data and label matrices do not agree")
            
            colnames(ytrain) <- "Labels"
            
            path<- paste0(pathWrite,feature,"/LogisticRegression/",Oversampling,"/")
            dir.create(path,recursive = TRUE,showWarnings = FALSE)
            if(ncol(xtrain)==1){
              colnames(xtrain)<-'xtrain'
              DATA=cbind(xtrain, ytrain)
              DATA<-as.data.frame(DATA)
              
              names(DATA)
              data.rose <- ROSE(Labels~xtrain, data = DATA, p=parameter)$data  #creates with prob xx oversampling disease
              
              model <-lm(Labels~xtrain, data =  data.rose)
              new <- data.frame(xtrain = xtest) 
              glmpredict<- predict(model, newdata = new)
              ntrain <- data.frame(xtrain = xtrain) 
              colnames(ntrain)<-'xtrain'
              glmpredictrain<-predict(model, newdata=ntrain)
              label_test <- factor(ytest)
              sum<-sum+auc(roc(glmpredict,label_test))#BestOversampling  "ROSExx"
              
            }
            else{
              DATA=cbind(xtrain, ytrain)
              DATA<-as.data.frame(DATA)
              
              xtest<-as.data.frame(xtest)
              xtrain<-as.data.frame(xtrain)
              forrose<-array()
              for (vec in 1:ncol(xtrain)) {
                forrose[vec]<-paste0('V',vec)
              }
              colnames(xtest)<-forrose
              colnames(xtrain)<-forrose
              
              forrose[vec+1]='Labels'
              f <-paste ("Labels~", paste(sprintf("%s",  forrose), collapse="+"))
              colnames(DATA)<-forrose
              data.rose <- ROSE(as.formula(f), data = DATA, p=parameter)$data  #creates with prob xx oversampling disease
              
              model <- glm(Labels~., data = data.rose, family = "gaussian")
              glmpredict<- predict(model, newdata = xtest, type = "response")
              glmpredictrain<-predict(model, newdata=xtrain,  type="response")
              label_test <- factor(ytest)
              sum<-sum+auc(roc(glmpredict,label_test))#BestOversampling 
              
            }
            write.table(glmpredict, file = paste0(path,disease,dataset, j,"Probabilities.txt"),append = FALSE, row.names = FALSE, sep = " ")
            
            write.table(glmpredictrain, file = paste0(path,disease,dataset, j,"trainProbabilities.txt"),append = FALSE, row.names = FALSE, sep = " ")
            write(testIndices, file = paste0(path,disease,dataset, j,"indicesofTestSet.txt"),append = FALSE, sep = "\n")
            write(ytest, file = paste0(path,disease,dataset, j,"Labels.txt"),append = FALSE, sep = "\n")
            write.table(ytrain, file = paste0(path,disease,dataset, j,"TrainLabels.txt"),append = FALSE,quote=FALSE, row.names = FALSE, sep = " ")
          }
          else{
            empty=empty+1
            empties<-as.data.frame(rbind(empties,cbind(disease,dataset,j,feature)))
            print(empties)
            
          }
        }
        
        result[para]<-sum/k
        
        DATA <- matrix(c(parameter,result[para],feature),ncol=3,byrow=TRUE)
        
        write.table(DATA, file = paste0(pathWrite,disease,dataset,"MRMRROSEParameterScores.txt"),quote=FALSE,append = TRUE,col.names = FALSE, row.names = FALSE, sep = " ")
        
      }
    }
    
  }
  
}


