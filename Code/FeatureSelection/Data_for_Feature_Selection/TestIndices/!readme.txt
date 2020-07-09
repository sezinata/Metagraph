Traditional cv i applied. 2018 OMIM LABELs. starts with 1

fold <- do.stratified.cv.data(1:numberoflines, ind.pos, k = k, seed = init.seed)
i=1
for (i in 1:k) {
  x <- c(fold$fold.positives[[i]], fold$fold.non.positives[[i]])
  core.pos <- integer(0)
  for (j in 1:k) 
    if (j == i) {
   # core.pos <- c(core.pos, fold$fold.positives[[j]])
      test<-x
      testfold<-paste0("testfold", j)
      assign(testfold,test)
      write(test, file = paste0("C:/Users/sezin/My Cloud/SameWorksDoneDifferentLabelsandTraditionalCV/TestIndices/",dataset,disease,j,"indicesofTestSet.txt"),
            append = FALSE, sep = "\n")
      break;
    }
}