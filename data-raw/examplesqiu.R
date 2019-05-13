library(htmltab)
library(usethis)

url = "http://users.phhp.ufl.edu/pqiu/research/book/spc/data/"
tableExamples <- htmltab(url, which=1, rm_nodata_cols = T)
namesFileExamples <- tableExamples$Name[grepl(".dat", tableExamples$Name , fixed=TRUE)]

spc_in_spm_examples = c("example82", "example84", "example87", "example65", "example49", "example81", "example83", "example85", "example91", "example93", "example71", "example74a", "example74b", "example74c")
noHeader = c("example49", "example91")
for(nameFile in namesFileExamples){
  fileName = gsub(".dat", "", nameFile)
  if(fileName %in% spc_in_spm_examples){
    header = TRUE
    if(fileName %in% noHeader){
      header = FALSE
    }
    urlFile = paste(url,nameFile,sep="",collapse = "")
    dfExample = read.csv(urlFile,sep="", header=header)
    src = paste("data-raw/",fileName,".csv",sep="",collapse="")

    if (fileName %in% c("example82","example84")){
      df = dfExample[,1:5]
      X.id = rep(seq(1:nrow(df)), ncol(df))
      X = as.numeric(unlist(df)) #convert dataset to vector
      X = X[order(X.id)] #order according the id
      X.id = X.id[order(X.id)] #order the id
      dfExample = data.frame(X = X, X.id = X.id)
    }else if (fileName == "example87"){
      df = dfExample$x
      Y = df[1:14] #first 14 observation in control
      X = df[15:length(df)]
      X.id = 1:length(X)
      Y = c(Y, rep(NA, (length(X)-length(Y))))
      dfExample = data.frame(X = X, X.id = X.id, Y = Y)
    }else if (fileName == "example49") {
      nb = 10 #number of batches
      n = 5 #numbers of observations in each batch
      df1 = dfExample$V1
      Y1 = df1[1:(nb*n)] #N(0,1)
      X1 = df1[(nb*n+1):(2*nb*n)] #N(1,1)

      df2 = dfExample$V2
      Y2 = df2[1:(nb*n)] #N(0,1)
      X2 = df2[(nb*n+1):(2*nb*n)] #N(0,22)

      X.id = rep(seq(1:nb), n)
      X.id = X.id[order(X.id)] #order the id
      dfExample = data.frame(X1 = X1, X2 = X2, Y1 = Y1,Y2 = Y2, X.id = X.id)
    }else if(fileName %in% c("example81", "example85") ){
      df = dfExample[,1:10]
      X.id = rep(seq(1:nrow(df)), ncol(df))
      X = as.numeric(unlist(df)) #convert dataset to vector
      X = X[order(X.id)] #order according the id
      X.id = X.id[order(X.id)] #order the id
      dfExample = data.frame(X = X, X.id = X.id)
    }else if(fileName == "example83"){
      df = dfExample[,1:6]
      X.id = rep(seq(1:nrow(df)), ncol(df))
      X = as.numeric(unlist(df)) #convert dataset to vector
      X = X[order(X.id)] #order according the id
      X.id = X.id[order(X.id)] #order the id
      dfExample = data.frame(X = X, X.id = X.id)
    }else if(fileName == "example91"){
      nb = 20 #number of batches
      n = 50 #numbers of observations in each batch
      df1 = dfExample$V1
      X1 = df1[1:(nb*n)] #N(0,1)
      X.id = rep(seq(1:nb), each=n)
      X1 = X1[order(X.id)] #order according the id

      df2 = dfExample$V2
      X2 = df2[1:(nb*n)] #N(0,1)
      X2 = X2[order(X.id)] #order according the id

      X.id = X.id[order(X.id)] #order the id
      dfExample = data.frame(X1 = X1, X2 = X2, X.id = X.id)
    }else if(fileName == "example93"){
      nb = 20 #number of batches
      n = 10 #numbers of observations in each batch
      df1 = dfExample$X
      X1 = df1[1:(nb*n)] #N(0,1)
      X.id = rep(seq(1:nb), each=n)
      X1 = X1[order(X.id)] #order according the id

      df2 = dfExample$X.1
      X2 = df2[1:(nb*n)] #N(0,1)
      X2 = X2[order(X.id)] #order according the id

      X.id = X.id[order(X.id)] #order the id
      dfExample = data.frame(X1 = X1, X2 = X2, X.id = X.id)
    }else if(fileName %in% c("example71", "example74a", "example74b", "example74c")){
      nb = 30 #number of batches
      n = 1 #numbers of observations in each batch
      df1 = dfExample$X
      X1 = df1[1:(nb*n)] #N(0,1)
      X.id = rep(seq(1:nb), each=n)
      X1 = X1[order(X.id)] #order according the id

      df2 = dfExample$X.1
      X2 = df2[1:(nb*n)] #N(0,1)
      X2 = X2[order(X.id)] #order according the id

      df3 = dfExample$X.2
      X3 = df3[1:(nb*n)] #N(0,1)
      X3 = X3[order(X.id)] #order according the id

      X.id = X.id[order(X.id)] #order the id
      dfExample = data.frame(X1 = X1, X2 = X2, X3 = X3, X.id = X.id)
    }
    write.csv(dfExample, src)
    assign(fileName, dfExample)
    do.call("use_data", list(as.name(fileName), overwrite = TRUE))
  }
}
