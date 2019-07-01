# psych
preGramm<-function(A,B,metaNor=TRUE){
    Adata<-as.data.frame(A@assays@data@listData$counts)
    Aname<-as.data.frame(A@elementMetadata@listData$X)
    Bdata<-as.data.frame(B@assays@data@listData$counts)
    Bname<-as.data.frame(B@elementMetadata@listData$X)
    Bdata[Bdata==0]<-1  #microbiome
    Adata_n<-Adata
    if(metaNor){
    sumA<-apply(Adata,2,sum)
    for(i in seq_len(nrow(Adata))){
    for(j in seq_len(ncol(Adata))){
    Adata_n[i,j]<-Adata[i,j]/sumA[j]*30000
    }
    }
}
    else{
    Adata_n<-Adata_n
    }
    Bdata_n<-Bdata
    sumB<-apply(Bdata,2,sum)
    for(i in seq_len(nrow(Bdata))){
    for(j in seq_len(ncol(Bdata))){
    Bdata_n[i,j]<-Bdata[i,j]/sumB[j]*30000
    }
}
    Alog<-log(Adata_n+1)
    gB<-psych::geometric.mean(Bdata_n)
    Blog<-Bdata_n
    for(i in seq_len(nrow(Bdata_n))){
    for(j in seq_len(ncol(Bdata_n))){
    Blog[i,j]<-log(Bdata_n[i,j])/log(gB[j])
    }
}
    Alog<-cbind(Aname,Alog)
    Blog<-cbind(Bname,Blog)
    results1<-list()
    results1[["x"]]<-Alog
    results1[["y"]]<-Blog
    return(results1)
}

