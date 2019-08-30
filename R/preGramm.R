# psych
# new DMwR phyloseq
preGramm<-function(A,B,metaNor=TRUE,rarefaction=FALSE){
    Adata<-as.data.frame(A@assays@data@listData$counts)#metabolites
    Aname<-as.data.frame(A@elementMetadata@listData$X)
    Bdata<-as.data.frame(B@assays@data@listData$counts)#microbe
    Bname<-as.data.frame(B@elementMetadata@listData$X)

    if(rarefaction){
      Bdata_number<-nrow(Bdata)*ncol(Bdata)
      taxmat = matrix(sample(letters, Bdata_number, replace = TRUE), nrow = nrow(Bdata),
                      ncol = 7)
      rownames(taxmat) <- rownames(Bdata)
      colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family",
                            "Genus", "Species")
      OTU_B = phyloseq::otu_table(Bdata, taxa_are_rows = TRUE)
      TAX_B = phyloseq::tax_table(taxmat)
      physeq_B = phyloseq(OTU_B, TAX_B)
     # set.seed(100)
      raretest<-phyloseq::rarefy_even_depth(physeq_B,rngseed=100)
      Bdata<-raretest@otu_table@.Data



    }
    else{
      Bdata<-Bdata
    }


    Bdata[Bdata==0]<-1  #microbiome
    Adata_n<-Adata
    Adata[Adata==0]<-NA
    Adata<-DMwR::knnImputation(Adata,k=3,scale = TRUE)
    if(metaNor){

    # sumA<-apply(Adata,2,sum)
    #
    #
    # for(i in seq_len(nrow(Adata))){
    # for(j in seq_len(ncol(Adata))){
    # Adata_n[i,j]<-Adata[i,j]/sumA[j]*30000
    # }
    # }

      Adata_n<- apply(Adata,2,function(each_col){
      the_sum  <-  sum(each_col,na.rm = TRUE)
      each_col <-  each_col/the_sum*30000
       return(each_col)
    })



}
    else{
    Adata_n<-Adata_n
    }
    Bdata_n<-Bdata
#     sumB<-apply(Bdata,2,sum)
#     for(i in seq_len(nrow(Bdata))){
#     for(j in seq_len(ncol(Bdata))){
#     Bdata_n[i,j]<-Bdata[i,j]/sumB[j]*30000
#     }
# }
    Bdata_n<- apply(Bdata,2,function(each_col){
      the_sum  <-  sum(each_col,na.rm = TRUE)
      each_col <-  each_col/the_sum*30000
      return(each_col)
    })



    Alog<-log(Adata_n+1)
#     gB<-psych::geometric.mean(Bdata_n)
#     Blog<-Bdata_n
#
#     for(i in seq_len(nrow(Bdata_n))){
#     for(j in seq_len(ncol(Bdata_n))){
#     Blog[i,j]<-log(Bdata_n[i,j])/log(gB[j])
#     }
# }
    Blog<-apply(Bdata_n,2,function(each_col){
      log_col<-psych::geometric.mean(each_col)
      each_col<-log(each_col)/log(log_col)

    })



    Alog<-cbind(Aname,Alog)
    Blog<-cbind(Bname,Blog)
    results1<-list()
    results1[["x"]]<-Alog
    results1[["y"]]<-Blog
    return(results1)
}

