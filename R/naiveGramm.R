#minerva
naiveGramm<-function(x,y,z,r=0.5,alpha=0.05){
  #  colnames(y)<-colnames(x)<-colnames(z)
    options(warn = -1)
    xdata<-as.data.frame(x@assays@data@listData$counts)
    xname<-as.character(x@elementMetadata@listData$X)
    ydata<-as.data.frame(y@assays@data@listData$counts)
    yname<-as.character(y@elementMetadata@listData$X)
    if(is.data.frame(z)){

      zdata<-as.data.frame(z@assays@data@listData$counts)

    rresult<-presult<-ltresult<-lnresult<-matrix(0,nrow = nrow(xdata),
                                                 ncol = nrow(ydata))
    rownames(presult) <-rownames(rresult)<-rownames(lnresult)<-xname
    colnames(presult)<-colnames(rresult)<-colnames(lnresult)<-yname
    for(i in seq_len(nrow(xdata))){
    for(j in seq_len(nrow(ydata))){
        options(warn = -1)
        x1<-t(xdata[i,])
        y1<-t(rbind(ydata[j,],zdata))
        lmx<-lm(x1~y1)
        x11<-y1[,1]-lmx$coefficients[3]*y1[,2]
        x2<-t(rbind(xdata[i,],zdata))
        y2<-t(ydata[j,])
        lmy<-lm(y2~x2)
        y22<-x2[,1]-lmy$coefficients[3]*x2[,2]
        summ<-summary(lmx)
        pvalue<-summ$coefficients[2,4]
        rvalue<-lmx$coefficients[2]*sd(y1[,1])/sd(x1)
        r2value<-rvalue^2
        if(pvalue<alpha&rvalue>r){
        presult[i,j]<-pvalue
        rresult[i,j]<-rvalue
        ltresult[i,j]<-r2value
        lnresult[i,j]<-"linear"
        }
        else{
          #MIC

        micxy<-minerva::mine(y22,x11)
        micr<-micxy$MIC
        micp<-0
        for(ii in seq_len(101))
        {
            #bootstrap
            options(warn = -1)
            bootx<-matrix(sample(x11,replace = TRUE))
            booty<-matrix(sample(y22,replace = TRUE))
            if(sd(bootx)==0){
            bootx<-matrix(sample(x11,replace = FALSE))
            }
            if(sd(booty)==0){
            booty<-matrix(sample(y22,replace = FALSE))
            }
            tmp<-minerva::mine(booty,bootx)
            MICtp<-tmp$MIC
            tempM<-ifelse(micr<=MICtp,1,0)
            micp<-tempM+micp
            }
            micp<-micp/101
            prsmic<-cor.test(y22,x11,method="pearson")
            prsr<-prsmic$estimate
            ltmic<-1-micr+prsr^2
            presult[i,j]<-micp
            rresult[i,j]<-micr
          #ltresult[i,j]<-ltmic
            lnresult[i,j]<-"nonlinear"
        }

        }

    }
}
    if(!is.data.frame(z)){

    options(warn = -1)
rresult<-presult<-matrix(0,
                  nrow = nrow(xdata),ncol = nrow(ydata))
ltresult<-lnresult<-matrix(0,
                           nrow = nrow(xdata),ncol = nrow(ydata))
    rownames(presult) <-rownames(rresult)<-rownames(lnresult)<-xname
    colnames(presult)<-colnames(rresult)<-colnames(lnresult)<-yname

    for(i in seq_len(nrow(xdata))){
    for(j in seq_len(nrow(ydata))){
        x1<-t(xdata[i,])
        y1<-t(ydata[j,])
        lmx<-lm(x1~y1)
        summ<-summary(lmx)
        pvalue<-summ$coefficients[2,4]
        rvalue<-lmx$coefficients[2]*sd(y1)/sd(x1)
        if(pvalue<alpha&rvalue>r){
        presult[i,j]<-pvalue
        rresult[i,j]<-rvalue
        r2value<-rvalue^2
        ltresult[i,j]<-r2value
        lnresult[i,j]<-"linear"
        }
        else{
          #MIC
        warnings('off')
        micxy<-minerva::mine(y1,x1)
        micr<-micxy$MIC
        micp<-0
        for(ii in seq_len(101))
        {
            bootx<-matrix(sample(x1,replace = FALSE))
            booty<-matrix(sample(y1,replace = FALSE))
            tmp<-minerva::mine(bootx,booty)
            MICtp<-tmp$MIC
            tempM<-ifelse(micr<=MICtp,1,0)
            micp<-tempM+micp
            }
            micp<-micp/101
        prsmic<-cor.test(y1,x1,method="pearson")
        prsr<-prsmic$estimate
        ltmic<-1-micr+prsr^2
        presult[i,j]<-micp
        rresult[i,j]<-micr
        lnresult[i,j]<-"nonlinear"
        }
        }
        }
        }
    results2<-list()
    results2[["r"]]<-rresult
    results2[["p"]]<-presult
    results2[["type"]]<-lnresult
    return(results2)
}

