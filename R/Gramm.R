Gramm<-function(A,B,C,metaNor=TRUE,r=0.5,alpha=0.05){
  xdatag<-as.data.frame(A@assays@data@listData$counts)
  xnameg<-as.character(A@elementMetadata@listData$X)
  AAA<-cbind(xnameg,xdatag)

  ydatag<-as.data.frame(B@assays@data@listData$counts)
  ynameg<-as.character(B@elementMetadata@listData$X)
  BBB<-cbind(ynameg,ydatag)

    preresult<-preGramm(A,B,metaNor)
    xpos<-preresult$x
    ypos<-preresult$y
    nlfitGramstop<-function(X,Y){
    Xdata<-X[,-1]
    Ydata<-Y[,-1]
    Xname<-X[,1]
    Yname<-Y[,1]
    modellist<-c("line2P","line3P","log2P","exp2P","exp3P","power2P","power3P")
    for(i in seq_len(nrow(Xdata))){
    for(j in seq_len(nrow(Ydata))){
        xx<-as.numeric(Xdata[i,])
        yy<-as.numeric(Ydata[j,])
        yadj <- yy - min(yy) + 1
        zzz <- data.frame(xx, yadj)
        n = length(xx)
        k = 3
        r2<--2
        tryCatch({
        fit <- nls(yadj ~ SSlogis(xx,a,b,c),data=zzz)
        sum.exp3P <- summary(fit)
        ss.res <- sum((residuals(fit))^2)
        ss.total.uncor <- sum(yy^2)
        ss.total.cor <- sum((yy - mean(yy))^2)
        ss.reg <- ss.total.cor - ss.res
        dfR = k - 1
        dfE = n - k
        Fval = (ss.reg/dfR)/(ss.res/dfE)
        pval = pf(Fval, dfR, dfE, lower.tail = FALSE)
        pval <- unname(pval)
        RSE <- sum.exp3P$sigma
        SSE <- (RSE^2) * (n - 1)
        adjr2 <- 1 - SSE/((var(yy)) * (n - 1))
        r2 <- 1 - (1 - adjr2) * ((n - k)/(n - 1))
        r2 = format(r2, digits = 5)
        r2 = as.numeric(r2)
        },error=function(e){cat("S curve ERROR :",conditionMessage(e),"\n")})

        if(r2<0){
        fitr2<-matrix(NA,2,7)
        for(m in seq_len(7)){
        tryCatch({
        t1<-basicTrendline::trendline_summary(as.numeric(Xdata[i,]),
                      as.numeric(Ydata[j,]),model =modellist[m] )
        fitr2[1,m]<-t1$R.squared
        fitr2[2,m]<-m
        },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
        }
        maxr2<-which.max(fitr2[1,])
        tryCatch(
            {
            basicTrendline::trendline(as.numeric(Xdata[i,]),as.numeric(Ydata[j,]),
            model = modellist[maxr2],xlab = Xname[i],ylab = Yname[j])

            },error=function(e){cat("ERROR :",conditionMessage(e),"\n") }
        )
        }
        else{
        fitr2<-matrix(NA,2,7)
        for(m in seq_len(7)){
        tryCatch({
        t1<-basicTrendline::trendline_summary(as.numeric(Xdata[i,]),
                             as.numeric(Ydata[j,]),model =modellist[m] )
        fitr2[1,m]<-t1$R.squared
        fitr2[2,m]<-m
        },error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
        }
          #contrast s with 7 curves r2
        if(r2>max(fitr2[1,],na.rm = TRUE)){
        param<- vector("expression", 2)
        aa<-coef(fit)[1]
        aa<-unname(aa)
            aa<-format(aa,digits = 5)
            aa<-as.character(aa)
            bb<-coef(fit)[2]
            bb<-unname(bb)
            bb<-format(bb,digits = 5)
            bb<-as.character(bb)
            cc<-coef(fit)[3]
            cc<-unname(cc)
            cc<-format(cc,digits = 5)
            cc<-as.character(cc)
            pval<-format(pval, digits = 5)
            pval <- paste("=", unname(pval))
            rval <- paste("=", unname(r2))
            expression(italic("y") ==frac(aa,1+ ~ italic("e")^(frac(bb~-x,cc))))
param[1]<-bquote(expression(italic("y") ==frac(.(aa),
                            1+ ~ italic("e")^(frac(.(bb)~-x,.(cc))))))[2]
param[2]  <- bquote(expression(italic("R")^2 ==.(r2)* "," ~ ~italic("p") ~ ~.(pval) ))[2]
            #substitute(expression(italic("R")^2 == r2 * "," ~ ~italic("p") ~ ~pval ))[2]

            investr::plotFit(fit,interval = "confidence", shade = TRUE,
                             col.fit = "blue",xlab = Xname[i],ylab = Yname[j])

            legend("topleft",  legend = param,  cex = 1, bty = "n")

          }
            else{
            maxr2<-which.max(fitr2[1,])
            tryCatch(
              {
             basicTrendline::trendline(as.numeric(Xdata[i,]),as.numeric(Ydata[j,]),
             model = modellist[maxr2],xlab = Xname[i],ylab = Yname[j])

             },error=function(e){cat("ERROR :",conditionMessage(e),"\n") }
            )
          }
        }
      }
    }




  }
    #minerva
    naiveGramm1<-function(x,y,z,r=0.5,alpha=0.05){
      colnames(y)<-colnames(x)<-colnames(z)
      options(warn = -1)
      if(is.data.frame(z)){

        xdata<-x[,-1]
        ydata<-y[,-1]
        xname<-x[,1]
        yname<-y[,1]
        zdata<-z[,-1]
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
        xdata<-x[,-1]
        ydata<-y[,-1]
        xname<-x[,1]
        yname<-y[,1]
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


    uniresult<-naiveGramm1(xpos,ypos,C,r,alpha)
    result_score_loading<-list()
    unir<-uniresult$r
    max10<-arrayInd(order(unir,decreasing=TRUE)[seq_len(10)],dim(unir))
    pdf("R value top 10 pairs.pdf")
    for(top in seq_len(10)){
    xtop<-AAA[max10[top,1],]
    ytop<-BBB[max10[top,2],]
    nlfitGramstop(xtop,ytop)
}
    dev.off()
    resultsall<-list()
    resultsall[["pretreatment"]]<-preresult
    resultsall[["correlation"]]<-uniresult
}
