#' Draw Kaplan-Meier Survival Plot with an object of class 'gam'
#' @param model a fitted gam object as produced by gam()
#' @param data A data frame or list containing the values of the model covariates at which predictions are required
#' @param np Numeric Number of time intervals
#' @param timevar character Name of variable containing time
#' @param until Maximum value of time
#' @param id Optional labels of preset values
#' @return A ggplot
#' @export
#' @examples
#' require(mgcv)
#' data(pbc,package="survival")
#' pbc$status1 <- as.numeric(pbc$status==2)
#' pbc$stage <- factor(pbc$stage)
#' model = gam(time ~ trt +sex + s(sqrt(protime))+s(platelet)+ s(age)+s(bili)+s(albumin),
#'    weights=status1, family=cox.ph, data=pbc, method="REML")
#' dfSex=averageData(pbc,list(sex=c("m","f")))
#' drawSurv(model,data=dfSex,id=list(sex=c("m","f")))
#' dfBili=averageData(pbc,list(bili=c(1,10)))
#' drawSurv(model,data=dfBili,id=list(bili=c(1,10)))
drawSurv=function(model,data,np=100,timevar="time",until=NULL,id=list()){

    if(is.null(until)) until=max(model$model[[timevar]],na.rm=TRUE)

    if(length(id)==0) id=list(id=1:nrow(data))

    for( i in 1:nrow(data)){
        newd <- data.frame(matrix(0,np,0))
        for (n in names(data)) newd[[n]] <- rep(data[[n]][i],np)
        newd$time <- seq(0,until,length=np)
        fv <- predict(model,newdata=newd,type="response",se=TRUE)
        newd$fit=fv$fit
        # newd$ymax=fv$fit+se*fv$se.fit
        # newd$ymin=fv$fit-se*fv$se.fit
        se <- fv$se.fit/fv$fit
        newd$ymax=exp(log(fv$fit)+se)
        newd$ymin=exp(log(fv$fit)-se)

        idname=names(id)[1]
        newd[[idname]]=id[[1]][i]
        if(i==1){
            final=newd
        } else{
            final=rbind(final,newd)
        }
    }
    final[[idname]]=factor(final[[idname]])
    final
    ggplot(data=final,aes_string(x="time",y="fit",fill=idname,group=idname))+
        geom_line(aes_string(color=idname))+
        geom_ribbon(aes_string(ymax="ymax",ymin="ymin"),alpha=0.3)+
        ylab("cumulative survival")+xlab("days")+
        theme_bw()+
        theme(legend.position = "top")
}

#' Make average data from a data.frame
#' @param data A data.frame
#' @param newValue Optional list of preset values
#' @return A data.frame
#' @export
#' @examples
#' data(pbc,package="survival")
#' dfSex=averageData(pbc,list(sex=c("m","f")))
#' dfBili=averageData(pbc,list(bili=c(1,10)))
averageData=function(data,newValue=list()){
    newd=list()
    for(i in 1:ncol(data)){
        if(is.numeric(data[[i]])) {
            newd[[i]]=mean(data[[i]],na.rm=TRUE)
        } else if(is.factor(data[[i]])){
            newd[[i]]=levels(data[[i]])[1]
        } else{
            newd[[i]]=sort(unique(data[[i]]))[1]
        }
    }
    names(newd)=names(data)
    df=as.data.frame(newd)
    df
    if(length(newValue)>0){
        no=length(newValue[[1]])
        for(i in 1:no){
            if(i==1) {
                final=df
            } else{
                final=rbind(final,df)
            }
        }
        final[[names(newValue)[1]]]=newValue[[1]]
        df=final
    }
    df
}
