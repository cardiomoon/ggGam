#'Extract independent variable name from formula
#'@param formula formula
#'@export
formula2vars=function(formula){
        #formula=model$formula
        temp=deparse(formula)

        temp=paste0(temp,collapse="")
        temp
        temp1=unlist(strsplit(temp,"~"))[2]
        temp2=unlist(strsplit(temp1,"\\+"))
        temp2=gsub(" ","",temp2)
        temp2
        temp3=gsub(",bs=\".*\"|^s\\(|^ti\\(|\\)$","",temp2)
        temp3
        temp4=gsub("by=|,k=.*$|,sp=.*$","",temp3)
        unique(unlist(strsplit(temp4,",")))

}

#' Make new data for predict
#' @param model An object
#' @param length numeric length of continuous variable to to predict
#' @param by character optional factor variable
#' @param type character type argument to be passed to predict.gam
#' @importFrom stats plogis
#' @importFrom predict3d restoreData2 restoreData3
#' @export
makeNewData=function(model,length=100,by=NULL,type="response"){
                 # length=100;by=NULL;type="response"
     xvars=formula2vars(model$formula)
     xvars
     if(!is.null(by)) {
             xvars2=setdiff(xvars,by)
     } else{
        xvars2=xvars
     }
     xvars2

     for(j in 1:length(xvars)){
        result=list()
        for(i in 1:length(xvars)){
               x=model$model[[xvars[i]]]
               if(is.numeric(x)) {  ## numeric
                 if(i==j){
                   result[[i]]=seq(from=min(x,na.rm=TRUE),to=max(x,na.rm=TRUE),length.out=length)
                 } else{
                   result[[i]]=mean(x,na.rm=TRUE)
                 }
               } else if(is.factor(x)) {
                    if(i==j) {
                         result[[i]]=levels(x)
                    } else if(xvars[i] %in% by){
                         result[[i]]=levels(x)
                    } else{
                         result[[i]]=levels(x)[1]
                    }
               } else {
                   if(i==j) {
                     result[[i]]=unique(x)
                   } else if(xvars[i] %in% by){
                     result[[i]]=unique(x)
                   } else{
                     result[[i]]=names(which.max(table(x)))[1]
                 }

               }
          }
          result
          names(result)=xvars
          newdata=as.data.frame(expand.grid(result))
          newdata=restoreData2(newdata)
          newdata=restoreData3(newdata)
          names(newdata)
          if(model$family$family=="binomial"){
                  df1=as.data.frame(predict(model,newdata=newdata,type="link",se.fit=TRUE))
                  df1$ymax=df1$fit+1.96*df1$se.fit
                  df1$ymin=df1$fit-1.96*df1$se.fit
                  if(type=="response") df1[]=lapply(df1,plogis)
          } else if(model$family$family=="Cox PH"){
                  df1=as.data.frame(predict(model,newdata=newdata,type="link",se.fit=TRUE))
                  df1$ymax=df1$fit+1.96*df1$se.fit
                  df1$ymin=df1$fit-1.96*df1$se.fit
                  if(type=="response") df1[]=lapply(df1,plogis)
          }else{

            df1=as.data.frame(predict(model,newdata=newdata,type=type,se.fit=TRUE))
            df1$ymax=df1$fit+1.96*df1$se.fit
            df1$ymin=df1$fit-1.96*df1$se.fit
          }
          df=cbind(newdata,df1)
          df$xvar=xvars[j]
          if(j==1) {
               final=df
          } else{
               final=rbind(final,df)
          }
     }
     final
}


#'Make list of args with call string
#'@param string A character
#'@return A list
call2vars=function(string){
  temp1 = gsub("^.*\\(|\\)$| |\"", "", string)
  temp2=unlist(strsplit(temp1,","))
  temp3=unlist(strsplit(temp2,"="))
  result=list()
  for(i in 1:(length(temp3)/2)){
    result[[temp3[2*i-1]]]=temp3[2*i]
  }
  result
}

#' Draw a ggplot with an object of class gam
#' @param model An object of class gam
#' @param select numeric Choices of dependent variables to plot
#' @param point logical Whether or not draw point
#' @param se logical Whether or not draw confidence interval
#' @param by NULL or character optional name of factor variable
#' @param scales Should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y")?
#' @param type character type argument to be passed to predict.gam
#' @param byauto logical Whether or not choose variabels to facet automatically
#' @param facet logical Whether or not make facetted plot
#' @param fillcolor Character Name of fillcolor
#' @param pointalpha,fillalpha Numeric alpha value
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line geom_ribbon ylab facet_wrap xlab scale_x_continuous
#' @examples
#' require(mgcv)
#' model <- gam(mpg ~ s(wt), data = mtcars, method = "REML")
#' plot(model,shift=coef(model)[1],pages=1,all.terms=TRUE,shade=TRUE,seWithMean=TRUE,residuals=TRUE)
#' ggGam(model)
#' mtcars$am=factor(mtcars$am,labels=c("automatic","manual"))
#' model <- gam(mpg ~ s(wt)+am, data = mtcars, method = "REML")
#' plot(model,shift=coef(model)[1],pages=1,all.terms=TRUE,shade=TRUE,seWithMean=TRUE,residuals=TRUE)
#' ggGam(model)
#' ggGam(model,by=am)
#' ggGam(model,by=am,facet=TRUE)
#' model <- gam(mpg ~ s(wt,by=am)+ am, data = mtcars, method = "REML")
#' plot(model,shift=coef(model)[1],pages=1,all.terms=TRUE,shade=TRUE,seWithMean=TRUE,residuals=TRUE)
#' ggGam(model)
#' ggGam(model,by=am)
#' ggGam(model,by=am,facet=TRUE)
#' ggGam(model,by=am,point=FALSE)
#' ggGamCat(model)
#' data(mpg,package="gamair")
#' model <- gam(hw.mpg ~ s(weight, fuel, bs = "fs"),data = mpg,method = "REML")
#' ggGam(model)
#' ggGam(model,by=fuel)
#' ggGam(model,by=fuel,facet=TRUE)
#' model2 <- gam(hw.mpg ~ s(weight) + s(length) + s(price) + fuel + drive + style,
#'    data=mpg, method="REML")
#' plot(model2,shift=coef(model)[1],pages=1,all.terms=TRUE,shade=TRUE,seWithMean=TRUE,residuals=TRUE)
#' ggGam(model2,se=TRUE)
#' ggGam(model2,se=TRUE,by=fuel)
#' ggGam(model2,se=TRUE,by=fuel,select=1)
#' ggGam(model2,se=FALSE,by=style)
#' ggGam(model2,se=FALSE,by=drive,point=FALSE)
#' model6 <- gam(hw.mpg ~ s(weight, by=fuel), data=mpg, method="REML")
#' ggGam(model6)
ggGam=function(model,select=NULL,point=TRUE,se=TRUE,by=NULL,scales="free_x",type=NULL,byauto=FALSE,facet=FALSE,
               fillcolor="red",pointalpha=0.3,fillalpha=0.3){

       # model=m1;
      # select=NULL;by=NULL;point=TRUE;se=TRUE;scales="free_x";type=NULL;byauto=FALSE;facet=FALSE;
      # pointalpha=0.3;fillalpha=0.3

     temp=deparse(match.call())
     res=call2vars(temp)
     by=res$by

     byall=names(model$var.summary)[sapply(model$var.summary,is.factor)]
     if(length(byall)==1) {
       if(byauto && is.null(by)) {
         by=byall[1]
       }
     }

     if(is.null(type)) {
        if(model$family$family %in% c("Cox PH","binomial")) {
          type="link"
        } else{
          type="response"
        }

     }
     df1=makeNewData(model,by=by,type=type)

     table(df1$xvar)
     xvars=formula2vars(model$formula)
     xvars
     byall
     xvars2=setdiff(xvars,byall)
     yvar=names(model$model)[1]
     xvars2
     df1
     select
     if(!is.null(select)) {
          xvars2=xvars2[select]
          if(is.na(xvars2)) xvars2=xvars[select]
          df1=df1[df1$xvar %in% xvars2,]
     }
     df2<-tidyr::pivot_longer(df1,cols=all_of(xvars2))
     df3<-df2[df2[["xvar"]]==df2[["name"]],]
     df3$name=factor(df3$name,levels=xvars2)
     if(is.null(by)) {
          fillvar=NULL
     } else {
          fillvar=by[1]
     }
     model$model
     xvars2
     fillvar
     longdf<-pivot_longer(model$model,cols=all_of(xvars2))
     longdf$name=factor(longdf$name,levels=xvars2)
     longdf
     if(point) {
          p=ggplot(data=longdf,aes_string(x="value",fill=fillvar,group=fillvar))+
               geom_point(aes_string(y=yvar,color=fillvar),alpha=pointalpha)
     } else{
          if(!is.null(fillvar)){
          p<-ggplot(df3,aes_string(x="value",fill=fillvar,group=fillvar))
          } else{
          p<-ggplot(df3,aes_string(x="value",group="1"))
          }
     }
     p<-p+ geom_line(data=df3,aes_string(y="fit",color=fillvar))


     if(se) {
         if(is.null(fillvar)) {
           p<- p+geom_ribbon(data=df3,aes_string(y="fit",ymax="ymax",ymin="ymin"),fill=fillcolor,alpha=fillalpha)
         } else{
          p<- p+geom_ribbon(data=df3,aes_string(y="fit",ymax="ymax",ymin="ymin"),alpha=fillalpha)
         }
     }

     if(model$family$family=="Cox PH") {
         yvar=ifelse(type=="link","Hazard Ratio","Probabilty")
     } else if(model$family$family=="binomial"){
         yvar=ifelse(type=="link","Odds Ratio","Probability")
     }

     p<-p+ylab(yvar)+scale_x_continuous(guide=guide_axis(n.dodge=2))
     p
     facet
     select
     fillvar
     xvars2
     if(length(xvars2)>1) {
             p<-p+facet_wrap("name",scales=scales)+xlab("")+
               theme(legend.position="top")
     } else if(facet){
           if(is.null(select)){
             if(is.null(fillvar)){
               p<-p+xlab(xvars2)
             } else{
            p<-p+facet_wrap(fillvar,scales=scales)+xlab("")+
               theme(legend.position="none")
             }
           } else if(length(select)>1){
             p<-p+facet_wrap(fillvar,scales=scales)+xlab("")+
               theme(legend.position="none")
           } else{
             p<-p+xlab(xvars2)
           }
     } else if(!is.null(fillvar)){
            p<-p+theme(legend.position="top")+xlab(xvars2)
      }else{
             p<-p+xlab(xvars2)
     }
     p
}



#' Draw a ggplot summarizing categorical variables with an object of class gam
#' @param model An object of class gam
#' @param select numeric Choices of dependent categorical variables to plot
#' @param scales Should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y")?
#' @param type character type argument to be passed to predict.gam
#' @importFrom ggplot2 ggplot geom_errorbar scale_x_discrete labs guide_axis facet_wrap theme
#' @importFrom dplyr vars
#' @export
#' @examples
#' require(mgcv)
#' data(Wage,package="ISLR")
#' model=gam(wage~s(age,by=education)+education+jobclass,data=Wage,method="REML")
#' plot(model,shift=coef(model)[1],pages=1,all.terms=TRUE,shade=TRUE,
#'     seWithMean=TRUE,residuals=TRUE)
#' ggGamCat(model)
ggGamCat=function(model,select=NULL,scales="free_x",type=NULL){
        # scales="free_x";select=c(2,3)

      if(is.null(type)) {
          if(model$family$family %in% c("Cox PH","binomial")) {
            type="link"
          } else{
            type="response"
          }

      }
        df=makeNewData(model,type=type)
        catVars=names(model$var.summary)[sapply(model$var.summary,is.factor)]
        catVars
        if(!is.null(select)) catVars=catVars[select]
        df=df[df$xvar %in% catVars,]
        df2<-tidyr::pivot_longer(df,cols=all_of(catVars))
        df2
        df2<-df2[df2[["xvar"]]==df2[["name"]],]
        yvar=names(model$model)[1]
        df2

        if(model$family$family=="Cox PH") {
          yvar=ifelse(type=="link","Hazard Ratio","Probabilty")
        } else if(model$family$family=="binomial"){
          yvar=ifelse(type=="link","Odds Ratio","Probability")
        }
        df2$name=factor(df2$name,levels=catVars)
        ggplot(df2, aes_string(x="value",y="fit")) +
                geom_point() +
                geom_errorbar(aes_string(ymin="fit-2*se.fit",ymax="fit+2*se.fit"),width=0.25) +
                scale_x_discrete(guide=guide_axis(n.dodge=2))+
                facet_wrap("name",scales=scales)+
                labs(x="",y=yvar)
}



#'Show basis function
#' @param model An object of class gam
#' @param which NULL or one of 1:4
#' @param ... Further arguments to be passed to plot
#' @importFrom graphics lines matplot par plot
#' @importFrom stats predict
#' @export
#' @examples
#' require(mgcv)
#' set.seed(1)
#' x <- seq(0, pi * 2, 0.1)
#' sin_x <- sin(x)
#' y <- sin_x + rnorm(n = length(x), mean = 0, sd = sd(sin_x / 2))
#' data1 <- data.frame(y,x)
#' model=gam(y~s(x),method="REML")
#' basisFun(model)
#' basisFun(model,which=1)
#' basisFun(model,which=2)
#' basisFun(model,which=3)
#' basisFun(model,which=4)
basisFun=function(model,which=NULL,...){


y=model$model[[1]]
x=model$model[[2]]
model_matrix=predict(model,type="lpmatrix")

result=model_matrix
for(i in 1:nrow(result)){
        for(j in 1:ncol(result)){
                result[i,j]=model_matrix[i,j]*model$coeff[j]
        }
}
result
res1=apply(result,1,sum)
if(is.null(which)) {
        par(mfrow=c(2,2))
        plot(y~x,main="Raw Data",pch=21,bg='gray',col=NA,...)
        plot(y~x,main="Raw Data + Basis Functions",pch=21,bg='gray',col=NA,...)
        matplot(x,model_matrix,add=T,type="l",...)
        plot(y~x,main="Basis Functions*Coeff",pch=21,bg='gray',col=NA,...)
        matplot(x,result,add=T,type="l",...)
        plot(y~x,main="Basis Functions*Coeff with Prediction",pch=21,bg='gray',col=NA,...)
        matplot(x,result,add=T,type="l",...)
        newx=seq(min(x),max(x),length.out=100)
        newdata=data.frame(newx)
        names(newdata)=names(model$model)[2]
        yhat=predict(model,newdata=newdata)
        lines(yhat~newx,lwd=2)
        par(mfrow=c(1,1))

} else if(which==1) {
        plot(y~x,main="Raw Data",pch=21,bg='gray',col=NA,...)
} else if(which==2){
        plot(y~x,main="Raw Data + Basis Functions",pch=21,bg='gray',col=NA,...)
        matplot(x,model_matrix,add=T,type="l",...)
} else if(which==3){
        plot(y~x,main="Basis Functions*Coeff",pch=21,bg='gray',col=NA,...)
        matplot(x,result,add=T,type="l",...)
} else{
       plot(y~x,main="Basis Functions*Coeff with Prediction",pch=21,bg='gray',col=NA,...)
       matplot(x,result,add=T,type="l",...)
       newx=seq(min(x),max(x),length.out=100)
       newdata=data.frame(newx)
       names(newdata)=names(model$model)[2]
       yhat=predict(model,newdata=newdata)
       lines(yhat~newx,lwd=2)
}

}
