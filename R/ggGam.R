#'Extract independent variable name from formula
#'@param formula formula
formula2vars=function(formula){
        temp=deparse(formula)
        temp=paste0(temp,collapse="")

        temp1=unlist(strsplit(temp,"~"))[2]
        temp2=unlist(strsplit(temp1,"\\+"))
        temp3=gsub("^.*\\(|\\)| ","",temp2)
        gsub(",.*$","",temp3)
}

#' Make new data for predict
#' @param model An object
#' @param length numeric length of continuous variable to to predict
#' @param by character optional factor variable
#' @export
makeNewData=function(model,length=100,by=NULL){
     xvars=formula2vars(model$formula)
     xvars
     if(!is.null(by)) xvars=setdiff(xvars,by)
     xvars
     x=model$var.summary

     names(x)

     for(j in 1:length(xvars)){

          select=which(names(x)==xvars[j])
          select
          result=list()

          for(i in 1:length(x)){
               if(is.factor(x[[i]])) {
                    if(i==select) {
                         result[[i]]=levels(model$model[[names(x)[i]]])
                    } else if(names(x)[i] %in% by){
                         result[[i]]=levels(model$model[[names(x)[i]]])
                    } else{
                         result[[i]]=x[[i]][1]
                    }
               } else {  ## numeric
                    if(i==select){
                         result[[i]]=seq(from=x[[i]][1],to=x[[i]][3],length.out=length)
                    } else{
                         result[[i]]=x[[i]][2]
                    }
               }
          }
          result
          names(result)=names(x)
          newdata=as.data.frame(expand.grid(result))
          predict(model,newdata=newdata,type="response",se.fit=TRUE,na.action="na.omit")
          df1=as.data.frame(predict(model,newdata=newdata,type="response",se.fit=TRUE))
          df1$ymax=df1$fit+1.96*df1$se.fit
          df1$ymin=df1$fit-1.96*df1$se.fit
          df=cbind(newdata,df1)
          df$xvar=names(x)[select]
          if(j==1) {
               final=df
          } else{
               final=rbind(final,df)
          }
     }
     final
}

#' Draw a ggplot with an object of class gam
#' @param model An object of class gam
#' @param select numeric Choices of dependent variables to plot
#' @param point logical Whether or not draw point
#' @param se logical Whether or not draw confidence interval
#' @param by character optional name of factor variable
#' @param byauto logical Whether or not select categorical variables automatically
#' @param scales Should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y")?
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line geom_ribbon ylab facet_wrap xlab
#' @examples
#' require(mgcv)
#' mtcars$am=factor(mtcars$am,labels=c("automatic","manual"))
#' model <- gam(mpg ~ s(wt,by=am)+am, data = mtcars, method = "REML")
#' plot(model,shift=coef(model)[1],pages=1,all.terms=TRUE,shade=TRUE,seWithMean=TRUE,residuals=TRUE)
#' ggGam(model)
#' ggGam(model,point=TRUE)
ggGam=function(model,select=NULL,point=FALSE,se=TRUE,by=NULL,byauto=TRUE,scales="free_x"){

     # select=NULL;point=FALSE;se=TRUE;by=NULL;byauto=FALSE;scales="free";

     byall=names(model$var.summary)[sapply(model$var.summary,is.factor)]
     if(length(byall)>1) {
          if(byauto & is.null(by)) {
                  cat("Only one factor variable is used as fill variable")
                  by=byall[1]
          }
     } else if(length(byall)==1){
             if(byauto & is.null(by)) by=byall
     } else{
          by=NULL
     }

     df1=makeNewData(model,by=by)
     xvars=formula2vars(model$formula)
     xvars2=setdiff(xvars,byall)
     yvar=names(model$model)[1]
     xvars2
     df1
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
     longdf<-pivot_longer(model$model,cols=all_of(xvars2))
     longdf$name=factor(longdf$name,levels=xvars2)
     if(point) {
          p=ggplot(data=longdf,aes_string(x="value",fill=fillvar,group=fillvar))+
               geom_point(aes_string(y=yvar,color=fillvar),alpha=0.3)
     } else{
          if(!is.null(fillvar)){
          p<-ggplot(df3,aes_string(x="value",fill=fillvar,group=fillvar))
          } else{
          p<-ggplot(df3,aes_string(x="value",group="1"))
          }
     }

     p <- p+geom_line(data=df3,aes_string(y="fit"),color="blue")
     if(se) p<-p+
          geom_ribbon(data=df3,aes_string(y="fit",ymax="ymax",ymin="ymin"),alpha=0.3)

     p<-p+ylab(yvar)
     if(length(xvars2)>1) {
             p<-p+facet_wrap("name",scales=scales)+xlab("")
     } else{
             p<-p+xlab(xvars2)
     }
     p
}



#' Draw a ggplot summarizing categorical variables with an object of class gam
#' @param model An object of class gam
#' @param scales Should scales be fixed ("fixed"), free ("free"), or free in one dimension ("free_x", "free_y")?
#' @importFrom ggplot2 ggplot geom_errorbar scale_x_discrete labs guide_axis facet_wrap
#' @importFrom dplyr vars
#' @export
#' @examples
#' require(mgcv)
#' data(Wage,package="ISLR")
#' model=gam(wage~s(age,by=education)+education+jobclass,data=Wage,method="REML")
#' plot(model,shift=coef(model)[1],pages=1,all.terms=TRUE,shade=TRUE,seWithMean=TRUE,residuals=TRUE)
#' ggGamCat(model)
ggGamCat=function(model,scales="free_x"){
        # scales="free_x"
        # require(ggplot2)
        df=makeNewData(model)
        catVars=names(model$var.summary)[sapply(model$var.summary,is.factor)]
        df=df[df$xvar %in% catVars,]
        df2<-tidyr::pivot_longer(df,cols=all_of(catVars))
        df2
        df2<-df2[df2[["xvar"]]==df2[["name"]],]
        yvar=names(model$model)[1]
        df2
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
        lines(yhat~newx,col="red",lwd=2)
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
       lines(yhat~newx,col="red",lwd=2)
}

}
