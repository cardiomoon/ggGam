#' Producing Poisson model data frame
#' @param x,t Numeric vector giving the coordinates of the points to be interpolated.
#' @param to an optional set of numeric values specifying where interpolation is to take place
#' @export
#' @importFrom stats approx
app <- function(x,t,to) {
    ## wrapper to approx for calling from apply...
    y <- if (sum(!is.na(x))<1) {
        rep(NA,length(to))
    } else {
        approx(t,x,to,method="constant",rule=2)$y
    }
    if (is.factor(x)) {
        factor(levels(x)[y],levels=levels(x)) }
    else {
        y
    }
} ## app



#' Convert data.frame into a data.frame of artificial Poisson data with the interpolated Time varying covariate values
#' @param dat A data.frame
#' @param event Name for Poisson response
#' @param et Name of event time
#' @param t Name of observation time
#' @param status Numeric 1 for death 0 otherwise
#' @param id Patient id
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @examples
#' library(survival)
#' pbcseq$status1 <- as.numeric(pbcseq$status==2) ## deaths
#' pb <- tdpois(pbcseq) ## conversion
#' @return A data.frame of artificial Poisson data
tdpois <- function(dat,event="z",et="futime",t="day",status="status1",
                   id="id") {


    if (event %in% names(dat)) warning("event name in use")

    te <- sort(unique(dat[[et]][dat[[status]]==1])) ## event times
    sid <- unique(dat[[id]])
    inter <- interactive()
    if (inter) prg <- txtProgressBar(min = 0, max = length(sid), initial = 0,
                                     char = "=",width = NA, title="Progress", style = 3)
    ## create dataframe for poisson model data
    dat[[event]] <- 0; start <- 1
    dap <- dat[rep(1:length(sid),length(te)),]
    for (i in 1:length(sid)) { ## work through patients
        di <- dat[dat[[id]]==sid[i],] ## ith patient's data
        tr <- te[te <= di[[et]][1]] ## times required for this patient
        ## Now do the interpolation of covariates to event times...
        um <- data.frame(lapply(X=di,FUN=app,t=di[[t]],to=tr))
        ## Mark the actual event...
        if (um[[et]][1]==max(tr)&&um[[status]][1]==1) um[[event]][nrow(um)] <- 1
        um[[et]] <- tr ## reset time to relevant event times
        dap[start:(start-1+nrow(um)),] <- um ## copy to dap
        start <- start + nrow(um)
        if (inter) setTxtProgressBar(prg, i)
    }
    if (inter) close(prg)
    dap[1:(start-1),]
}



#' Extract x and coordinate with stepfun
#' @param x An object of class stepfun
#' @return A data.frame
#' @importFrom stats knots median
#' @export
stepfun2df=function(x)
{
    knF <- knots(x)
    xval <- knF
    rx <- range(xval)
    dr <- if (length(xval) > 1L) {
            max(0.08 * diff(rx), median(diff(xval)))
         } else {
             abs(xval)/16
         }
    xlim <- rx + dr * c(-1, 1)
        dr <- if (length(xval) > 1L) {
        max(0.08 * diff(rx), median(diff(xval)))
    }  else {
        abs(xval)/16
    }
    ti <- c(xlim[1L] - dr, xval, xlim[2L] + dr)
    ti.l <- ti[-length(ti)]
    ti.r <- ti[-1L]
    y <- x(0.5 * (ti.l + ti.r))
    x=c(rbind(ti.l,ti.r))
    y1=rep(y,each=2)
    df=data.frame(x=x,y=y1)
    df
}


#'Draw Kaplan-Meier Survival Curve For Follow-up Data
#'@param model An object of class gam
#'@param data A data.frame
#'@importFrom stats coef vcov stepfun
#'@importFrom ggplot2 ggplot
#'@export
#'@examples
#' library(survival)
#' library(mgcv)
#'pbcseq$status1 <- as.numeric(pbcseq$status==2) ## deaths
#'pb <- tdpois(pbcseq) ## conversion
#'pb$tf <- factor(pb$futime) ## add factor for event time
#'b <- bam(z ~ tf - 1 + trt + s(sqrt(protime)) + s(platelet) + s(age) + s(bili) +
#'s(albumin) + s(sqrt(ast)),family=poisson,data=pb,discrete=TRUE,nthreads=2)
#'data=pbcseq[pbcseq$id %in% c(10,25,66),]
#'data=pbcseq[pbcseq$id %in% c(25),]
#'drawFUSurv(b,data)
drawFUSurv=function(model,data){
    te <- as.numeric(levels(model$model$tf))
    select=unique(data$id)

    for(i in 1:length(select)){
        di <- data[data$id==select[i],]
        pd <- data.frame(lapply(X=di,FUN=app,t=di$day,to=te))
        pd$tf <- factor(te)

        X <- predict(model,newdata=pd,type="lpmatrix")
        eta <- drop(X%*%coef(model))
        H <- cumsum(exp(eta))
        J <- apply(exp(eta)*X,2,cumsum)
        se <- diag(J%*%vcov(model)%*%t(J))^.5

        df=stepfun2df(stepfun(te,c(1,exp(-H))))
        df1=stepfun2df(stepfun(te,c(1,exp(-H+se))))
        df2=stepfun2df(stepfun(te,c(1,exp(-H-se))))
        df$ymax=df1$y
        df$ymin=df2$y
        df$id=select[i]
        if(i==1){
            final=df
        } else{
            final=rbind(final,df)
        }

    }

    final$id=factor(final$id)
    final$x[final$x<0]=0
    p<- ggplot(final,aes_string(x="x"))+
        geom_line(aes_string(y="y",group="id",color="id"))+
        geom_ribbon(aes_string(ymax="ymax",ymin="ymin",fill="id"),alpha=0.4)+
        theme_bw()+
        labs(y="Cumulative Survival",x="days")+
        theme(legend.position="top")

    if(length(select)==1){
        fu=data.frame(x=data$day[data$id==select[1]]) ## measurement times
        p<-p+geom_rug(data=fu)
    }
    p
}


#'Draw follow-up plot for laboratory data
#'@param model An object of class gam
#'@param data A data.frame
#'@param which A character vector Names of columns to drae
#'@param point logical Whether or not draw point
#'@param line logical Whether or not draw line
#'@importFrom tidyr pivot_longer
#'@importFrom  ggplot2 geom_rug
#'@export
#'@examples
#' library(survival)
#' library(mgcv)
#'pbcseq$status1 <- as.numeric(pbcseq$status==2) ## deaths
#'pb <- tdpois(pbcseq) ## conversion
#'pb$tf <- factor(pb$futime) ## add factor for event time
#'b <- bam(z ~ tf - 1 + trt + s(sqrt(protime)) + s(platelet) + s(age) + s(bili) +
#'s(albumin) + s(sqrt(ast)),family=poisson,data=pb,discrete=TRUE,nthreads=2)
#'data=pbcseq[pbcseq$id==25,]
#'drawFUData(b,data,which=c("protime","platelet","bili","albumin","ast"))
drawFUData=function(model,data,which=c("protime","platelet","bili","albumin","ast"),
                    point=TRUE,line=TRUE){

    te <- as.numeric(levels(model$model$tf))
    pd <- data.frame(lapply(X=data,FUN=app,t=data$day,to=te))
    df=data.frame(x=data$day)
    df1=data.frame(x=te)
    for(i in 1:length(which)){
        df[[which[i]]]=data[[which[i]]]
        df1[[which[i]]]=pd[[which[i]]]
    }
    longdf=pivot_longer(df,cols=all_of(which))
    longdf$name=factor(longdf$name,levels=which)
    longdf1=pivot_longer(df1,cols=all_of(which))
    longdf1$name=factor(longdf1$name,levels=which)
    p<-ggplot(data=longdf,aes_string(x="x",y="value"))
    if(point) p <-p+geom_point()
    if(line)  p <-p+geom_line(data=longdf1)
    p<-p+facet_wrap("name",scales="free_y")+
        labs(x="days",y="")
    p
}
