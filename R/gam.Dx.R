#' Draw diagnostic plots for a fitted gam model
#' @param b a fitted gam object as produced by gam().
#' @param type type of residuals
#' @param rep arguments passed to qq.gam()
#' @param level arguments passed to qq.gam()
#' @param select numeric choices of plots
#' @importFrom ggplot2 ggplot geom_point labs geom_hline geom_histogram theme_bw aes_string
#' @importFrom ggpubr ggarrange
#' @importFrom grDevices nclass.Sturges
#' @importFrom stats approx fitted napredict quantile residuals
#' @export
#' @examples
#' require(mgcv)
#' data(pbc,package="survival")
#' pbc$status1 <- as.numeric(pbc$status==2)
#' pbc$stage <- factor(pbc$stage)
#' b = gam(time ~ trt +sex + s(sqrt(protime))+s(platelet)+ s(age)+s(bili)+s(albumin),
#'    weights=status1, family=cox.ph, data=pbc, method="REML")
#' gam.Dx(b)
#' gam.Dx(b,select=4)
gam.Dx=function (b, type = c("deviance", "pearson","response"),rep=0,level=0.9,select=NULL)
{
    if(missing(select)) select=1:4

    type=match.arg(type)
    resid <- residuals(b, type = type)
    linpred <- if (is.matrix(b$linear.predictors) && !is.matrix(resid)) {
        napredict(b$na.action, b$linear.predictors[, 1])
    } else {
        napredict(b$na.action, b$linear.predictors)
    }
    p1=qqgam(b, rep = rep, level = level, type = type)+theme_bw()

    df=data.frame(x=linpred,y=resid)
    p2=ggplot(df,aes_string(x="x",y="y"))+geom_point(alpha=0.6)+labs(title="Residuals vs. linear pred.", x = "linear predictor",  y = "residuals")+geom_hline(yintercept=0,color="red")+theme_bw()

    df <- data.frame(residuals = residuals(b, type = type))
    n_bins=nclass.Sturges(df[["residuals"]])+2
    p3=ggplot(df,aes_string(x="residuals"))+
        geom_histogram(bins=n_bins,color="black",fill="grey80")+
        labs(title="Histogram of residuals",x="Deviance residuals",y="Frequency")+
        theme_bw()
    fv <- if (inherits(b$family, "extended.family")) {
        predict(b, type = "response")
    } else {
        fitted(b)
    }
    if (is.matrix(fv) && !is.matrix(b$y))
        fv <- fv[, 1]
    df=data.frame(x=fv,y=napredict(b$na.action, b$y))
    p4=ggplot(df,aes_string(x="x",y="y"))+
        geom_point(alpha=0.6)+
        labs(x = "Fitted Values", y = "Response", title = "Response vs. Fitted Values")+
        theme_bw()

    plotlist=list()
    selected=1
    if(1 %in% select){
          plotlist[[selected]]=p1
          selected=selected+1
    }
    if(2 %in% select){
        plotlist[[selected]]=p2
        selected=selected+1
    }
    if(3 %in% select){
        plotlist[[selected]]=p3
        selected=selected+1
    }
    if(4 %in% select){
        plotlist[[selected]]=p4
        selected=selected+1
    }

    if(length(select)==1){ncol=1;nrow=1}
    else if(length(select)==2){ncol=2;nrow=1}
    else {ncol=2;nrow=2}
    if(length(select)==1){
        plotlist[[1]]
    } else{
       ggpubr::ggarrange(plotlist=plotlist, ncol = ncol,nrow=nrow)
    }

}


#'QQ plots for gam model residuals
#'@param object	a fitted gam object as produced by gam()
#'@param rep How many replicate datasets to generate to simulate quantiles of the residual distribution. 0 results in an efficient simulation free method for direct calculation, if this is possible for the object family.
#'@param level	If simulation is used for the quantiles, then reference intervals can be provided for the QQ-plot, this specifies the level. 0 or less for no intervals, 1 or more to simply plot the QQ plot for each replicate generated.
#'@param s.rep how many times to randomize uniform quantiles to data under direct computation.
#'@param type what sort of residuals should be plotted?
#'@importFrom ggplot2 stat_qq stat_qq_line geom_abline ylim
#'@importFrom mgcv fix.family.qf fix.family.rd
#'@return A ggplot
qqgam=function (object, rep = 0, level = 0.9, s.rep = 10,
                type = c("deviance","pearson", "response"))
{
    # rep = 0; level = 0.9; s.rep = 10; type = "deviance"
    type=match.arg(type)
    ylab <- paste(type, "residuals")
    if (inherits(object, c("glm", "gam"))) {
        if (is.null(object$sig2))
            object$sig2 <- summary(object)$dispersion
    } else stop("object is not a glm or gam")
    object$na.action <- NULL
    D <- residuals(object, type = type)
    if (object$method %in% c("PQL", "lme.ML", "lme.REML", "lmer.REML",
                             "lmer.ML", "glmer.ML")) {
        df=data.frame(y=D)
        ggplot(df, aes_string(sample = "y"))+
            stat_qq(alpha=0.6) + stat_qq_line()

        return()
    }
    lim <- Dq <- NULL
    if (rep == 0) {
        fam <- fix.family.qf(object$family)
        if (is.null(fam$qf))
            rep <- 50
        level <- 0
    }
    n <- length(D)
    if (rep > 0) {
        fam <- fix.family.rd(object$family)
        if (!is.null(fam$rd)) {
            dm <- matrix(0, n, rep)
            for (i in 1:rep) {
                yr <- fam$rd(object$fitted.values, object$prior.weights,
                             object$sig2)
                object$y <- yr
                dm[, i] <- sort(residuals(object, type = type))
            }
            Dq <- quantile(as.numeric(dm), (1:n - 0.5)/n)
            alpha <- (1 - level)/2
            if (alpha > 0.5 || alpha < 0)
                alpha <- 0.05
            if (level > 0 && level < 1)
                lim <- apply(dm, 1, FUN = quantile, p = c(alpha,
                                                          1 - alpha))
            else if (level >= 1)
                lim <- level
        }
    } else {
        U <- (1:n - 0.5)/n
        if (!is.null(fam$qf)) {
            dm <- matrix(0, n, s.rep)
            for (i in 1:s.rep) {
                U <- sample(U, n)
                q0 <- fam$qf(U, object$fitted.values, object$prior.weights,
                             object$sig2)
                object$y <- q0
                dm[, i] <- sort(residuals(object, type = type))
            }
            Dq <- sort(rowMeans(dm))
        }
    }
    Dq
    if (!is.null(Dq)) {
        sx=sort(Dq)
        sy=sort(D)
        lenx <- length(sx)
        leny <- length(sy)
        if (leny < lenx)
            sx <- approx(1L:lenx, sx, n = leny)$y
        if (leny > lenx)
            sy <- approx(1L:leny, sy, n = lenx)$y
        df=data.frame(x=sx,y=sy)
        ggplot(df,aes_string(x="x",y="y"))+geom_point(alpha=0.6)+
            labs(x="Theoretical quantiles",y="Deviance residuals",
                 title="QQ plot of residuals")+
            ylim(range(c(lim, D)))+
            geom_abline(slope=1,intercept=0,color="red")

    } else {
        df=data.frame(y=D)
        ggplot(df, aes_string(sample = "y"))+
            stat_qq(alpha=0.6) + stat_qq_line(color="red")

    }
}

