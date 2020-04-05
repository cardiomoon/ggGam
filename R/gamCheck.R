#'Perform diagnostics for a fitted gam model
#'@param b	A fitted gam object as produced by gam().
#'@param k.sample Newmeric Above this k testing uses a random sub-sample of data.
#'@param k.rep Numeric How many re-shuffles to do to get p-value for k testing.
#'@param print logical Whether or not print result
#'@importFrom stats printCoefmat
#'@importFrom mgcv k.check
#'@export
#'@examples
#'library(mgcv)
#'set.seed(0)
#'dat <- gamSim(1,n=200)
#'b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat)
#'gamCheck(b)
gamCheck=function (b, k.sample = 5000, k.rep = 200,print=TRUE)
{
    gamm <- !(b$method %in% c("GCV", "GACV", "UBRE", "REML",
                              "ML", "P-ML", "P-REML", "fREML"))
    cat2=function(...){
        if(print) cat(...)
    }
    print2=function(...){
        if(print) print(...)
    }

    if (gamm) {
        cat2("\n'gamm' based fit - care required with interpretation.")
        cat2("\nChecks based on working residuals may be misleading.")
    } else {
        cat2("\nMethod:", b$method, "  Optimizer:", b$optimizer)
        if (!is.null(b$outer.info)) {
            if (b$optimizer[2] %in% c("newton", "bfgs")) {
                boi <- b$outer.info
                cat2("\n", boi$conv, " after ", boi$iter, " iteration",
                    sep = "")
                if (boi$iter == 1)
                    cat2(".")
                else cat2("s.")
                cat2("\nGradient range [", min(boi$grad), ",",
                    max(boi$grad), "]", sep = "")
                cat2("\n(score ", b$gcv.ubre, " & scale ", b$sig2,
                    ").", sep = "")
                ev <- eigen(boi$hess)$values
                if (min(ev) > 0)
                    cat2("\nHessian positive definite, ")
                else cat2("\n")
                cat2("eigenvalue range [", min(ev), ",", max(ev),
                    "].\n", sep = "")
            }
            else {
                cat2("\n")
                print2(b$outer.info)
            }
        }
        else {
            if (length(b$sp) == 0)
                cat2("\nModel required no smoothing parameter selection")
            else {
                cat2("\nSmoothing parameter selection converged after",
                    b$mgcv.conv$iter, "iteration")
                if (b$mgcv.conv$iter > 1)
                    cat2("s")
                if (!b$mgcv.conv$fully.converged)
                    cat2(" by steepest\ndescent step failure.\n")
                else cat2(".\n")
                cat2("The RMS", b$method, "score gradient at convergence was",
                    b$mgcv.conv$rms.grad, ".\n")
                if (b$mgcv.conv$hess.pos.def)
                    cat2("The Hessian was positive definite.\n")
                else cat2("The Hessian was not positive definite.\n")
            }
        }
        if (!is.null(b$rank)) {
            cat2("Model rank = ", b$rank, "/", length(b$coefficients),
                "\n")
        }
    }
    cat2("\n")
    kchck <- k.check(b, subsample = k.sample, n.rep = k.rep)
    if (!is.null(kchck)) {
        cat2("Basis dimension (k) checking results. Low p-value (k-index<1) may\n")
        cat2("indicate that k is too low, especially if edf is close to k'.\n\n")
        if(print) printCoefmat(kchck, digits = 3)
    }
    invisible(kchck)

}
