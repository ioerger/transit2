import scipy
import numpy
import heapq
import statsmodels.stats.multitest

import time
import sys
import collections
import functools

hasR = False
try:
    import rpy2.robjects
    hasR = True
except Exception as e:
    hasR = False

if hasR:
    from rpy2.robjects import r, DataFrame, globalenv, IntVector, FloatVector, StrVector, packages as rpackages

from pytransit.analysis import base
import pytransit
import pytransit.transit_tools as transit_tools
import pytransit.tnseq_tools as tnseq_tools
import pytransit.norm_tools as norm_tools

############# GUI ELEMENTS ##################

short_name = "ZINB"
long_name = "ZINB"
short_desc = "Perform ZINB analysis"
long_desc = """Perform ZINB analysis"""
EOL = "\n"
DEBUG = False
GENE = None
SEPARATOR = '\1' # for making names that combine conditions and interactions; try not to use a char a user might have in a condition name

transposons = ["", ""]
columns = []

class ZinbAnalysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(self, short_name, long_name, short_desc, long_desc, transposons, ZinbMethod)
def main():
    print("ZINB example")

class ZinbMethod(base.MultiConditionMethod):
    """
    Zinb
    """
    def __init__(self, combined_wig, metadata, annotation, normalization, output_file, ignored_conditions=[], included_conditions=[], winz=False, nterm=5.0, cterm=5.0, condition="Condition", covars=[], interactions = []):
        base.MultiConditionMethod.__init__(self, short_name, long_name, short_desc, long_desc, combined_wig, metadata, annotation, output_file,
                normalization=normalization, ignored_conditions=ignored_conditions, included_conditions=included_conditions, nterm=nterm, cterm=cterm)
        self.winz = winz
        self.covars = covars
        self.interactions = interactions
        self.condition = condition

    @classmethod
    def fromargs(self, rawargs):
        if not hasR:
            print("Error: R and rpy2 (~= 3.0) required to run ZINB analysis.")
            print("After installing R, you can install rpy2 using the command \"pip install 'rpy2~=3.0'\"")
            sys.exit(0)

        (args, kwargs) = transit_tools.cleanargs(rawargs)

        if (kwargs.get('-help', False) or kwargs.get('h', False)):
            print(ZinbMethod.usage_string())
            sys.exit(0)

        if (kwargs.get('v', False)):
            global DEBUG
            DEBUG = True

        if (kwargs.get('-gene', False)):
            global GENE
            GENE = kwargs.get('-gene', None)

        combined_wig = args[0]
        metadata = args[1]
        annotation = args[2]
        output_file = args[3]
        normalization = kwargs.get("n", "TTR")
        NTerminus = float(kwargs.get("iN", 5.0))
        CTerminus = float(kwargs.get("iC", 5.0))
        condition = kwargs.get("-condition", "Condition")
        covars = list(filter(None, kwargs.get("-covars", "").split(",")))
        interactions = list(filter(None, kwargs.get("-interactions", "").split(",")))
        winz = True if "w" in kwargs else False
        ignored_conditions = list(filter(None, kwargs.get("-ignore-conditions", "").split(",")))
        included_conditions = list(filter(None, kwargs.get("-include-conditions", "").split(",")))
        if len(included_conditions) > 0 and len(ignored_conditions) > 0:
            self.transit_error("Cannot use both include-conditions and ignore-conditions flags")
            print(ZinbMethod.usage_string())
            sys.exit(0)

        return self(combined_wig, metadata, annotation, normalization, output_file, ignored_conditions, included_conditions, winz, NTerminus, CTerminus, condition, covars, interactions)

    def wigs_to_conditions(self, conditionsByFile, filenamesInCombWig):
        """
            Returns list of conditions corresponding to given wigfiles.
            ({FileName: Condition}, [FileName]) -> [Condition]
            Condition :: [String]
        """
        return [conditionsByFile.get(f, self.unknown_cond_flag) for f in filenamesInCombWig]

    def wigs_to_covariates(self, covariatesMap, filenamesInCombWig):
        """
            Returns list of covariate lists. Each covariate list consists of covariates corresponding to given wigfiles.
            ([{FileName: Covar}], [FileName]) -> [[Covar]]
            Covar :: String
        """
        try:
            return [[covarsByFile[f]
                        for f in filenamesInCombWig]
                        for covarsByFile in covariatesMap]
        except KeyError:
            self.transit_error("Error: Covariates not found for file {0}".format(f))
            sys.exit(0)

    def wigs_to_interactions(self, interactionsMap, filenamesInCombWig):
        """
            Returns list of interaction lists. Each interaction list consists of covariates corresponding to given wigfiles.
            ([{FileName: Interaction}], [FileName]) -> [[Interaction]]
            Interaction :: String
        """
        try:
            return [[covarsByFile[f]
                        for f in filenamesInCombWig]
                        for covarsByFile in interactionsMap]
        except KeyError:
            self.transit_error("Error: Interaction var not found for file {0}".format(f))
            sys.exit(0)

    def stats_for_gene(self, siteIndexes, groupWigIndexMap, data):
        """
            Returns a dictionary of {Group: {Mean, NzMean, NzPerc}}
            ([SiteIndex], [Condition], [WigData]) -> [{Condition: Number}]
            SiteIndex :: Number
            WigData :: [Number]
            Group :: String (combination of '<interaction>_<condition>')
        """
        nonzero = lambda xs: xs[numpy.nonzero(xs)]
        nzperc = lambda xs: numpy.count_nonzero(xs)/float(xs.size)

        means = {}
        nz_means = {}
        nz_percs = {}

        for (group, wigIndexes) in groupWigIndexMap.items():
            if (len(siteIndexes) == 0): # If no TA sites, write 0
                means[group] = 0
                nz_means[group] = 0
                nz_percs[group] = 0
            else:
                arr = data[wigIndexes][:, siteIndexes]
                means[group] = numpy.mean(arr) if len(arr) > 0 else 0
                nonzero_arr = nonzero(arr)
                nz_means[group] = numpy.mean(nonzero_arr) if len(nonzero_arr) > 0 else 0
                nz_percs[group] = nzperc(arr)

        return {'mean': means, 'nz_mean': nz_means, 'nz_perc': nz_percs}

    def stats_by_rv(self, data, RvSiteindexesMap, genes, conditions, interactions):
        """
            Returns Dictionary of Stats by condition for each Rv
            ([[Wigdata]], {Rv: SiteIndex}, [Gene], [Condition], [Interaction]) -> {Rv: {Condition: Number}}
            Wigdata :: [Number]
            SiteIndex :: Number
            Gene :: {start, end, rv, gene, strand}
            Condition :: String
            Interaction :: String
        """

        ## Group wigfiles by (interaction, condition) pair
        ## {'<interaction>_<condition>': [Wigindexes]}
        groupWigIndexMap = collections.defaultdict(lambda: [])
        for i, conditionForWig in enumerate(conditions):
            if (len(interactions) > 0):
                for interaction in interactions:
                    groupName = conditionForWig + SEPARATOR + interaction[i] 
                    groupWigIndexMap[groupName].append(i)
            else:
                groupName = conditionForWig
                groupWigIndexMap[groupName].append(i)

        statsByRv = {}
        for gene in genes:
            Rv = gene["rv"]
            statsByRv[Rv] = self.stats_for_gene(RvSiteindexesMap[Rv], groupWigIndexMap, data)

        ## TODO :: Any ordering to follow?
        statGroupNames = groupWigIndexMap.keys()
        return statsByRv, statGroupNames

    def global_stats_for_rep(self, data):
        """
        Returns the logit zero percentage and nz_mean for each replicate.
            [[WigData]] -> [[Number] ,[Number]]
        """

        logit_zero_perc = []
        nz_mean = []
        for wig in data:
            zero_perc = (wig.size - numpy.count_nonzero(wig))/float(wig.size)
            logit_zero_perc.append(numpy.log(zero_perc/(1 - zero_perc)))
            nz_mean.append(numpy.mean(wig[numpy.nonzero(wig)]))
        return [numpy.array(logit_zero_perc), numpy.array(nz_mean)]

    def group_by_condition(self, wigList, conditions):
        """
            Returns array of datasets, where each dataset corresponds to one condition.
            ([[Wigdata]], [Condition]) -> [[DataForCondition]]
            Wigdata :: [Number]
            Condition :: String
            DataForCondition :: [Number]
        """
        countsByCondition = collections.defaultdict(lambda: [])
        for i, c in enumerate(conditions):
          countsByCondition[c].append(wigList[i])

        return [numpy.array(v).flatten() for v in countsByCondition.values()]

    def melt_data(self, readCountsForRv, conditions, covariates, interactions, NZMeanByRep, LogZPercByRep):
        rvSitesLength = len(readCountsForRv[0])
        repeatAndFlatten = lambda xs: numpy.repeat(xs, rvSitesLength)

        return [
                numpy.concatenate(readCountsForRv).astype(int),
                repeatAndFlatten(conditions),
                list(map(repeatAndFlatten, covariates)),
                list(map(repeatAndFlatten, interactions)),
                repeatAndFlatten(NZMeanByRep),
                repeatAndFlatten(LogZPercByRep)
               ]

    def def_r_zinb_signif(self):
        r('''
	    model_offset_2 <- function(x, terms = NULL, offset = TRUE)
	    ## allow optionally different terms
	    ## potentially exclude "(offset)"
	    {
	      if(is.null(terms)) terms <- attr(x, "terms")
	      offsets <- attr(terms, "offset")
	      if(length(offsets) > 0) {
		ans <- if(offset) x$"(offset)" else NULL
		if(is.null(ans)) ans <- 0
		for(i in offsets) ans <- ans + x[[deparse(attr(terms, "variables")[[i + 1]])]]
		ans
	      }
	      else {
		ans <- if(offset) x$"(offset)" else NULL
	      }
	      if(!is.null(ans) && !is.numeric(ans)) stop("'offset' must be numeric")
	      ans
	    }

            ######################################
            # this code was copied from the pscl library (v1.5.2) by Simon Jackman
            #   and updated how control args such as upper/lower are passed to the optimizer
            # https://github.com/cran/pscl/tree/master/R

	    zeroinfl_sid = function (formula, data, subset, na.action, weights, offset,
                                     dist = c("poisson", "negbin", "geometric"), 
                                     link = c("logit", "probit", "cloglog", "cauchit", "log"), control = zeroinfl.control(...),
                                     model = TRUE, y = TRUE, x = FALSE, ...)
	    {
	      ziPoisson <- function(parms) {
		mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
		phi <- as.vector(linkinv(Z %*% parms[(kx + 1):(kx + kz)] +
					   offsetz))
		loglik0 <- log(phi + exp(log(1 - phi) - mu))
		loglik1 <- log(1 - phi) + dpois(Y, lambda = mu, log = TRUE)
		loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] *
								 loglik1[Y1])
		loglik
	      }
	      ziNegBin <- function(parms) {
		mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
		phi <- as.vector(linkinv(Z %*% parms[(kx + 1):(kx + kz)] +
					   offsetz))
		theta <- exp(parms[(kx + kz) + 1])
		loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0,
										 size = theta, mu = mu, log = TRUE))))
		loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y,
								   size = theta, mu = mu, log = TRUE))
		loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] *
								 loglik1[Y1])
		loglik
	      }
	      ziGeom <- function(parms) ziNegBin(c(parms, 0))
	      gradPoisson <- function(parms) {
		eta <- as.vector(X %*% parms[1:kx] + offsetx)
		mu <- exp(eta)
		etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
		muz <- linkinv(etaz)
		clogdens0 <- -mu
		dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) +
							    clogdens0)
		wres_count <- ifelse(Y1, Y - mu, -exp(-log(dens0) + log(1 -
									  muz) + clogdens0 + log(mu)))
		wres_zero <- ifelse(Y1, -1/(1 - muz) * linkobj$mu.eta(etaz),
				    (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
		colSums(cbind(wres_count * weights * X, wres_zero * weights *
				Z))
	      }
	      gradGeom <- function(parms) {
		eta <- as.vector(X %*% parms[1:kx] + offsetx)
		mu <- exp(eta)
		etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
		muz <- linkinv(etaz)
		clogdens0 <- dnbinom(0, size = 1, mu = mu, log = TRUE)
		dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) +
							    clogdens0)
		wres_count <- ifelse(Y1, Y - mu * (Y + 1)/(mu + 1), -exp(-log(dens0) +
									   log(1 - muz) + clogdens0 - log(mu + 1) + log(mu)))
		wres_zero <- ifelse(Y1, -1/(1 - muz) * linkobj$mu.eta(etaz),
				    (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
		colSums(cbind(wres_count * weights * X, wres_zero * weights *
				Z))
	      }
	      gradNegBin <- function(parms) {
		eta <- as.vector(X %*% parms[1:kx] + offsetx)
		mu <- exp(eta)
		etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
		muz <- linkinv(etaz)
		theta <- exp(parms[(kx + kz) + 1])
		clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
		dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) +
							    clogdens0)
		wres_count <- ifelse(Y1, Y - mu * (Y + theta)/(mu + theta),
				     -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) -
					    log(mu + theta) + log(mu)))
		wres_zero <- ifelse(Y1, -1/(1 - muz) * linkobj$mu.eta(etaz),
				    (linkobj$mu.eta(etaz) - exp(clogdens0) * linkobj$mu.eta(etaz))/dens0)
		wres_theta <- theta * ifelse(Y1, digamma(Y + theta) -
					       digamma(theta) + log(theta) - log(mu + theta) + 1 -
					       (Y + theta)/(mu + theta), exp(-log(dens0) + log(1 -
												 muz) + clogdens0) * (log(theta) - log(mu + theta) +
															1 - theta/(mu + theta)))
		colSums(cbind(wres_count * weights * X, wres_zero * weights *
				Z, wres_theta))
	      }
	      dist <- match.arg(dist)
	      loglikfun <- switch(dist, poisson = ziPoisson, geometric = ziGeom,
				  negbin = ziNegBin)
	      gradfun <- switch(dist, poisson = gradPoisson, geometric = gradGeom,
				negbin = gradNegBin)
	      linkstr <- match.arg(link)
	      linkobj <- make.link(linkstr)
	      linkinv <- linkobj$linkinv
	      if (control$trace)
		cat("Zero-inflated Count Model\n", paste("count model:",
							 dist, "with log link\n"), paste("zero-inflation model: binomial with",
											 linkstr, "link\n"), sep = "")
	      cl <- match.call()
	      if (missing(data))
		data <- environment(formula)
	      mf <- match.call(expand.dots = FALSE)
	      m <- match(c("formula", "data", "subset", "na.action", "weights",
			   "offset"), names(mf), 0)
	      mf <- mf[c(1, m)]
	      mf$drop.unused.levels <- TRUE
	      if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]],
							as.name("|"))) {
		ff <- formula
		formula[[3]][1] <- call("+")
		mf$formula <- formula
		ffc <- . ~ .
		ffz <- ~.
		ffc[[2]] <- ff[[2]]
		ffc[[3]] <- ff[[3]][[2]]
		ffz[[3]] <- ff[[3]][[3]]
		ffz[[2]] <- NULL
	      }
	      else {
		ffz <- ffc <- ff <- formula
		ffz[[2]] <- NULL
	      }
	      if (inherits(try(terms(ffz), silent = TRUE), "try-error")) {
		ffz <- eval(parse(text = sprintf(paste("%s -", deparse(ffc[[2]])),
						 deparse(ffz))))
	      }
	      mf[[1]] <- as.name("model.frame")
	      mf <- eval(mf, parent.frame())
	      mt <- attr(mf, "terms")
	      mtX <- terms(ffc, data = data)
	      X <- model.matrix(mtX, mf)
	      mtZ <- terms(ffz, data = data)
	      mtZ <- terms(update(mtZ, ~.), data = data)
	      Z <- model.matrix(mtZ, mf)
	      Y <- model.response(mf, "numeric")
	      if (length(Y) < 1)
		stop("empty model")
	      if (all(Y > 0))
		stop("invalid dependent variable, minimum count is not zero")
	      if (!isTRUE(all.equal(as.vector(Y), as.integer(round(Y +
								   0.001)))))
		stop("invalid dependent variable, non-integer values")
	      Y <- as.integer(round(Y + 0.001))
	      if (any(Y < 0))
		stop("invalid dependent variable, negative counts")
	      if (control$trace) {
		cat("dependent variable:\n")
		tab <- table(factor(Y, levels = 0:max(Y)), exclude = NULL)
		names(dimnames(tab)) <- NULL
		print(tab)
	      }
	      n <- length(Y)
	      kx <- NCOL(X)
	      kz <- NCOL(Z)
	      Y0 <- Y <= 0
	      Y1 <- Y > 0
	      weights <- model.weights(mf)
	      if (is.null(weights))
		weights <- 1
	      if (length(weights) == 1)
		weights <- rep.int(weights, n)
	      weights <- as.vector(weights)
	      names(weights) <- rownames(mf)
	      offsetx <- model_offset_2(mf, terms = mtX, offset = TRUE)
	      if (is.null(offsetx))
		offsetx <- 0
	      if (length(offsetx) == 1)
		offsetx <- rep.int(offsetx, n)
	      offsetx <- as.vector(offsetx)
	      offsetz <- model_offset_2(mf, terms = mtZ, offset = FALSE)
	      if (is.null(offsetz))
		offsetz <- 0
	      if (length(offsetz) == 1)
		offsetz <- rep.int(offsetz, n)
	      offsetz <- as.vector(offsetz)
	      start <- control$start
	      if (!is.null(start)) {
		valid <- TRUE
		if (!("count" %in% names(start))) {
		  valid <- FALSE
		  warning("invalid starting values, count model coefficients not specified")
		  start$count <- rep.int(0, kx)
		}
		if (!("zero" %in% names(start))) {
		  valid <- FALSE
		  warning("invalid starting values, zero-inflation model coefficients not specified")
		  start$zero <- rep.int(0, kz)
		}
		if (length(start$count) != kx) {
		  valid <- FALSE
		  warning("invalid starting values, wrong number of count model coefficients")
		}
		if (length(start$zero) != kz) {
		  valid <- FALSE
		  warning("invalid starting values, wrong number of zero-inflation model coefficients")
		}
		if (dist == "negbin") {
		  if (!("theta" %in% names(start)))
		    start$theta <- 1
		  start <- list(count = start$count, zero = start$zero,
				theta = as.vector(start$theta[1]))
		}
		else {
		  start <- list(count = start$count, zero = start$zero)
		}
		if (!valid)
		  start <- NULL
	      }
	      if (is.null(start)) {
		if (control$trace)
		  cat("generating starting values...")
		model_count <- glm.fit(X, Y, family = poisson(), weights = weights,
				       offset = offsetx)
		model_zero <- glm.fit(Z, as.integer(Y0), weights = weights,
				      family = binomial(link = linkstr), offset = offsetz)
		start <- list(count = model_count$coefficients, zero = model_zero$coefficients)
		if (dist == "negbin")
		  start$theta <- 1
		if (control$EM & dist == "poisson") {
		  mui <- model_count$fitted
		  probi <- model_zero$fitted
		  probi <- probi/(probi + (1 - probi) * dpois(0, mui))
		  probi[Y1] <- 0
		  ll_new <- loglikfun(c(start$count, start$zero))
		  ll_old <- 2 * ll_new
		  while (abs((ll_old - ll_new)/ll_old) > control$reltol) {
		    ll_old <- ll_new
		    model_count <- glm.fit(X, Y, weights = weights *
					     (1 - probi), offset = offsetx, family = poisson(),
					   start = start$count)
		    model_zero <- suppressWarnings(glm.fit(Z, probi,
							   weights = weights, offset = offsetz, family = binomial(link = linkstr),
							   start = start$zero))
		    mui <- model_count$fitted
		    probi <- model_zero$fitted
		    probi <- probi/(probi + (1 - probi) * dpois(0,
								mui))
		    probi[Y1] <- 0
		    start <- list(count = model_count$coefficients,
				  zero = model_zero$coefficients)
		    ll_new <- loglikfun(c(start$count, start$zero))
		  }
		}
		if (control$EM & dist == "geometric") {
		  mui <- model_count$fitted
		  probi <- model_zero$fitted
		  probi <- probi/(probi + (1 - probi) * dnbinom(0,
								size = 1, mu = mui))
		  probi[Y1] <- 0
		  ll_new <- loglikfun(c(start$count, start$zero))
		  ll_old <- 2 * ll_new
		  while (abs((ll_old - ll_new)/ll_old) > control$reltol) {
		    ll_old <- ll_new
		    model_count <- suppressWarnings(glm.fit(X, Y,
							    weights = weights * (1 - probi), offset = offsetx,
							    family = MASS::negative.binomial(1), start = start$count))
		    model_zero <- suppressWarnings(glm.fit(Z, probi,
							   weights = weights, offset = offsetz, family = binomial(link = linkstr),
							   start = start$zero))
		    start <- list(count = model_count$coefficients,
				  zero = model_zero$coefficients)
		    mui <- model_count$fitted
		    probi <- model_zero$fitted
		    probi <- probi/(probi + (1 - probi) * dnbinom(0,
								  size = 1, mu = mui))
		    probi[Y1] <- 0
		    ll_new <- loglikfun(c(start$count, start$zero))
		  }
		}
		if (control$EM & dist == "negbin") {
		  mui <- model_count$fitted
		  probi <- model_zero$fitted
		  probi <- probi/(probi + (1 - probi) * dnbinom(0,
								size = start$theta, mu = mui))
		  probi[Y1] <- 0
		  ll_new <- loglikfun(c(start$count, start$zero, log(start$theta)))
		  ll_old <- 2 * ll_new
		  offset <- offsetx
		  while (abs((ll_old - ll_new)/ll_old) > control$reltol) {
		    ll_old <- ll_new
		    model_count <- suppressWarnings(glm.nb(Y ~ 0 +
							     X + offset(offset), weights = weights * (1 -
													probi), start = start$count, init.theta = start$theta))
		    model_zero <- suppressWarnings(glm.fit(Z, probi,
							   weights = weights, offset = offsetz, family = binomial(link = linkstr),
							   start = start$zero))
		    start <- list(count = model_count$coefficients,
				  zero = model_zero$coefficients, theta = model_count$theta)
		    mui <- model_count$fitted
		    probi <- model_zero$fitted
		    probi <- probi/(probi + (1 - probi) * dnbinom(0,
								  size = start$theta, mu = mui))
		    probi[Y1] <- 0
		    ll_new <- loglikfun(c(start$count, start$zero,
					  log(start$theta)))
		  }
		}
		if (control$trace)
		  cat("done\n")
	      }
	      if (control$trace)
		cat("calling optim() for ML estimation:\n")
	      method <- control$method
	      hessian <- control$hessian
	      ocontrol <- control
	      control$method <- control$hessian <- control$EM <- control$start <- NULL
	      fit <- optim(fn = loglikfun, gr = gradfun, par = c(start$count,
								 start$zero, if (dist == "negbin") log(start$theta) else NULL),
                           # this is the only line that changed...
                           # for generality, should also check whether control$upper/lower are defined
			   #method = method, hessian = hessian, control = control)
			   method = method, hessian = hessian, control = control, upper = control$upper, lower=control$lower)
	      if (fit$convergence > 0)
		warning("optimization failed to converge")
	      coefc <- fit$par[1:kx]
	      names(coefc) <- names(start$count) <- colnames(X)
	      coefz <- fit$par[(kx + 1):(kx + kz)]
	      names(coefz) <- names(start$zero) <- colnames(Z)
	      vc <- -solve(as.matrix(fit$hessian))
	      if (dist == "negbin") {
		np <- kx + kz + 1
		theta <- as.vector(exp(fit$par[np]))
		SE.logtheta <- as.vector(sqrt(diag(vc)[np]))
		vc <- vc[-np, -np, drop = FALSE]
	      }
	      else {
		theta <- NULL
		SE.logtheta <- NULL
	      }
	      colnames(vc) <- rownames(vc) <- c(paste("count", colnames(X),
						      sep = "_"), paste("zero", colnames(Z), sep = "_"))
	      mu <- exp(X %*% coefc + offsetx)[, 1]
	      phi <- linkinv(Z %*% coefz + offsetz)[, 1]
	      Yhat <- (1 - phi) * mu
	      res <- sqrt(weights) * (Y - Yhat)
	      nobs <- sum(weights > 0)
	      rval <- list(coefficients = list(count = coefc, zero = coefz),
			   residuals = res, fitted.values = Yhat, optim = fit, method = method,
			   control = ocontrol, start = start, weights = if (identical(as.vector(weights),
										      rep.int(1L, n))) NULL else weights, offset = list(count = if (identical(offsetx,
																			      rep.int(0, n))) NULL else offsetx, zero = if (identical(offsetz,
																										      rep.int(0, n))) NULL else offsetz), n = nobs, df.null = nobs -
			     2, df.residual = nobs - (kx + kz + (dist == "negbin")),
			   terms = list(count = mtX, zero = mtZ, full = mt), theta = theta,
			   SE.logtheta = SE.logtheta, loglik = fit$value, vcov = vc,
			   dist = dist, link = linkstr, linkinv = linkinv, converged = fit$convergence <
			     1, call = cl, formula = ff, levels = .getXlevels(mt,
									      mf), contrasts = list(count = attr(X, "contrasts"),
												    zero = attr(Z, "contrasts")))
	      if (model)
		rval$model <- mf
	      if (y)
		rval$y <- Y
	      if (x)
		rval$x <- list(count = X, zero = Z)
	      class(rval) <- "zeroinfl"
	      return(rval)
	    }

            # end of modified pscl code
            ######################################

            zinb_signif = function(df,
                zinbMod1,
                zinbMod0,
                nbMod1,
                nbMod0, DEBUG = F) {
              suppressMessages(require(pscl))
              suppressMessages(require(MASS))
              melted = df

              # filter out genes that have low saturation across all conditions, since pscl sometimes does not fit params well (resulting in large negative intercepts and high std errors)
              NZpercs = aggregate(melted$cnt,by=list(melted$cond),FUN=function(x) { sum(x>0)/length(x) })
              # if (max(NZpercs$x)<=0.15) { return(c(pval=1,status="low saturation (near-essential) across all conditions, not analyzed")) }

              sums = aggregate(melted$cnt,by=list(melted$cond),FUN=sum)
              # to avoid model failing due to singular condition, add fake counts of 1 to all conds if any cond is all 0s
              if (0 %in% sums[,2]) {
                # print("adding pseudocounts")
                for (i in 1:length(sums[,1])) {
                  subset = melted[melted$cond==sums[i,1],]
                  newvec = subset[1,]
                  newvec$cnt = 1 # note: NZmean and NZperc are copied from last dataset in condition
                  #newvec$cnt = as.integer(mean(subset$cnt))+1 # add the mean for each condition as a pseudocount
                  melted = rbind(melted,newvec) }
              }
              status = "-"
              minCount = min(melted$cnt)
              f1 = ""
              mod1 = tryCatch(
                {
                  if (minCount == 0) {
                    f1 = zinbMod1
                    #mod = zeroinfl(as.formula(zinbMod1),data=melted,dist="negbin")
                    # to do: check whether upper/lower are defined when call optim() above...
                    mod = zeroinfl_sid(as.formula(zinbMod1),data=melted,dist="negbin", method="L-BFGS-B", upper=5, lower=-5) # also consider: trace=T
                    coeffs = summary(mod)$coefficients
                    # [,1] is col of parms, [,2] is col of stderrs, assume Intercept is always first
                    if (coeffs$count[,2][1]>0.5) { status = 'warning: high stderr on Intercept for mod1' }
                    coeffs1 = summary(mod)$coefficients
                    a1 = max(abs(c(coeffs1$count[,1],coeffs1$zero[,1])))
                    a2 = max(coeffs1$count[,2],coeffs1$zero[,2])
                    status = paste("A1", round(a1, 2), "A2", round(a2, 2))
                    mod
                  } else {
                    f1 = nbMod1
                    glm.nb(as.formula(nbMod1),data=melted)
                  }
                },
                error=function(err) {
                  status <<- err$message
                  return(NULL)
                })
              f0 = ""
              mod0 = tryCatch( # null model, independent of conditions
                {
                  if (minCount == 0) {
                    f0 = zinbMod0
                    #mod = zeroinfl(as.formula(zinbMod0),data=melted,dist="negbin")
                    mod = zeroinfl_sid(as.formula(zinbMod0),data=melted,dist="negbin", method="L-BFGS-B", upper=5, lower=-5)
                    coeffs = summary(mod)$coefficients
                    # [,1] is col of parms, [,2] is col of stderrs, assume Intercept is always first
                    #if (coeffs$count[,2][1]>0.5) { status = 'warning: high stderr on Intercept for mod0' }
                    coeffs0 = summary(mod)$coefficients
                    b1 = max(abs(c(coeffs0$count[,1],coeffs0$zero[,1])))
                    b2 = max(coeffs0$count[,2],coeffs0$zero[,2])
                    status = paste(status, "B1", round(b1, 2), "B2", round(b2, 2))
                    mod
                  } else {
                    f0 = nbMod0
                    glm.nb(as.formula(nbMod0),data=melted)
                  }
                },
                error=function(err) {
                  status <<- err$message
                  return(NULL)
                })
              if (DEBUG) {
                  print("Model 1:")
                  print(f1)
                  print(summary(mod1))
                  print("Model 0:")
                  print(f0)
                  print(summary(mod0))
              }

              if (is.null(mod1) | is.null(mod0)) { return (c(1, paste0("Model Error. ", status))) }
              if ((minCount == 0) && (sum(is.na(coef(summary(mod1))$count[,4]))>0)) { return(c(1, "Model has coefs, but some pvals are NAs")) } # rare failure mode - has coefs, but pvals are NA
              df1 = attr(logLik(mod1),"df"); df0 = attr(logLik(mod0),"df") # should be (2*ngroups+1)-3
              pval = pchisq(2*(logLik(mod1)-logLik(mod0)),df=df1-df0,lower.tail=F) # alternatively, could use lrtest()
              # this gives same answer, but I would need to extract the Pvalue...
              #require(lmtest)
              #print(lrtest(mod1,mod0))
              tryCatch(
               if (DEBUG) {
                coeffs1 = summary(mod1)$coefficients
                a1 = max(abs(c(coeffs1$count[,1],coeffs1$zero[,1])))
                a2 = max(coeffs1$count[,2],coeffs1$zero[,2])
                coeffs0 = summary(mod0)$coefficients
                b1 = max(abs(c(coeffs0$count[,1],coeffs0$zero[,1])))
                b2 = max(coeffs0$count[,2],coeffs0$zero[,2])
                print(sprintf("coefficients: max(abs(mod1))=%0.2f, max(stderr(mod1)=%0.2f, max(abs(mod0))=%0.2f, max(stderr(mod0))=%0.2f", a1,a2,b1,b2)) }
               ,error=function(err) {
                  status <<- err$message
                  return(NULL)
               })
              return (c(pval, status))
            }
        ''')
        return globalenv['zinb_signif']

    def def_r_zinb_signif_tom(self):
        r('''
            zinb_signif = function(df,
                zinbMod1,
                zinbMod0,
                nbMod1,
                nbMod0, DEBUG = F) {
              suppressMessages(require(pscl))
              suppressMessages(require(MASS))
              melted = df

              # filter out genes that have low saturation across all conditions, since pscl sometimes does not fit params well (resulting in large negative intercepts and high std errors)
              NZpercs = aggregate(melted$cnt,by=list(melted$cond),FUN=function(x) { sum(x>0)/length(x) })
              if (max(NZpercs$x)<=0.15) { return(c(pval=1,status="low saturation (near-essential) across all conditions, not analyzed")) }

              sums = aggregate(melted$cnt,by=list(melted$cond),FUN=sum)
              # to avoid model failing due to singular condition, add fake counts of 1 to all conds if any cond is all 0s
              if (0 %in% sums[,2]) {
                # print("adding pseudocounts")
                for (i in 1:length(sums[,1])) {
                  subset = melted[melted$cond==sums[i,1],]
                  newvec = subset[1,]
                  newvec$cnt = 1 # note: NZmean and NZperc are copied from last dataset in condition
                  #newvec$cnt = as.integer(mean(subset$cnt))+1 # add the mean for each condition as a pseudocount
                  melted = rbind(melted,newvec) }
              }
              status = "-"
              minCount = min(melted$cnt)
              f1 = ""
              mod1 = tryCatch(
                {
                  if (minCount == 0) {
                    f1 = zinbMod1
                    #mod = zeroinfl(as.formula(zinbMod1),data=melted,dist="negbin")
                    #mod = zeroinfl(as.formula(zinbMod1),data=melted,dist="negbin",EM=TRUE)
                    #mod = zeroinfl(as.formula(zinbMod1),data=melted,dist="negbin",start=list(count=c(-1,0,0,0,0,0),zero=c(-1,0,0,0,0,0)))
                    mod = zeroinfl(as.formula(zinbMod1),data=melted,dist="negbin",control=list(method="L-BFGS-B",trace=T))#,lower=rep(-3,13),upper=rep(3,13)))
                    coeffs = summary(mod)$coefficients
                    # [,1] is col of parms, [,2] is col of stderrs, assume Intercept is always first
                    if (coeffs$count[,2][1]>0.5) { status = 'warning: high stderr on Intercept for mod1' }
                    mod
                  } else {
                    f1 = nbMod1
                    glm.nb(as.formula(nbMod1),data=melted)
                  }
                },
                error=function(err) {
                  status <<- err$message
                  return(NULL)
                })
              f0 = ""
              mod0 = tryCatch( # null model, independent of conditions
                {
                  if (minCount == 0) {
                    f0 = zinbMod0
                    mod = zeroinfl(as.formula(zinbMod0),data=melted,dist="negbin")
                    coeffs = summary(mod)$coefficients
                    # [,1] is col of parms, [,2] is col of stderrs, assume Intercept is always first
                    #if (coeffs$count[,2][1]>0.5) { status = 'warning: high stderr on Intercept for mod0' }
                    mod
                  } else {
                    f0 = nbMod0
                    glm.nb(as.formula(nbMod0),data=melted)
                  }
                },
                error=function(err) {
                  status <<- err$message
                  return(NULL)
                })
              if (DEBUG) {
                  print("Model 1:")
                  print(f1)
                  print(summary(mod1))
                  print("Model 0:")
                  print(f0)
                  print(summary(mod0))
              }

              if (is.null(mod1) | is.null(mod0)) { return (c(1, paste0("Model Error. ", status))) }
              if ((minCount == 0) && (sum(is.na(coef(summary(mod1))$count[,4]))>0)) { return(c(1, "Has Coefs, pvals are NAs")) } # rare failure mode - has coefs, but pvals are NA
              df1 = attr(logLik(mod1),"df"); df0 = attr(logLik(mod0),"df") # should be (2*ngroups+1)-3
              pval = pchisq(2*(logLik(mod1)-logLik(mod0)),df=df1-df0,lower.tail=F) # alternatively, could use lrtest()
              # this gives same answer, but I would need to extract the Pvalue...
              #require(lmtest)
              #print(lrtest(mod1,mod0))
              return (c(pval, status))
            }
        ''')

        return globalenv['zinb_signif']

    def winsorize(self, data):
        unique_counts = numpy.unique(numpy.concatenate(data))
        if (len(unique_counts) < 2):
            return data
        else:
            n, n_minus_1 = unique_counts[heapq.nlargest(2, range(len(unique_counts)), unique_counts.take)]
            result = [[ n_minus_1 if count == n else count
                        for count in wig] for wig in data]
        return numpy.array(result)

    def is_number(self, s):
        try:
            float(s)
            return True
        except ValueError:
            return False

    def run_zinb(self, data, genes, NZMeanByRep, LogZPercByRep, RvSiteindexesMap, conditions, covariates, interactions):
        """
            Runs Zinb for each gene across conditions and returns p and q values
            ([[Wigdata]], [Gene], [Number], [Number], {Rv: [SiteIndex]}, [Condition], [Covar], [Interaction]) -> Tuple([Number], [Number], [Status])
            Wigdata :: [Number]
            Gene :: {start, end, rv, gene, strand}
            SiteIndex: Integer
            Condition :: String
            Covar :: String
            Interaction :: String
            Status :: String
        """

        count = 0
        self.progress_range(len(genes))
        pvals,Rvs, status = [],[], []
        r_zinb_signif = self.def_r_zinb_signif()
        if (self.winz):
            self.transit_message("Winsorizing and running analysis...")

        self.transit_message("Condition: %s" % self.condition)

        comp1a = "1+cond"
        comp1b = "1+cond"

        # include cond in mod0 only if testing interactions
        comp0a = "1" if len(self.interactions)==0 else "1+cond"
        comp0b = "1" if len(self.interactions)==0 else "1+cond"
        for I in self.interactions: comp1a += "*"+I; comp1b += "*"+I; comp0a += "+"+I; comp0b += "+"+I
        for C in self.covars: comp1a += "+"+C; comp1b += "+"+C; comp0a += "+"+C; comp0b += "+"+C
        zinbMod1 = "cnt~%s+offset(log(NZmean))|%s+offset(logitZperc)" % (comp1a,comp1b)
        zinbMod0 = "cnt~%s+offset(log(NZmean))|%s+offset(logitZperc)" % (comp0a,comp0b)

        nbMod1 = "cnt~%s" % (comp1a)
        nbMod0 = "cnt~%s" % (comp0a)
        toRFloatOrStrVec = lambda xs: FloatVector([float(x) for x in xs]) if self.is_number(xs[0]) else StrVector(xs)

        for gene in genes:
            count += 1
            Rv = gene["rv"]
            ## Single gene case for debugging
            if (GENE):
                Rv = None
                if GENE in RvSiteindexesMap:
                    Rv = GENE
                else:
                    for g in genes:
                        if (g['gene'] == GENE):
                            Rv = g["rv"]
                            break
                if not Rv:
                    self.transit_error("Cannot find gene: {0}".format(GENE))
                    sys.exit(0)

            if (len(RvSiteindexesMap[Rv]) <= 1):
                status.append("TA sites <= 1")
                pvals.append(1)
            else:
                # For winsorization
                # norm_data = self.winsorize((map(
                #     lambda wigData: wigData[RvSiteindexesMap[Rv]], data))) if self.winz else list(map(lambda wigData: wigData[RvSiteindexesMap[Rv]], data))
                norm_data = list(map(lambda wigData: wigData[RvSiteindexesMap[Rv]], data))
                ([ readCounts,
                   condition,
                   covarsData,
                   interactionsData,
                   NZmean,
                   logitZPerc]) = self.melt_data(
                           norm_data,
                           conditions, covariates, interactions, NZMeanByRep, LogZPercByRep)
                if (numpy.sum(readCounts) == 0):
                    status.append("No counts in all conditions")
                    pvals.append(1)
                else:
                    df_args = {
                        'cnt': IntVector(readCounts),
                        'cond': toRFloatOrStrVec(condition),
                        'NZmean': FloatVector(NZmean),
                        'logitZperc': FloatVector(logitZPerc)
                        }
                    ## Add columns for covariates and interactions if they exist.
                    df_args.update(list(map(lambda t_ic: (t_ic[1], toRFloatOrStrVec(covarsData[t_ic[0]])), enumerate(self.covars))))
                    df_args.update(list(map(lambda t_ic: (t_ic[1], toRFloatOrStrVec(interactionsData[t_ic[0]])), enumerate(self.interactions))))

                    melted = DataFrame(df_args)
                    # r_args = [IntVector(readCounts), StrVector(condition), melted, map(lambda x: StrVector(x), covars), FloatVector(NZmean), FloatVector(logitZPerc)] + [True]
                    debugFlag = True if DEBUG or GENE else False
                    pval, msg = r_zinb_signif(melted, zinbMod1, zinbMod0, nbMod1, nbMod0, debugFlag)
                    status.append(msg)
                    pvals.append(float(pval))
                if (DEBUG or GENE):
                    self.transit_message("Pval for Gene {0}: {1}, status: {2}".format(Rv, pvals[-1], status[-1]))
                if (GENE):
                    self.transit_message("Ran for single gene. Exiting...")
                    sys.exit(0)
            Rvs.append(Rv)
            # Update progress
            text = "Running ZINB Method... %5.1f%%" % (100.0*count/len(genes))
            self.progress_update(text, count)

        pvals = numpy.array(pvals)
        mask = numpy.isfinite(pvals)
        qvals = numpy.full(pvals.shape, numpy.nan)
        qvals[mask] = statsmodels.stats.multitest.fdrcorrection(pvals)[1] # BH, alpha=0.05

        p,q,statusMap = {},{},{}
        for i,rv in enumerate(Rvs):
            p[rv],q[rv],statusMap[rv] = pvals[i],qvals[i],status[i]
        return (p, q, statusMap)

    def Run(self):
        self.transit_message("Starting ZINB analysis")
        start_time = time.time()
        packnames = ("MASS", "pscl")
        r_packages_needed = [x for x in packnames if not rpackages.isinstalled(x)]
        if (len(r_packages_needed) > 0):
            self.transit_error(
                    "Error: Following R packages are required: %(0)s. From R console, You can install them using install.packages(c(%(0)s))"
                    % ({'0': '"{0}"'.format('", "'.join(r_packages_needed))}))
            sys.exit(1)


        self.transit_message("Getting Data")
        (sites, data, filenamesInCombWig) = tnseq_tools.read_combined_wig(self.combined_wig)

        self.transit_message("Normalizing using: %s" % self.normalization)
        (data, factors) = norm_tools.normalize_data(data, self.normalization)

        condition_name = self.condition
        conditionsByFile, covariatesByFileList, interactionsByFileList, orderingMetadata = tnseq_tools.read_samples_metadata(self.metadata, self.covars, self.interactions, condition_name=condition_name)

        ## [Condition] in the order of files in combined wig
        conditions = self.wigs_to_conditions(
            conditionsByFile,
            filenamesInCombWig)
        ## [Covariate] in the order of files in combined wig
        covariates = self.wigs_to_covariates(
            covariatesByFileList,
            filenamesInCombWig)
        ## [Interaction] in the order of files in combined wig
        interactions = self.wigs_to_interactions(
            interactionsByFileList,
            filenamesInCombWig)
        data, conditions, covariates, interactions = self.filter_wigs_by_conditions(
                data,
                conditions,
                covariates = covariates,
                interactions = interactions,
                ignored_conditions = self.ignored_conditions,
                included_conditions = self.included_conditions)

        genes = tnseq_tools.read_genes(self.annotation_path)

        TASiteindexMap = {TA: i for i, TA in enumerate(sites)}
        RvSiteindexesMap = tnseq_tools.rv_siteindexes_map(genes, TASiteindexMap, nterm=self.NTerminus, cterm=self.CTerminus)
        # statsByRv, statGroupNames = self.stats_by_rv(data, RvSiteindexesMap, genes, conditions, interactions)
        LogZPercByRep, NZMeanByRep = self.global_stats_for_rep(data)

        self.transit_message("Running ZINB")
        pvals, qvals, run_status = self.run_zinb(data, genes, NZMeanByRep, LogZPercByRep, RvSiteindexesMap, conditions, covariates, interactions)

        def orderStats(x, y):
            ic1 = x.split(SEPARATOR)
            ic2 = y.split(SEPARATOR)
            c1, i1 = (ic1[0], ic1[1]) if len(ic1) > 1 else (ic1[0], None)
            c2, i2 = (ic2[0], ic2[1]) if len(ic2) > 1 else (ic2[0], None)

            if len(self.included_conditions) > 0:
                condDiff = (self.included_conditions.index(c1) - self.included_conditions.index(c2))
                ## Order by interaction, if stat belongs to same condition
                if condDiff == 0 and i1 is not None and i2 is not None:
                    return (orderingMetadata['interaction'].index(i1) - orderingMetadata['interaction'].index(i2))
                return condDiff

            ## Order by samples metadata, if include flag not provided.
            condDiff = (orderingMetadata['condition'].index(c1) - orderingMetadata['condition'].index(c2))
            if condDiff == 0 and i1 is not None and i2 is not None:
                return (orderingMetadata['interaction'].index(i1) - orderingMetadata['interaction'].index(i2))
            return condDiff

        # orderedStatGroupNames = sorted(statGroupNames, key=functools.cmp_to_key(orderStats))
        # headersStatGroupNames = [x.replace(SEPARATOR,'_') for x in orderedStatGroupNames]

        self.transit_message("Adding File: %s" % (self.output))
        file = open(self.output,"w")
        head = ("Rv Gene TAs".split() +
                # list(map(lambda v: "Mean_" + v, headersStatGroupNames)) +
                # list(map(lambda v: "NZmean_" + v, headersStatGroupNames)) +
                # list(map(lambda v: "NZperc_" + v, headersStatGroupNames)) +
                "pval padj".split() + ["status"])

        file.write("#Console: python %s\n" % " ".join(sys.argv))
        file.write('\t'.join(head)+EOL)
        for gene in genes:
            Rv = gene["rv"]
            vals = ([Rv, gene["gene"], str(len(RvSiteindexesMap[Rv]))] +
                    # ["%0.2f" % statsByRv[Rv]['mean'][group] for group in orderedStatGroupNames] +
                    # ["%0.2f" % statsByRv[Rv]['nz_mean'][group] for group in orderedStatGroupNames] +
                    # ["%0.2f" % statsByRv[Rv]['nz_perc'][group] for group in orderedStatGroupNames] +
                    ["%f" % x for x in [pvals[Rv], qvals[Rv]]]) + [run_status[Rv]]
            file.write('\t'.join(vals)+EOL)
        file.close()
        self.transit_message("Finished Zinb analysis")
        self.transit_message("Time: %0.1fs\n" % (time.time() - start_time))

    @classmethod
    def usage_string(self):
        return """python %s zinb <combined wig file> <samples_metadata file> <annotation .prot_table> <output file> [Optional Arguments]

        Optional Arguments:
        -n <string>         :=  Normalization method. Default: -n TTR
        --ignore-conditions <cond1,cond2> :=  Comma separated list of conditions to ignore, for the analysis.
        --include-conditions <cond1,cond2> :=  Comma separated list of conditions to include, for the analysis. Conditions not in this list, will be ignored.
        -iN <float>     :=  Ignore TAs occuring within given percentage (as integer) of the N terminus. Default: -iN 5
        -iC <float>     :=  Ignore TAs occuring within given percentage (as integer) of the C terminus. Default: -iC 5
        --condition     :=  columnname (in samples_metadata) to use as the Condition. Default: "Condition"
        --covars <covar1,covar2...>     :=  Comma separated list of covariates (in metadata file) to include, for the analysis.
        --interactions <covar1,covar2...>     :=  Comma separated list of covariates to include, that interact with the condition for the analysis. Must be factors
        --gene <RV number or Gene name> := Run method for one gene and print model output.

        """ % (sys.argv[0])

if __name__ == "__main__":
    main()

