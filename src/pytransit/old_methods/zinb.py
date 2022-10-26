from pytransit.components.parameter_panel import panel, progress_update
import time
import sys
import collections
import functools
import heapq

import numpy

from pytransit.old_methods import analysis_base as base
from pytransit.specific_tools.transit_tools import HAS_R
from pytransit.specific_tools import transit_tools, tnseq_tools, norm_tools, logging
from pytransit.generic_tools import misc

############# GUI ELEMENTS ##################

short_name = "ZINB"
long_name = "ZINB"
short_desc = "Perform ZINB analysis"
long_desc = """Perform ZINB analysis"""
GENE = None
DEBUG = True

transposons = ["", ""]
columns = []


class Analysis(base.TransitAnalysis):
    def __init__(self):
        base.TransitAnalysis.__init__(
            self, short_name, long_name, short_desc, long_desc, transposons, ZinbMethod
        )


def main():
    print("ZINB example")


class Method(base.MultiConditionMethod):
    """
    Zinb
    """

    def __init__(
        self,
        covars=[],
        interactions=[],
        pseudocount=1,
        refs=[],
        prot_table=None,
    ):
        base.MultiConditionMethod.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            combined_wig_path,
            metadata_path,
            annotation_path,
            output_path,
            normalization=normalization,
            excluded_conditions=excluded_conditions,
            included_conditions=included_conditions,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
        )
        self.winz = winz
        self.covars = covars
        self.interactions = interactions
        self.group_by = group_by
        self.pseudocount = pseudocount
        self.refs = refs

        if prot_table == None:
            self.gene_name_to_description = None
        else:
            self.gene_name_to_description = {}
            with open(prot_table) as file:
                for line in file:
                    w = line.rstrip().split("\t")
                    rv, descr = w[8], w[0]
                    self.gene_name_to_description[rv] = descr

    @classmethod
    def from_args(self, args, kwargs):
        global DEBUG

        if kwargs.get("-help", False) or kwargs.get("h", False):
            print(ZinbMethod.usage_string)
            sys.exit(0)

        if kwargs.get("v", False):
            DEBUG = True

        if kwargs.get("-gene", False):
            global GENE
            GENE = kwargs.get("-gene", None)

        combined_wig_path = args[0]
        metadata_path = args[1]
        annotation_path = args[2]
        output_path = args[3]
        normalization = kwargs.get("n", "TTR")
        n_terminus = float(kwargs.get("iN", 5.0))
        c_terminus = float(kwargs.get("iC", 5.0))
        pseudocount = float(kwargs.get("PC", 5.0))
        group_by = kwargs.get("-group-by", kwargs["-condition"])
        covars = list(filter(None, kwargs.get("-covars", "").split(",")))
        interactions = list(filter(None, kwargs.get("-interactions", "").split(",")))
        refs = kwargs.get(
            "-ref", []
        )  # list of condition names to use a reference for calculating LFCs
        if refs != []:
            refs = refs.split(",")
        winz = True if "winz" in kwargs else False
        excluded_conditions = list(
            filter(None, kwargs.get("-exclude-conditions", "").split(","))
        )
        included_conditions = list(
            filter(None, kwargs.get("-include-conditions", "").split(","))
        )
        prot_table = kwargs.get("-prot_table", None)

        # check for unrecognized flags
        flags = "-n --exclude-conditions --include-conditions -iN -iC -PC --condition --covars --interactions --gene --ref --prot_table -winz".split()
        from pytransit.specific_tools import console_tools
        console_tools.handle_unrecognized_flags(
            flags,
            kwargs,
            self.usage_string,
        )

        return self(
            combined_wig_path,
            metadata_path,
            annotation_path,
            normalization,
            output_path,
            excluded_conditions,
            included_conditions,
            winz,
            n_terminus,
            c_terminus,
            condition,
            covars,
            interactions,
            pseudocount,
            refs,
            prot_table=prot_table,
        )

    @staticmethod
    def global_stats_for_rep(data):
        """
        Returns the logit zero percentage and nz_mean for each replicate.
            [[WigData]] -> [[Number] ,[Number]]
        note: these are not winsorized even if temp.winz==True
        """

        logit_zero_perc = []
        nz_mean = []
        for wig in data:
            zero_perc = (wig.size - numpy.count_nonzero(wig)) / float(wig.size)
            logit_zero_perc.append(numpy.log(zero_perc / (1 - zero_perc)))
            nz_mean.append(numpy.mean(wig[numpy.nonzero(wig)]))
        return [numpy.array(logit_zero_perc), numpy.array(nz_mean)]

    @staticmethod
    def melt_data(
        readCountsForRv,
        conditions,
        covariates,
        interactions,
        non_zero_mean_by_replicate,
        log_z_perc_by_replicate,
    ):
        rvSitesLength = len(readCountsForRv[0])
        repeatAndFlatten = lambda xs: numpy.repeat(xs, rvSitesLength)

        return [
            numpy.concatenate(readCountsForRv).astype(int),
            repeatAndFlatten(conditions),
            list(map(repeatAndFlatten, covariates)),
            list(map(repeatAndFlatten, interactions)),
            repeatAndFlatten(non_zero_mean_by_replicate),
            repeatAndFlatten(log_z_perc_by_replicate),
        ]
    
    @staticmethod
    def def_r_zinb_signif():
        from pytransit.specific_tools.transit_tools import r, globalenv
        r(
            """
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
              if (max(NZpercs$x)<0.15) { return(c(pval=1,status="low saturation (<15%) across all conditions (pan-growth-defect) - not analyzed")) }

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
                    mod = zeroinfl(as.formula(zinbMod1),data=melted,dist="negbin")
                    coeffs = summary(mod)$coefficients
                    # [,1] is col of parms, [,2] is col of stderrs, assume Intercept is always first
                    #if (coeffs$count[,2][1]>0.5) { status = 'warning: high stderr on Intercept for mod1' }
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
              if ((minCount == 0) && (sum(is.na(coef(summary(mod1))$count[,4]))>0)) { return(c(1, "Has Coefs, but Pvals are NAs (model failure)")) } # rare failure mode - has coefs, but pvals are NA
              df1 = attr(logLik(mod1),"df"); df0 = attr(logLik(mod0),"df") # should be (2*ngroups+1)-3
              if (DEBUG) print(sprintf("delta_log_likelihood=%f",logLik(mod1)-logLik(mod0)))
              pval = pchisq(2*(logLik(mod1)-logLik(mod0)),df=df1-df0,lower.tail=F) # alternatively, could use lrtest()
              # this gives same answer, but I would need to extract the Pvalue...
              #require(lmtest)
              #print(lrtest(mod1,mod0))
              return (c(pval, status))
            }
        """
        )

        return globalenv["zinb_signif"]
    
    @staticmethod
    def calculate(
            *,
            group_by,
            interactions,
            covars,
            winz,
            
            data,
            genes,
            non_zero_mean_by_replicate,
            log_z_perc_by_replicate,
            rv_site_indexes_map,
            filtered_conditions_from_comwig,
            filtered_covariates_from_comwig,
            interactions_from_wigs,
        ):
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
            global DEBUG
            from pytransit.specific_tools.transit_tools import DataFrame, IntVector, FloatVector, StrVector
            import statsmodels.stats.multitest

            count = 0
            
            pvals, Rvs, status = [], [], []
            r_zinb_signif = Method.def_r_zinb_signif()
            logging.log("Running analysis...")
            if winz:
                logging.log("Winsorizing insertion count data")

            logging.log(f"Grouping by: {group_by}")

            comp1a = "1+cond"
            comp1b = "1+cond"

            # include cond in mod0 only if testing interactions
            comp0a = "1" if len(interactions) == 0 else "1+cond"
            comp0b = "1" if len(interactions) == 0 else "1+cond"
            for I in interactions:
                comp1a += "*" + I
                comp1b += "*" + I
                comp0a += "+" + I
                comp0b += "+" + I
            for C in covars:
                comp1a += "+" + C
                comp1b += "+" + C
                comp0a += "+" + C
                comp0b += "+" + C
            zinbMod1 = "cnt~%s+offset(log(NZmean))|%s+offset(logitZperc)" % (comp1a, comp1b)
            zinbMod0 = "cnt~%s+offset(log(NZmean))|%s+offset(logitZperc)" % (comp0a, comp0b)

            nbMod1 = "cnt~%s" % (comp1a)
            nbMod0 = "cnt~%s" % (comp0a)
            toRFloatOrStrVec = (
                lambda xs: FloatVector([float(x) for x in xs])
                if misc.str_is_float(xs[0])
                else StrVector(xs)
            )

            for gene in genes:
                count += 1
                Rv = gene["rv"]
                ## Single gene case for debugging
                if GENE:
                    Rv = None
                    if GENE in rv_site_indexes_map:
                        Rv = GENE
                    else:
                        for g in genes:
                            if g["gene"] == GENE:
                                Rv = g["rv"]
                                break
                    if not Rv:
                        logging.error("Cannot find gene: {0}".format(GENE))

                if DEBUG:
                    logging.log("======================================================================")
                    logging.log(gene["rv"] + " " + gene["gene"])

                if len(rv_site_indexes_map[Rv]) <= 1:
                    status.append("TA sites <= 1, not analyzed")
                    pvals.append(1)
                else:
                    norm_data = list(
                        map(lambda wigData: wigData[rv_site_indexes_map[Rv]], data)
                    )
                    if winz:
                        norm_data = tnseq_tools.winsorize(norm_data)
                    (
                        [
                            readCounts,
                            condition,
                            covarsData,
                            interactionsData,
                            NZmean,
                            logitZPerc,
                        ]
                    ) = Method.melt_data(
                        readCountsForRv=norm_data,
                        conditions=filtered_conditions_from_comwig,
                        covariates=filtered_covariates_from_comwig,
                        interactions=interactions_from_wigs,
                        non_zero_mean_by_replicate=non_zero_mean_by_replicate,
                        log_z_perc_by_replicate=log_z_perc_by_replicate,
                    )
                    if numpy.sum(readCounts) == 0:
                        status.append(
                            "pan-essential (no counts in all conditions) - not analyzed"
                        )
                        pvals.append(1)
                    else:
                        df_args = {
                            "cnt": IntVector(readCounts),
                            "cond": toRFloatOrStrVec(condition),
                            "NZmean": FloatVector(NZmean),
                            "logitZperc": FloatVector(logitZPerc),
                        }
                        ## Add columns for covariates and interactions if they exist.
                        df_args.update(
                            list(
                                map(
                                    lambda t_ic: (
                                        t_ic[1],
                                        toRFloatOrStrVec(covarsData[t_ic[0]]),
                                    ),
                                    enumerate(covars),
                                )
                            )
                        )
                        df_args.update(
                            list(
                                map(
                                    lambda t_ic: (
                                        t_ic[1],
                                        toRFloatOrStrVec(interactionsData[t_ic[0]]),
                                    ),
                                    enumerate(interactions),
                                )
                            )
                        )

                        melted = DataFrame(df_args)
                        # r_args = [IntVector(readCounts), StrVector(condition), melted, map(lambda x: StrVector(x), covars), FloatVector(NZmean), FloatVector(logitZPerc)] + [True]
                        debugFlag = True if DEBUG or GENE else False
                        pval, msg = r_zinb_signif(
                            melted, zinbMod1, zinbMod0, nbMod1, nbMod0, debugFlag
                        )
                        status.append(msg)
                        pvals.append(float(pval))
                    if DEBUG or GENE:
                        logging.log(
                            "Pval for Gene {0}: {1}, status: {2}".format(
                                Rv, pvals[-1], status[-1]
                            )
                        )
                    if GENE:
                        logging.log("Ran for single gene. Exiting...")
                        sys.exit(0)
                Rvs.append(Rv)
                # Update progress
                percentage = (100.0 * count / len(genes))
                text = "Running ZINB Method... %5.1f%%" % percentage
                progress_update(text, percentage)

            pvals = numpy.array(pvals)
            mask = numpy.isfinite(pvals)
            qvals = numpy.full(pvals.shape, numpy.nan)
            qvals[mask] = statsmodels.stats.multitest.fdrcorrection(pvals)[
                1
            ]  # BH, alpha=0.05

            p, q, statusMap = {}, {}, {}
            for i, rv in enumerate(Rvs):
                p[rv], q[rv], statusMap[rv] = pvals[i], qvals[i], status[i]
        return (p, q, statusMap)

    # pairs is a list of (var,val); samples is a set; varsByFileList is a list of dictionaries mapping values to samples for each var (parallel to vars)
    # recursive: keep calling till vars reduced to empty

    def Run(self):
        transit_tools.require_r_to_be_installed()
        from pytransit.specific_tools.transit_tools import EOL, SEPARATOR, rpackages
        logging.log("Starting ZINB analysis")
        start_time = time.time()
        packnames = ("MASS", "pscl")
        r_packages_needed = [x for x in packnames if not rpackages.isinstalled(x)]
        if len(r_packages_needed) > 0:
            logging.error(
                "Error: Following R packages are required: %(0)s. From R console, You can install them using install.packages(c(%(0)s))"
                % ({"0": '"{0}"'.format('", "'.join(r_packages_needed))})
            )
    
        logging.log("Getting Data")
        (sites, data, filenames_in_comb_wig) = tnseq_tools.CombinedWigData.load(
            self.combined_wig_path
        )

        logging.log("Normalizing using: %s" % self.normalization)
        (data, factors) = norm_tools.normalize_data(data, self.normalization)

        condition_name = self.group_by
        # if a covar is not found, this crashes; check for it?
        # read it first with no condition specified, to get original Condition names
        (
            conditions_by_file1,
            covariatesByFileList1,
            interactionsByFileList1,
            orderingMetadata1,
        ) = tnseq_tools.CombinedWigMetadata.read_condition_data(
            self.metadata_path, self.covars, self.interactions
        )  # without specifiying condition
        (
            conditions_by_file,
            covariatesByFileList,
            interactionsByFileList,
            orderingMetadata,
        ) = tnseq_tools.CombinedWigMetadata.read_condition_data(
            self.metadata_path, self.covars, self.interactions, condition_name=condition_name
        )
        
        # 
        ## [Condition] in the order of files in combined wig
        # 
        conditions_from_comwig = [
            conditions_by_file.get(f, "FLAG-UNMAPPED-CONDITION-IN-WIG") for f in filenames_in_comb_wig
        ]
        # 
        ## [Covariate] in the order of files in combined wig
        # 
        try:
            covariates_from_comwig = [
                [covarsByFile.get(f, "?") for f in filenames_in_comb_wig]
                for covarsByFile in covariatesByFileList
            ]
        except KeyError:
            for f in filenames_in_comb_wig:
                if f not in covariatesByFileList[0]:
                    print(f)
            logging.error("Error: Covariates not found for sample:")
        # 
        ## [Interaction] in the order of files in combined wig
        # 
        try:
            interactions_from_wigs = [
                [covarsByFile.get(f, "?") for f in filenames_in_comb_wig]
                for covarsByFile in interactionsByFileList
            ]
        except KeyError:
            for f in filenames_in_comb_wig:
                if f not in interactionsByFileList[0]:
                    print(f)
            logging.error("Error: Interaction var not found for sample")

        metadata = transit_tools.get_samples_metadata(metadata_path)
        conditionNames = metadata[
            "Condition"
        ]  # original Condition names for each sample, as ordered list
        fileNames = metadata["Filename"]

        # this is the new way to filter samples, where --include/exclude-conditions refers to Condition column in metadata, regardless of whether --condition was specified
        (
            data,
            fileNames,
            conditionNames,
            filtered_conditions_from_comwig,
            filtered_covariates_from_comwig,
            interactions_from_wigs,
        ) = transit_tools.filter_wigs_by_conditions2(
            data,
            wig_fingerprints=fileNames,
            condition_names=conditionNames,  # original Condition column in samples metadata file
            included_cond=self.included_conditions,
            excluded_cond=self.excluded_conditions,
            conditions=conditions_from_comwig,  # user might have specified a column other than Condition
            covariates=covariates_from_comwig,
            interactions=interactions_from_wigs,
        )
        samples_used = set(fileNames)

        # show the samples associated with each condition (and covariates or interactions, if defined), and count samples in each cross-product of vars

        vars = [condition_name] + self.covars + self.interactions
        vars2vals = {}
        vars2vals[condition_name] = misc.no_duplicates(filtered_conditions_from_comwig)
        for i, var in enumerate(self.covars):
            vars2vals[var] = misc.no_duplicates(filtered_covariates_from_comwig[i])
        for i, var in enumerate(self.interactions):
            vars2vals[var] = list(set(interactions_from_wigs[i]))
        varsByFileList = (
            [conditions_by_file] + covariatesByFileList + interactionsByFileList
        )
        for i, var in enumerate(vars):
            print("\nsamples for Condition/Covariate/Interaction: %s" % vars[i])
            filesByVar = misc.invert_dict(varsByFileList[i])
            for k, v in filesByVar.items():
                samples = list(samples_used.intersection(set(v)))
                if k in vars2vals.get(var, []):
                    print("%s: %s" % (k, " ".join(samples)))
        pairs = []
        print("\nsamples in cross-product:")
        any_empty = transit_tools.expand_var(
            [], vars, varsByFileList, vars2vals, set(samples_used)
        )
        if any_empty:
            print(
                "warning: ZINB requires samples in all combinations of conditions; the fact that one is empty could result in Model Errors"
            )

        genes = tnseq_tools.read_genes(self.annotation_path)

        ta_site_index_map = {TA: i for i, TA in enumerate(sites)}
        rv_site_indexes_map = tnseq_tools.rv_siteindexes_map(
            genes, ta_site_index_map, n_terminus=self.n_terminus, c_terminus=self.c_terminus
        )
        stats_by_rv, stat_group_names = transit_tools.get_stats_by_rv(
            data, rv_site_indexes_map, genes, filtered_conditions_from_comwig, interactions_from_wigs, self.winz
        )
        log_z_perc_by_replicate, non_zero_mean_by_replicate = Method.global_stats_for_rep(data)

        logging.log("Running ZINB")
        pvals, qvals, run_status = Method.calculate(
            group_by=self.group_by,
            interactions=self.interactions,
            covars=self.covars,
            winz=self.winz,
            data=data,
            genes=genes,
            non_zero_mean_by_replicate=non_zero_mean_by_replicate,
            log_z_perc_by_replicate=log_z_perc_by_replicate,
            rv_site_indexes_map=rv_site_indexes_map,
            filtered_conditions_from_comwig=filtered_conditions_from_comwig,
            filtered_covariates_from_comwig=filtered_covariates_from_comwig,
            interactions_from_wigs=interactions_from_wigs,
        )

        def order_stats(x, y):
            ic1 = x.split(SEPARATOR)
            ic2 = y.split(SEPARATOR)
            c1, i1 = (ic1[0], ic1[1]) if len(ic1) > 1 else (ic1[0], None)
            c2, i2 = (ic2[0], ic2[1]) if len(ic2) > 1 else (ic2[0], None)

            # use --include-conditions to determine order of columns in output file
            # this only works if an alternative --condition was not specified
            # otherwise don't try to order them this way because it gets too complicated
            # possibly should require --covars and --interactions to be unspecified too
            if (
                self.group_by == "Condition"
                and len(self.included_conditions) > 0
                and len(self.excluded_conditions) == 0
            ):
                condDiff = self.included_conditions.index(
                    c1
                ) - self.included_conditions.index(c2)
                ## Order by interaction, if stat belongs to same condition
                if condDiff == 0 and i1 is not None and i2 is not None:
                    return orderingMetadata["interaction"].index(i1) - orderingMetadata[
                        "interaction"
                    ].index(i2)
                return condDiff

            ## Order by samples metadata, if include flag not provided.
            condDiff = orderingMetadata["condition"].index(c1) - orderingMetadata[
                "condition"
            ].index(c2)
            if condDiff == 0 and i1 is not None and i2 is not None:
                return orderingMetadata["interaction"].index(i1) - orderingMetadata[
                    "interaction"
                ].index(i2)
            return condDiff

        orderedStatGroupNames = sorted(
            stat_group_names, key=functools.cmp_to_key(order_stats)
        )
        headersStatGroupNames = [
            x.replace(SEPARATOR, "_") for x in orderedStatGroupNames
        ]

        logging.log("Adding File: %s" % (self.output_path))
        file = open(self.output_path, "w")
        if len(headersStatGroupNames) == 2:
            lfcNames = ["LFC"]
        else:
            lfcNames = list(map(lambda v: "LFC_" + v, headersStatGroupNames))
        head = (
            "Rv Gene TAs".split()
            + list(map(lambda v: "Mean_" + v, headersStatGroupNames))
            + lfcNames
            + list(map(lambda v: "NZmean_" + v, headersStatGroupNames))
            + list(map(lambda v: "NZperc_" + v, headersStatGroupNames))
            + "pval padj".split()
            + ["status"]
        )

        file.write("#Console: python3 %s\n" % " ".join(sys.argv))
        file.write(
            "#parameters: normalization=%s, trimming=%s/%s%% (N/C), pseudocounts=%s\n"
            % (self.normalization, self.n_terminus, self.c_terminus, self.pseudocount)
        )
        if self.gene_name_to_description != None:
            head.append("annotation")
        file.write("#" + "\t".join(head) + EOL)
        for gene in genes:
            Rv = gene["rv"]
            means = [stats_by_rv[Rv]["mean"][group] for group in orderedStatGroupNames]
            pseudocount = self.pseudocount
            if len(means) == 2:
                LFCs = [
                    numpy.math.log((means[1] + pseudocount) / (means[0] + pseudocount), 2)
                ]  # still need to adapt this to use --ref if defined
            else:
                if len(self.refs) == 0:
                    m = numpy.mean(means)  # grand mean across all conditions
                else:
                    m = numpy.mean(
                        [stats_by_rv[Rv]["mean"][group] for group in self.refs]
                    )
                LFCs = [numpy.math.log((x + pseudocount) / (m + pseudocount), 2) for x in means]
            vals = (
                [Rv, gene["gene"], str(len(rv_site_indexes_map[Rv]))]
                + [
                    "%0.1f" % stats_by_rv[Rv]["mean"][group]
                    for group in orderedStatGroupNames
                ]
                + ["%0.3f" % x for x in LFCs]
                + [
                    "%0.1f" % stats_by_rv[Rv]["nz_mean"][group]
                    for group in orderedStatGroupNames
                ]
                + [
                    "%0.2f" % stats_by_rv[Rv]["nz_perc"][group]
                    for group in orderedStatGroupNames
                ]
                + ["%f" % x for x in [pvals[Rv], qvals[Rv]]]
            ) + [run_status[Rv]]
            if self.gene_name_to_description != None:
                vals.append(self.gene_name_to_description.get(Rv, "?"))
            file.write("\t".join(vals) + EOL)
        file.close()
        logging.log("Finished Zinb analysis")
        logging.log("Time: %0.1fs\n" % (time.time() - start_time))

ZinbMethod = Method