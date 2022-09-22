import sys
import os
import time
import ntpath
import math
import random
import datetime
import collections
import heapq

import numpy

import wx.grid
from pytransit.tools import logging, gui_tools
from pytransit.tools.transit_tools import wx
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.parameter_panel import panel,progress_update

from pytransit.tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.basics.lazy_dict import LazyDict
import pytransit.basics.csv as csv
import pytransit.basics.misc as misc
from pytransit.tools.transit_tools import wx, pub, basename, HAS_R, FloatVector, DataFrame, StrVector, EOL
from pytransit.universal_data import universal
from pytransit.components import file_display, results_area, parameter_panel, panel_helpers
from pytransit.tools.gui_tools import nice_error_log

command_name = sys.argv[0]


@misc.singleton
class Analysis:
    name = "PathwayEnrichment" # HANDLE_THIS
    identifier  = name
    short_name  = name.lower()
    long_name   = name.upper()
    short_desc  = f"Perform {name} analysis"
    long_desc   = f"""Perform {name} analysis"""
    rows = []

    inputs = LazyDict(
        resampling_file = None,
        associations_file = None,
        pathways_file = None,
        output_path= None,
    )
    
    valid_cli_flags = [
        "-M", 
        "-Pval_col",
        "-Qval_col",
        "-ranking",
        "-LFC_col",
        "-p",
        "-Nperm",
        "-PC"
    ]
    usage_string = """python3 %s pathway_enrichment <resampling_file> <associations> <pathways> <output_file> [-M <FET|GSEA|GO>] [-PC <int>] [-ranking SLPV|LFC] [-p <float>] [-Nperm <int>] [-Pval_col <int>] [-Qval_col <int>]  [-LFC_col <int>]

        Optional parameters:
        -M FET|GSEA|ONT:     method to use, FET for Fisher's Exact Test (default), GSEA for Gene Set Enrichment Analysis (Subramaniam et al, 2005), or ONT for Ontologizer (Grossman et al, 2007)
        -Pval_col <int>    : indicate column with *raw* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for sorting) (default: -2)
        -Qval_col <int>    : indicate column with *adjusted* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for significant cutoff) (default: -1)
        for GSEA...
        -ranking SLPV|LFC  : SLPV is signed-log-p-value (default); LFC is log2-fold-change from resampling 
        -LFC_col <int>     : indicate column with log2FC (starting with 0; can also be negative, i.e. -1 means last col) (used for ranking genes by SLPV or LFC) (default: 6)
        -p <float>         : exponent to use in calculating enrichment score; recommend trying 0 or 1 (as in Subramaniam et al, 2005)
        -Nperm <int>       : number of permutations to simulate for null distribution to determine p-value (default=10000)
        for FET...
        -PC <int>          :  pseudo-counts to use in calculating p-value based on hypergeometric distribution (default=2)
        """ % sys.argv[0]
    
    
    wxobj = None
    panel = None
    def __init__(self, *args, **kwargs):
        self.instance = self.method = self.gui = self # for compatibility with older code/methods
        self.full_name        = f"[{self.short_name}]  -  {self.short_desc}"
   
    
    def __str__(self):
        return f"""
            Analysis Method:
                Short Name:  {self.short_name}
                Long Name:   {self.long_name}
                Short Desc:  {self.short_desc}
                Long Desc:   {self.long_desc}
        """.replace('\n            ','\n').strip()
    
    def __repr__(self): return f"{self.inputs}"
    def __call__(self): return self

    #####################  need to pull out in other function #########################

    def create_input_field(self, panel, sizer, label, value,tooltip=None):
        get_text = panel_helpers.create_text_box_getter(
            panel,
            sizer,
            label_text=label,
            default_value=value,
            tooltip_text=tooltip,
        )
        return lambda *args: get_text()
    def create_int_field(self, panel, sizer, label, value,tooltip=None):
        get_text = panel_helpers.create_text_box_getter(
            panel,
            sizer,
            label_text=label,
            default_value=value,
            tooltip_text=tooltip,
        )
        return lambda *args: int(get_text())

    def call_from_results_panel(self, results_file):
        self.inputs.resampling_file = results_file
        self.define_panel()


    def define_panel(self,_=None):
        universal.interface ="gui" #without this, selected from Results Panel shows an error
        universal.selected_method = self #without this, selected from Results Panel shows an error
        panel.method_sizer_text.Show() #without this panel below behaves strangely when method selected from resultspanel


        # Ask for resampling file, if it did not come from the 
        if Analysis.inputs.resampling_file==None:
            Analysis.inputs.resampling_file=gui_tools.ask_for_file(
            message = "Select a File for PathWay Analysis",
            allowed_extensions='All files (*.*)|*.*',
            ) 

        self.panel = panel_helpers.make_panel()
        self.value_getters = LazyDict()
        main_sizer = wx.BoxSizer(wx.VERTICAL)


        self.value_getters.associations_file = panel_helpers.create_file_input(self.panel,main_sizer, \
          button_label="Associations file",default_file_name="sanger_assocations.txt",allowed_extensions="All files (*.*)|*.*", \
          popup_title="Choose Associations file", \
          tooltip_text="Must exist to run Pathway enrichment.")

        
        self.value_getters.pathways_file = panel_helpers.create_file_input(self.panel,main_sizer, \
          button_label="Pathways file",default_file_name="sanger_pathways.txt",allowed_extensions="All files (*.*)|*.*", \
          popup_title="Choose Pathways file", \
          tooltip_text="Must exist to run Pathway enrichment.")

        self.value_getters.method = panel_helpers.create_text_box_getter(self.panel, main_sizer, label_text="Method", default_value="FET", tooltip_text="method to use, FET for Fisher's Exact Test (default), GSEA for Gene Set Enrichment Analysis (Subramaniam et al, 2005), or ONT for Ontologizer (Grossman et al, 2007)")
        self.value_getters.pval_col = self.create_int_field(self.panel,main_sizer, value=-2, label="Pval Col", tooltip="indicate column with *raw* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for sorting)")
        self.value_getters.qval_col = self.create_int_field(self.panel,main_sizer, value=-1, label="Qval Col", tooltip="indicate column with *adjusted* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for significant cutoff)")
        self.value_getters.ranking = panel_helpers.create_text_box_getter(self.panel, main_sizer, label_text="ranking", default_value="SPLV", tooltip_text="SLPV is signed-log-p-value (default); LFC is log2-fold-change from resampling")
        self.value_getters.LFC_col = self.create_input_field(self.panel,main_sizer, value=6, label="LFC col", tooltip="indicate column with log2FC (starting with 0; can also be negative, i.e. -1 means last col) (used for ranking genes by SLPV or LFC)")
        self.value_getters.enrichment_exponent =  self.create_input_field(self.panel,main_sizer, label="Enrichment Exponent",value=1,tooltip="exponent to use in calculating enrichment score; recommend trying 0 or 1 (as in Subramaniam et al, 2005)")
        self.value_getters.num_permutations = self.create_input_field(self.panel,main_sizer, label="Number of Permutations",value=10000,tooltip="number of permutations to simulate for null distribution to determine p-value")
        self.value_getters.pseudocount = self.create_int_field(self.panel,main_sizer, label="Pseudocount",value=2,tooltip="pseudo-counts to use in calculating p-value based on hypergeometric distribution")

            
        panel_helpers.create_run_button(self.panel, main_sizer)
        parameter_panel.hide_all_options()    
        parameter_panel.set_panel(self.panel)
        self.panel.SetSizer(main_sizer)
        self.panel.Layout()
        main_sizer.Fit(self.panel)
        universal.frame.Layout()

    @classmethod
    def from_gui(cls,frame):            

        for each_key, each_getter in Analysis.value_getters.items():
            try:
                Analysis.inputs[each_key] = each_getter()
            except Exception as error:
                logging.error(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')

        Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Analysis.short_name}_output.dat",
            output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
        )   
        return Analysis.instance

    @staticmethod
    def from_args(args, kwargs):
        if len(args)!=4: logging.log(Analysis.usage_string); sys.exit(0) # use transit_error()?

        console_tools.handle_help_flag(kwargs, Analysis.usage_string)
        console_tools.handle_unrecognized_flags(Analysis.valid_cli_flags, kwargs, Analysis.usage_string)

        # save the data
        Analysis.inputs.update(dict(
            resampling_file = args[0],
            associations_file = args[1],
            pathways_file = args[2],
            output_path=args[3],
            method = kwargs.get("M", "FET"),
            pval_col = int(kwargs.get("Pval_col", "-2")),
            qval_col = int(kwargs.get("Qval_col", "-1")),
            ranking = kwargs.get("ranking", "SPLV"),
            LFC_col = int(kwargs.get("LFC_col", "6")),
            enrichment_exponent = kwargs.get("p", "1"),
            num_permutations = int(kwargs.get("Nperm", "10000")),
            pseudocount = int(kwargs.get("PC", "2")),
        ))
        
        return Analysis.instance
        
    def Run(self):
        with gui_tools.nice_error_log:
            from pytransit.tools import stat_tools
            logging.log(f"Starting {Analysis.identifier} analysis")
            start_time = time.time()

            
                #checking validation of inputs
            if self.inputs.method == "FET":
                self.fisher_exact_test()
            elif self.inputs.method == "GSEA":
                self.GSEA()
            elif self.inputs.method == "ONT":
                self.Ontologizer()
            else:
                self.inputs.method = "Not a valid method"
                progress_update("Not a valid method", 100)
        
        # 
        # process data
        # 
        # if True:
        #     rows, summary_info = stat_tools.{analysis_name}(**self.inputs) # HANDLE_THIS
        
        # 
        # write output
        # 
            if True:
                logging.log(f"Adding File: {self.inputs.output_path}")
                # 
                # write to file
                # 
                header = "pathway total_genes(M) genes_in_path(n) significant_genes(N) signif_genes_in_path(k) expected k+PC n_adj_by_PC enrichement pval qval description genes"
                transit_tools.write_result(
                    path=self.inputs.output_path, # path=None means write to STDOUT
                    file_kind=Analysis.identifier,
                    rows=self.rows,
                    column_names=header.split(),
                    extra_info=dict(
                        parameters=dict(
                            resampling_file = self.inputs.resampling_file,
                            associations_file = self.inputs.associations_file,
                            pathways_file = self.inputs.pathways_file,
                            output_path=self.inputs.output_path,
                            method = self.inputs.method,
                            pval_col = self.inputs.pval_col,
                            qval_col = self.inputs.qval_col,
                            ranking = self.inputs.ranking,
                            LFC_col = self.inputs.LFC_col,
                            enrichment_exponent = self.inputs.enrichment_exponent,
                            num_permutations = self.inputs.num_permutations,
                            pseudocount = self.inputs.pseudocount,
                        ),
                    ),
                )
                logging.log(f"Finished {Analysis.identifier} analysis in {time.time() - start_time:0.1f}sec")
            results_area.add(self.inputs.output_path)
        
    # def write(self, msg):
    #     self.output.write(msg + "\n")

    def read_resampling_file(self, filename):
        logging.log("Reading in Resampling File", filename)
        genes, hits, headers = [], [], []
        with open(filename) as file:
            for line in file:
                if line[0] == "#":
                    headers.append(line)
                    continue
                w = line.rstrip().split("\t")
                genes.append(w)
                qval = float(w[self.inputs.qval_col])
                if qval < 0.05:
                    hits.append(w[0])
        return genes, hits, headers
        # assume these are listed as pairs (tab-sep)
        # return bidirectional hash (genes->[terms], terms->[genes]; each can be one-to-many, hence lists)
        # filter could be a subset of genes we want to focus on (throw out the rest)

    def read_associations(self, filename, filter=None):
        logging.log("Reading in Associations File")
        associations = {}
        with open(filename) as file:
            for line in file:
                if line[0] == "#":
                    continue
                w = line.rstrip().split("\t")
                if filter != None and w[0] not in filter:
                    continue  # skip genes in association file that are not relevant (i.e. not in resampling file)
                # store mappings in both directions
                for (a, b) in [(w[0], w[1]), (w[1], w[0])]:
                    if a not in associations:
                        associations[a] = []
                    if b not in associations[a]:
                        associations[a].append(b)  # ignore duplicates
        return associations

    def read_pathways(self, filename):
        pathways = {}
        with open(filename) as file:
            for line in file:
                if line[0] == "#":
                    continue
                w = line.rstrip().split("\t")
                pathways[w[0]] = w[1]
        return pathways
     ############### GSEA ######################

    def makeindex(self, lst):
        index = {}
        for i in range(len(lst)):
            index[lst[i]] = i
        return index

    # based on GSEA paper (Subramanian et al, 2005, PNAS)
    # ranks and scores are hashes from genes into ranks and SLPV
    # when p=0, ES(S) reduces to the standard K-S statistic; p=1 is used in PNAS paper

    def enrichment_score(self, A, ranks, scores, p=0):
        n = len(ranks)
        n2 = int(n / 2.0)
        Aranks = [ranks.get(x, n2) for x in A]  # default to middle if not found
        Ascores = [scores.get(x, 0) for x in A]  # default to 0 if not found
        pairs = list(zip(Aranks, Ascores))
        pairs.sort()  # sort A by ranks
        Aranks, Ascores = [x[0] for x in pairs], [x[1] for x in pairs]
        powers = [math.pow(abs(x), p) for x in Ascores]
        NR = sum(powers)
        if NR == 0:
            return 0  # special case
        Nmiss = n - len(A)  # totalGenes-hits
        powersum, best = 0, -1
        for i in range(len(powers)):
            powersum += powers[i]
            Phit = powersum / float(NR)
            Pmiss = (Aranks[i] - i) / float(Nmiss)
            es = abs(Phit - Pmiss)  # looking for max deviation
            if es > best:
                best = es
        return best

    def mean_rank(self, A, orfs2ranks):
        n2 = len(orfs2ranks.keys()) / 2
        return round(numpy.mean([orfs2ranks.get(x, n2) for x in A]), 1)

    # during initialization, self.inputs.resampling_file etc have been set, and self.output has been opened

    def GSEA(self):
        from statsmodels.stats import multitest
        data, hits, headers = self.read_resampling_file(
            self.inputs.resampling_file
        )  # hits are not used in GSEA()
        orfs_in_resampling_file = [w[0] for w in data]
        headers = headers[-1].rstrip().split("\t")  # last line prefixed by '#'
        associations = self.read_associations(
            self.inputs.associations_file, filter=orfs_in_resampling_file
        )  # bidirectional map; includes term->genelist and gene->termlist
        # filter: project associations (of orfs to pathways) onto only those orfs appearing in the resampling file

        ontology = self.read_pathways(self.inputs.pathways_file)
        genenames = {}
        for gene in data:
            genenames[gene[0]] = gene[1]
        n2 = int(len(data) / 2)
        terms = list(ontology.keys())
        terms2orfs = associations
        allgenes = [x[0] for x in data]

        # self.rows.append("# method=GSEA, Nperm=%d, p=%d" % (self.Nperm, self.p))
        # self.rows("# ranking genes by %s" % self.ranking)
        # self.rowst("# total genes: %s, mean rank: %s" % (len(data), n2))

        # rank by SLPV=sign(LFC)*log10(pval)
        # note: genes with lowest p-val AND negative LFC have highest scores (like positive correlation)
        # there could be lots of ties with pval=0 or 1, but so be it (there are probably fewer such ties than with Qvals) (randomize order below)
        pairs = []  # pair are: rv and score (SLPV)
        for w in data:
            orf = w[0]
            if self.inputs.ranking == "SLPV":
                # Pval_col = headers.index("p-value")
                Pval = float(w[self.Pval_col])
                LFC = float(w[self.LFC_col])
                SLPV = (-1 if LFC < 0 else 1) * math.log(Pval + 0.000001, 10)
                pairs.append((orf, SLPV))
            elif self.inputs.ranking == "LFC":
                # LFC_col = headers.index("log2FC")
                LFC = float(w[self.LFC_col])
                pairs.append((orf, LFC))

        # pre-randomize ORFs, to avoid genome-position bias in case of ties in pvals (e.g. 1.0)
        indexes = range(len(pairs))
        indexes = numpy.random.permutation(indexes).tolist()
        pairs = [pairs[i] for i in indexes]

        pairs.sort(
            key=lambda x: x[1], reverse=True
        )  # emulate ranking genes with *higher* correlation at top
        orfs2rank, orfs2score = {}, {}
        for i, (orf, score) in enumerate(pairs):
            orfs2score[orf] = score
            orfs2rank[orf] = i

        Nperm = self.Nperm
        results, Total = [], len(terms)
        for i, term in enumerate(terms):
            sys.stdout.flush()
            orfs = terms2orfs.get(term, [])
            num_genes_in_pathway = len(orfs)
            if num_genes_in_pathway < 2:
                continue  # skip pathways with less than 2 genes
            mr = self.mean_rank(orfs, orfs2rank)
            es = self.enrichment_score(
                orfs, orfs2rank, orfs2score, p=self.p
            )  # always positive, even if negative deviation, since I take abs
            higher = 0
            for n in range(Nperm):
                perm = random.sample(
                    allgenes, num_genes_in_pathway
                )  # compare to enrichment score for random sets of genes of same size
                e2 = self.enrichment_score(perm, orfs2rank, orfs2score, p=self.p)
                if e2 > es:
                    higher += 1
                if n > 100 and higher > 10:
                    break  # adaptive: can stop after seeing 10 events (permutations with higher ES)
            pval = higher / float(n)
            vals = [
                "#",
                term,
                num_genes_in_pathway,
                mr,
                es,
                pval,
                ontology.get(term, "?"),
            ]
            # sys.stderr.write(' '.join([str(x) for x in vals])+'\n')
            percentage = (100.0 * i) / Total
            text = "Running Pathway Enrichment Method... %5.1f%%" % (percentage)
            progress_update(text, percentage)
            results.append((term, mr, es, pval))

        results.sort(key=lambda x: x[1])  # sort on mean rank
        pvals = [x[-1] for x in results]
        rej, qvals = multitest.fdrcorrection(pvals)
        results = [tuple(list(res) + [q]) for res, q in zip(results, qvals)]

        n2 = int(len(data) / 2)
        up, down = 0, 0
        for term, mr, es, pval, qval in results:
            if qval < 0.05:
                if mr < n2:
                    up += 1
                else:
                    down += 1

        # self.rows.append(
        #     "# significant pathways enriched for conditionally ESSENTIAL genes: %s (qval<0.05, mean_rank<%s) (includes genes that are MORE required in condition B than A)"
        #     % (up, n2)
        # )
        for term, mr, es, pval, qval in results:
            if qval < 0.05 and mr < n2:
                self.rows.append(
                    "#   %s %s (mean_rank=%s)" % (term, ontology.get(term, "?"), mr)
                )
        # self.rows.append(
        #     "# significant pathways enriched for conditionally NON-ESSENTIAL genes: %s (qval<0.05, mean_rank>%s) (includes genes that are LESS required in condition B than A)"
        #     % (down, n2)
        # )
        # for term, mr, es, pval, qval in results:
        #     if qval < 0.05 and mr > n2:
        #         self.rows.append(
        #             "#   %s %s (mean_rank=%s)" % (term, ontology.get(term, "?"), mr)
        #         )
        # # self.rows.append("# pathways sorted by mean_rank")

        # self.inputs.output.write(
        #     "\t".join(
        #         "#pathway description num_genes mean_rank GSEA_score pval qval genes".split()
        #     )
        #     + "\n"
        # )
        for term, mr, es, pval, qval in results:
            rvs = terms2orfs[term]
            rvinfo = [(x, genenames.get(x, "?"), orfs2rank.get(x, n2)) for x in rvs]
            rvinfo.sort(key=lambda x: x[2])
            rvs = ["%s/%s (%s)" % x for x in rvinfo]
            rvs = " ".join(rvs)
            vals = (
                [term, ontology.get(term, "?"), len(terms2orfs[term]), "%0.1f" % mr]
                + ["%0.6f" % x for x in [es, pval, qval]]
                + [rvs]
            )
            #self.inputs.outputs.write("\t".join([str(x) for x in vals]) + "\n")
            self.rows.append(vals)
        #self.inputs.outputs.close()

    # ########## Fisher Exact Test ###############

    # HYPERGEOMETRIC
    # scipy.stats.hypergeom.sf() is survival function (1-cdf), so only enriched genes will be significant
    # M = all genes
    # n = category members overall
    # N = sample size (resampling hits)
    # k = number of hits in category (intersection)

    def fisher_exact_test(self):
        import scipy.stats
        from statsmodels.stats import multitest
        
        genes, hits, headers = self.read_resampling_file(
            self.inputs.resampling_file
        )  # use self.Qval_col to determine hits
        associations = self.read_associations(self.inputs.associations_file)
        pathways = self.read_pathways(self.inputs.pathways_file)

        # how many genes are there, and how many have associations?
        # how many genes in associations are not listed in resampling file?
        # do all associations have a definition in pathways?
        # how many pathways have >1 gene? (out of total?) what is max?

        genes_with_associations = 0
        for gene in genes:
            orf = gene[0]
            if orf in associations:
                genes_with_associations += 1
        # self.rows.append("# method=FET, PC=%s" % self.inputs.pseudocount)
        # self.rows.append(
        #     "# genes with associations=%s out of %s total"
        #     % (genes_with_associations, len(genes))
        # )
        # self.rows.append("# significant genes (qval<0.05): %s" % (len(hits)))

        terms = list(pathways.keys())
        terms.sort()
        term_counts = [len(associations.get(term, [])) for term in terms]
        goodterms = []
        for term, cnt in zip(terms, term_counts):
            if cnt > 1:
                goodterms.append(term)
        # self.rows.append(
        #     "# %s out of %s pathways have >=1 gene; max has %s"
        #     % (
        #         len(goodterms),
        #         len(terms),
        #         term_counts[term_counts.index(max(term_counts))],
        #     )
        # )

        results = []
        for term in goodterms:
            n = len(associations[term])  # number of pathway members overall
            M = len(genes)  # total genes
            N = len(hits)  # number of resampling hits
            intersection = list(filter(lambda x: x in associations[term], hits))
            k = len(intersection)
            # add pseudo-counts
            pseudocount = self.inputs.pseudocount
            k_pseudocount = int(k + pseudocount)
            n_pseudocount = n + int(
                M * pseudocount / float(N)
            )  # add same proportion to overall, round it
            expected = round((N * n / float(M)), 2)
            enrichment = round((k + pseudocount) / (expected + pseudocount), 3)
            pval = scipy.stats.hypergeom.sf(k_pseudocount, M, n_pseudocount, N)
            results.append([term, M, n, N, k, expected, k_pseudocount, n_pseudocount, enrichment, pval])

        pvals = [x[-1] for x in results]
        rej, qvals = multitest.fdrcorrection(pvals)
        results = [x + [y] for x, y in zip(results, qvals)]

        genenames = {}
        for gene in genes:
            genenames[gene[0]] = gene[1]

        #header = "#pathway total_genes(M) genes_in_path(n) significant_genes(N) signif_genes_in_path(k) expected k+PC n_adj_by_PC enrichement pval qval description genes"
        #print("\t".join(header.split()))
        #self.rows.append(header.split())

        results.sort(key=lambda x: x[-2])  # pvals
        for res in results:
            vals = res
            term = res[0]
            vals.append(pathways[term])
            intersection = list(filter(lambda x: x in associations[term], hits))
            intersection = ["%s/%s" % (x, genenames[x]) for x in intersection]
            vals.append(" ".join(intersection))
            #print("\t".join([str(x) for x in vals]))
            self.rows.append(vals)

        #results_area.add(self.inputs.output.name)
        #self.finish()
        logging.log("Finished Pathway Enrichment Method")

    # ########## Ontologizer ###############

    # this method is restricted to GO terms

    # this implements the union method of:
    #  Grossman et al (2007). Improved Detection of overrepresentation of Gene-Ontology
    #  annotation with parent-child analysis. Bioinformatics, 23(22):3024-3031.

    def Ontologizer(self):
        import scipy.stats
        from statsmodels.stats import multitest
        
        def warning(s):
            sys.stderr.write("%s\n" % s)  # use self.warning? prepend method name?

        # returns True if ch is a descendant of GO (as GO terms, like "GO:0006810")

        def descendant_of(ch, GO):
            if ch == GO:
                return True
            for node in parents.get(ch, []):
                if descendant_of(node, GO) == True:
                    return True
            return False

        # visited is a hash
        # assume B is above A
        # do upward DFS following parent pointers till hit NULL
        # return all nodes on all paths? include A and B?

        def get_ancestors(A, visited):
            if A in visited:
                return visited
            visited[A] = 1
            for parent in parents.get(A, []):
                get_ancestors(parent, visited)
            return visited

        def depth(go):
            if go not in parents:
                return 0
            return 1 + min([depth(x) for x in parents[go]])

        ontology, parents = {}, {}

        GOannot = self.inputs.associations_file
        OBOfile = self.inputs.pathways_file

        with open(OBOfile) as file:
            for line in file:
                if line[:3] == "id:":
                    id = line[4:-1]
                if line[:5] == "is_a:":
                    parent = line.split()[1]
                    if id not in parents:
                        parents[id] = []
                    parents[id].append(parent)
                if len(line) < 2:
                    id = None
                if line[:5] == "name:":
                    ontology[id] = line[6:-1]

        rv2gos, go2rvs = {}, {}
        MINTERMS, MAXTERMS = 2, 300

        with open(GOannot) as file:
            for line in file:
                w = line.rstrip().split("\t")
                rv, go = w[0], w[1]
                if rv not in rv2gos:
                    rv2gos[rv] = []
                if go not in rv2gos[rv]:
                    rv2gos[rv].append(go)
                if go not in go2rvs:
                    go2rvs[go] = []
                if rv not in go2rvs[go]:
                    go2rvs[go].append(rv)
                # expand to all parents...
                BP, CC, MF = ["GO:0008150", "GO:0005575", "GO:0003674"]
                for g in get_ancestors(go, {}):
                    if g not in rv2gos[rv]:
                        rv2gos[rv].append(g)
                        if g not in go2rvs:
                            go2rvs[g] = []
                        go2rvs[g].append(rv)

        warning(
            "GO terms with at least one ORF: %s" % len(go2rvs.keys())
        )  # what about between MIN and MAX?
        for go in go2rvs.keys():
            if go not in ontology:
                warning("not found: %s" % go)  # also indicate which gene?

        # could use class method, but would have to adapt it:
        # genes,hits,headers = self.read_resampling_file(self.inputs.resampling_file)

        genes, pvals = {}, []
        allorfs, studyset = [], []
        with open(self.inputs.resampling_file) as file:
            for line in file:
                if line[0] == "#":
                    continue
                w = line.rstrip().split("\t")
                genes[w[0]] = w
                allorfs.append(w[0])
                # pval,qval = float(w[-2]),float(w[-1])
                pval, qval = float(w[self.inputs.pval_col]), float(w[self.inputs.qval_col])
                if qval < 0.05:
                    studyset.append(w[0])
            pvals.append((w[0], pval))
        pvals.sort(key=lambda x: x[1])
        ranks = {}
        for i, (rv, pval) in enumerate(pvals):
            ranks[rv] = i + 1
        #self.rows.append("# number of resampling hits (qval<0.05): %s" % len(studyset))

        counts = []
        n, a = len(allorfs), len(studyset)
        for go in go2rvs.keys():
            if go not in ontology:
                continue
            m = len(go2rvs[go])  # orfs in the GO terms in the genome overall
            if m >= MINTERMS and m <= MAXTERMS:
                b = len(
                    list(filter(lambda x: x in go2rvs[go], studyset))
                )  # orfs with GO term in studyset
                if b == 0:
                    continue
                # calc p=len(rvs for parents of go), q=len(subset of parent rvs in studyset)
                P = set()
                if go not in parents:
                    continue
                npar = len(parents[go])
                for par in parents[go]:
                    P = P.union(go2rvs[par])
                p, q = len(P), len(P.intersection(studyset))
                enrich = ((b + 1) / float(q + 1)) / (
                    (m + 1) / float(p + 1)
                )  # enrichment
                # parents: p overall, q in studyset; GO in studyset: b out of q (subset of studyset labeled with a parent GO term)
                # assume m orfs have GO term inside parent (same as overall)
                # subtract 1 from b possibly because we are using 1-cdf? this is needed to get same numeric results as 'phypergeometric' in Ontologizer
                pval = 1.0 - scipy.stats.hypergeom.cdf(b - 1, p, m, q)
                counts.append([go, n, m, p, q, a, b, npar, enrich, pval])

        pvals = [x[-1] for x in counts]
        qvals = multitest.fdrcorrection(pvals)[1]

        counts = [x + [y] for x, y in zip(counts, qvals)]
        counts.sort(key=lambda x: x[-1])

        # self.rows.append(
        #     "# number of GO terms that are significantly enriched (qval<0.05): %s"
        #     % len(list(filter(lambda x: x[-1] < 0.05, counts)))
        # )

        # 1-sided pvals, only report enriched terms, not depleted

        # self.rows.append(
        #     "\t".join(
        #         "GO_term description total_orfs orfs_in_GO orfs_in_parents hits_in_parents hits GO_in_hits num_parent_nodes enrichment pval qval genes_in_intersection_of_hits_and_GO_term".split()
        #     )
        # )  # assume self.inputs.outputFile has already been opened
        for (go, n, m, p, q, a, b, npar, enrich, pval, qval) in counts:
            hits = filter(lambda x: x in go2rvs[go], studyset)
            hits = [(x, genes[x][1], ranks[x]) for x in hits]
            hits.sort(key=lambda x: x[2])  # sort on ranks
            hits = ["%s/%s(%s)" % (rv, gene, rank) for (rv, gene, rank) in hits]
            hits = ",".join(hits)  # add gene names
            vals = [
                go,
                ontology[go],
                n,
                m,
                a,
                b,
                p,
                q,
                npar,
                round(enrich, 3),
                round(pval, 6),
                round(qval, 6),
                hits,
            ]
            #print("\t".join([str(x) for x in vals]))
            self.rows.append(vals)

        logging.log("Adding File: %s" % (self.inputs.output_path))
        #results_area.add(self.output.name)
        #self.finish()
        logging.log("Finished Pathway Enrichment Method")
    


@transit_tools.ResultsFile
class ResultFileType1:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Analysis.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Analysis.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Analysis.identifier,
                    heading="",
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[
                        "qval"
                    ],
                ).Show(),
            })
        )
        
        # 
        # get column names
        # 
        comments, headers, rows = csv.read(self.path, seperator="\t", skip_empty_lines=True, comment_symbol="#")
        if len(comments) == 0:
            raise Exception(f'''No comments in file, and I expected the last comment to be the column names, while to load Anova file "{self.path}"''')
        self.column_names = comments[-1].split("\t")
        
        # 
        # get rows
        #
        self.rows = []
        for each_row in rows:
            self.rows.append({
                each_column_name: each_cell
                    for each_column_name, each_cell in zip(self.column_names, each_row)
            })
        
        # 
        # get summary stats
        #
        self.values_for_result_table.update({
            # HANDLE_THIS (additional summary_info for results table)
            # examples:
                # f"Gene Count": len(self.rows),
                # f"Padj<{Analysis.significance_threshold}": len([
                #     1 for each in self.rows
                #         if each.get("Padj", 0) < Analysis.significance_threshold 
                # ]),
        })
    
    def __str__(self):
        return f"""
            File for {Analysis.short_name}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()

Analysis.filetypes = [ ResultFileType1, ] # HANDLE_THIS
Method = GUI = Analysis # for compatibility with older code/methods