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
from pytransit.basics.lazy_dict import LazyDict

import pytransit.tools.gui_tools as gui_tools
import pytransit.components.file_display as file_display
import pytransit.tools.transit_tools as transit_tools
import pytransit.tools.tnseq_tools as tnseq_tools
import pytransit.tools.norm_tools as norm_tools
import pytransit.tools.stat_tools as stat_tools
import pytransit.tools.console_tools as console_tools
import pytransit.basics.csv as csv
import pytransit.components.results_area as results_area
from pytransit.tools.transit_tools import wx, pub, basename, HAS_R, FloatVector, DataFrame, StrVector, EOL
from pytransit.universal_data import universal
from pytransit.components.parameter_panel import panel as parameter_panel
from pytransit.components.parameter_panel import panel, progress_update
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.panel_helpers import make_panel, create_run_button, create_normalization_input, create_reference_condition_input, create_include_condition_list_input, create_exclude_condition_list_input, create_n_terminus_input, create_c_terminus_input, create_pseudocount_input, create_winsorize_input, create_alpha_input, create_button
command_name = sys.argv[0]

class Analysis:
    identifier  = "Anova"
    short_name  = "anova"
    long_name   = "ANOVA"
    short_desc  = "Perform Anova analysis"
    long_desc   = """Perform Anova analysis"""
    transposons = [ "himar1", "tn5" ]
    
    significance_threshold = 0.05
    
    inputs = LazyDict(
        combined_wig=None,
        metadata=None,
        annotation=None,
        normalization=None,
        output_path=None,
        
        excluded_conditions=[],
        included_conditions=[],
        n_terminus=0.0,
        c_terminus=0.0,
        pseudocount=5,
        winz=False,
        refs=[],
        alpha=1000,
    )
    
    valid_cli_flags = [
        "-n",
        "--include-conditions",
        "--exclude-conditions",
        "--ref",
        "-iN",
        "-iC",
        "-PC",
        "-alpha",
        "-winz",
    ]
    usage_string = f"""
        Usage: python3 transit.py anova <combined wig file> <samples_metadata file> <annotation .prot_table> <output file> [Optional Arguments]
        Optional Arguments:
            -n <string>         :=  Normalization method. Default: -n TTR
            --include-conditions <cond1,...> := Comma-separated list of conditions to use for analysis (Default: all)
            --exclude-conditions <cond1,...> := Comma-separated list of conditions to exclude (Default: none)
            --ref <cond> := which condition(s) to use as a reference for calculating LFCs (comma-separated if multiple conditions)
            -iN <N> :=  Ignore TAs within given percentage (e.g. 5) of N terminus. Default: -iN 0
            -iC <N> :=  Ignore TAs within given percentage (e.g. 5) of C terminus. Default: -iC 0
            -PC <N> := pseudocounts to use for calculating LFCs. Default: -PC 5
            -alpha <N> := value added to mse in F-test for moderated anova (makes genes with low counts less significant). Default: -alpha 1000
            -winz   := winsorize insertion counts for each gene in each condition (replace max cnt with 2nd highest; helps mitigate effect of outliers)
    """.replace("\n        ", "\n")
    
    
    wxobj = None
    panel = None
    
    def __init__(self, *args, **kwargs):
        Analysis.instance = self
        self.full_name        = f"[{self.short_name}]  -  {self.short_desc}"
        self.transposons_text = transit_tools.get_transposons_text(self.transposons)
        self.filetypes        = [File]
        self.method           = Analysis # backwards compat
        self.gui              = self     # backwards compat
    
    def __str__(self):
        return f"""
            Analysis Method:
                Short Name:  {self.short_name}
                Long Name:   {self.long_name}
                Short Desc:  {self.short_desc}
                Long Desc:   {self.long_desc}
                GUI:         {self.gui}
        """.replace('\n            ','\n').strip()
    
    def __repr__(self):
        return f"{self.inputs}"

    def define_panel(self, _):
        self.panel = make_panel()

        # 
        # parameter inputs
        # 
        # --include-conditions <cond1,...> := Comma-separated list of conditions to use for analysis (Default: all)
        # --exclude-conditions <cond1,...> := Comma-separated list of conditions to exclude (Default: none)
        # --ref <cond> := which condition(s) to use as a reference for calculating lfc_s (comma-separated if multiple conditions)
        # -iN <N> :=  Ignore TAs within given percentage (e.g. 5) of N terminus. Default: -iN 0
        # -iC <N> :=  Ignore TAs within given percentage (e.g. 5) of C terminus. Default: -iC 0
        # -PC <N> := pseudocounts to use for calculating LFC. Default: -PC 5
        # -winz   := winsorize insertion counts for each gene in each condition (replace max cnt with 2nd highest; helps mitigate effect of outliers)
        self.value_getters = LazyDict()
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        if True:
            self.value_getters.included_conditions    = create_include_condition_list_input(self.panel, main_sizer)
            self.value_getters.excluded_conditions    = create_exclude_condition_list_input(self.panel, main_sizer)
            self.value_getters.reference_condition    = create_reference_condition_input(self.panel, main_sizer)
            self.value_getters.n_terminus             = create_n_terminus_input(self.panel, main_sizer)
            self.value_getters.c_terminus             = create_c_terminus_input(self.panel, main_sizer)
            self.value_getters.normalization          = create_normalization_input(self.panel, main_sizer)
            self.value_getters.pseudocount            = create_pseudocount_input(self.panel, main_sizer)
            self.value_getters.alpha                  = create_alpha_input(self.panel, main_sizer)
            self.value_getters.winz                   = create_winsorize_input(self.panel, main_sizer)
            self.value_getters.refs                   = lambda *args: [] if self.value_getters.reference_condition() == "[None]" else [ self.value_getters.reference_condition() ]
            
            create_run_button(self.panel, main_sizer)
            
        parameter_panel.set_panel(self.panel)
        self.panel.SetSizer(main_sizer)
        self.panel.Layout()
        main_sizer.Fit(self.panel)

    @classmethod
    def from_gui(cls, frame):
        with gui_tools.nice_error_log:
            # 
            # get wig files
            # 
            combined_wig = universal.session_data.combined_wigs[0]
            Analysis.inputs.combined_wig = combined_wig.main_path
            Analysis.inputs.metadata     = combined_wig.metadata.path
            
            # 
            # get annotation
            # 
            Analysis.inputs.annotation_path = universal.session_data.annotation_path
            # FIXME: enable this once I get a valid annotation file example
            # if not transit_tools.validate_annotation(Analysis.inputs.annotation):
            #     return None
            
            # 
            # setup custom inputs
            # 
            for each_key, each_getter in Analysis.instance.value_getters.items():
                try:
                    Analysis.inputs[each_key] = each_getter()
                except Exception as error:
                    raise Exception(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
            transit_tools.log("included_conditions", Analysis.inputs.included_conditions)
            # 
            # save result files
            # 
            Analysis.inputs.output_path = gui_tools.ask_for_output_file_path(
                default_file_name="test1_output.dat",
                output_extensions='Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*',
            )
            if not Analysis.inputs.output_path:
                return None

            return Analysis.instance

    @classmethod
    def from_args(cls, args, kwargs):
        console_tools.handle_help_flag(kwargs, cls.usage_string)
        console_tools.handle_unrecognized_flags(cls.valid_cli_flags, kwargs, cls.usage_string)

        combined_wig      = args[0]
        annotation_path   = args[2]
        metadata          = args[1]
        output_path       = args[3]
        normalization     = kwargs.get("n", "TTR")
        n_terminus        = float(kwargs.get("iN", 0.0))
        c_terminus        = float(kwargs.get("iC", 0.0))
        winz              = "winz" in kwargs
        pseudocount       = int(kwargs.get("PC", 5))
        alpha             = float(kwargs.get("alpha", 1000))
        refs              = kwargs.get("-ref", [])  # list of condition names to use a reference for calculating lfc_s
        if refs != []: refs = refs.split(",")
        excluded_conditions = list( filter(None, kwargs.get("-exclude-conditions", "").split(",")) )
        included_conditions = list( filter(None, kwargs.get("-include-conditions", "").split(",")) )

        # save all the data
        Analysis.inputs.update(dict(
            combined_wig=combined_wig,
            metadata=metadata,
            annotation_path=annotation_path,
            normalization=normalization,
            output_path=output_path,
            
            excluded_conditions=excluded_conditions,
            included_conditions=included_conditions,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            pseudocount=pseudocount,
            winz=winz,
            refs=refs,
            alpha=alpha,
        ))
        
        return Analysis.instance
        
    def means_by_condition_for_gene(self, sites, conditions, data):
        """
            Returns a dictionary of {Condition: Mean} for each condition.
            ([Site], [Condition]) -> {Condition: Number}
            Site :: Number
            Condition :: String
        """
        n_ta_sites = len(sites)
        wigs_by_conditions = collections.defaultdict(lambda: [])
        for i, c in enumerate(conditions):
            wigs_by_conditions[c].append(i)

        return {
            c: numpy.mean(transit_tools.winsorize(data[wigIndex][:, sites]))
            if n_ta_sites > 0
            else 0
            for (c, wigIndex) in wigs_by_conditions.items()
        }

    def means_by_rv(self, data, rv_site_indexes_map, genes, conditions):
        """
            Returns Dictionary of mean values by condition
            ([[Wigdata]], {rv: SiteIndex}, [Gene], [Condition]) -> {rv: {Condition: Number}}
            Wigdata :: [Number]
            SiteIndex :: Number
            Gene :: {start, end, rv, gene, strand}
            Condition :: String
        """
        means_by_rv = {}
        for gene in genes:
            rv = gene["rv"]
            means_by_rv[rv] = self.means_by_condition_for_gene(
                rv_site_indexes_map[rv], conditions, data
            )
        return means_by_rv

    def group_by_condition(self, wig_list, conditions):
        """
            Returns array of datasets, where each dataset corresponds to one condition.
            ([[Wigdata]], [Condition]) -> [[DataForCondition]]
            Wigdata :: [Number]
            Condition :: String
            DataForCondition :: [Number]
        """
        counts_by_condition = collections.defaultdict(lambda: [])
        count_sum = 0
        for i, c in enumerate(conditions):
            count_sum += numpy.sum(wig_list[i])
            counts_by_condition[c].append(wig_list[i])

        return (
            count_sum,
            [numpy.array(v).flatten() for v in counts_by_condition.values()],
        )

    def calculate_anova(self, data, genes, means_by_rv, rv_site_indexes_map, conditions):
        """
            Runs Anova (grouping data by condition) and returns p and q values
            ([[Wigdata]], [Gene], {rv: {Condition: Mean}}, {rv: [SiteIndex]}, [Condition]) -> Tuple([Number], [Number])
            Wigdata :: [Number]
            Gene :: {start, end, rv, gene, strand}
            Mean :: Number
            SiteIndex: Integer
            Condition :: String
        """
        import scipy
        import scipy.stats
        import statsmodels.stats.multitest
        
        count = 0

        msrs, mses, f_stats, pvals, rvs, status = [],[],[],[],[],[]
        for gene in genes:
            count += 1
            rv = gene["rv"]
            if len(rv_site_indexes_map[rv]) <= 1:
                status.append("TA sites <= 1")
                msr,mse,f_stat,pval = 0,0,-1,1
            else:
                count_sum, counts_vec = self.group_by_condition(
                    list(map(lambda wigData: wigData[rv_site_indexes_map[rv]], data)),
                    conditions,
                )
                if self.inputs.winz:
                    counts_vec = transit_tools.winsorize(counts_vec)

                if count_sum == 0:
                    msr,mse,f_stat,pval = 0,0,-1,1
                    status.append("No counts in all conditions")
                else:
                    f_stat,pval = scipy.stats.f_oneway(*counts_vec)
                    status.append("-")
                    # counts_vec is a list of numpy arrays, or could be a list of lists
                    # pooled counts for each condition, over TAs in gene and replicates
                    if isinstance(counts_vec[0],numpy.ndarray): 
                      counts_vec_as_arrays = counts_vec
                      counts_vecAsLists = [grp.tolist() for grp in counts_vec]
                    else:
                      counts_vec_as_arrays = [numpy.array(grp) for grp in counts_vec]
                      counts_vecAsLists = counts_vec
                    all_counts = [item for sublist in counts_vecAsLists for item in sublist]
                    grand_mean = numpy.mean(all_counts)
                    group_means = [numpy.mean(grp) for grp in counts_vec_as_arrays]
                    k,n = len(counts_vec),len(all_counts)
                    df_between,df_within = k-1,n-k
                    msr,mse = 0,0
                    for grp in counts_vec_as_arrays: msr += grp.size*(numpy.mean(grp)-grand_mean)**2/float(df_between)
                    for grp,mn in zip(counts_vec_as_arrays,group_means): mse += numpy.sum((grp-mn)**2) 
                    mse /= float(df_within)
                    mse = mse+self.inputs.alpha ### moderation
                    f_mod = msr/float(mse)
                    Pmod = scipy.stats.f.sf(f_mod, df_between, df_within)
                    f_stat,pval = f_mod,Pmod
            pvals.append(pval)   
            f_stats.append(f_stat) 
            msrs.append(msr)
            mses.append(mse)
            rvs.append(rv)

            # Update progress
            percentage = 100.0 * count / len(genes)
            progress_update(f"Running Anova Method... {percentage:5.1f}%", percentage)

        pvals = numpy.array(pvals)
        mask = numpy.isfinite(pvals)
        qvals = numpy.full(pvals.shape, numpy.nan)
        qvals[mask] = statsmodels.stats.multitest.fdrcorrection(pvals[mask])[1]  # BH, alpha=0.05

        msr, mse, f, p, q, status_map = {},{},{},{},{},{}
        for i,rv in enumerate(rvs):
            msr[rv], mse[rv], f[rv], p[rv], q[rv], status_map[rv] = msrs[i], mses[i], f_stats[i], pvals[i], qvals[i], status[i]
        return msr, mse, f, p, q, status_map
    
    def calc_lfcs(self, means, refs=[], pseudocount=5):
        if len(refs) == 0:
            refs = means  # if ref condition(s) not explicitly defined, use mean of all
        grand_mean = numpy.mean(refs)
        lfcs = [math.log((x + pseudocount) / float(grand_mean + pseudocount), 2) for x in means]
        return lfcs

    def Run(self):
        with gui_tools.nice_error_log:
            transit_tools.log("Starting Anova analysis")
            start_time = time.time()
            
            # 
            # get data
            # 
            transit_tools.log("Getting Data")
            if True:
                sites, data, filenames_in_comb_wig = tnseq_tools.read_combined_wig(self.inputs.combined_wig)
                
                transit_tools.log(f"Normalizing using: {self.inputs.normalization}")
                data, factors = norm_tools.normalize_data(data, self.inputs.normalization)
                
                if self.inputs.winz: transit_tools.log("Winsorizing insertion counts")
                conditions_by_file, _, _, ordering_metadata = tnseq_tools.read_samples_metadata(self.inputs.metadata)
                conditions = [ conditions_by_file.get(f, None) for f in filenames_in_comb_wig ]
                conditions_list = transit_tools.select_conditions(
                    conditions=conditions,
                    included_conditions=self.inputs.included_conditions,
                    excluded_conditions=self.inputs.excluded_conditions,
                    ordering_metadata=ordering_metadata,
                )

                condition_names = [conditions_by_file[f] for f in filenames_in_comb_wig]

                (
                    data,
                    file_names,
                    condition_names,
                    conditions,
                    _,
                    _,
                ) = transit_tools.filter_wigs_by_conditions3(
                    data,
                    file_names=filenames_in_comb_wig, # it looks like file_names and condition_names have to be parallel to data (vector of wigs)
                    condition_names=condition_names, # original Condition column in samples metadata file
                    included_cond=self.inputs.included_conditions,
                    excluded_cond=self.inputs.excluded_conditions,
                    conditions=condition_names,
                ) # this is kind of redundant for ANOVA, but it is here because condition, covars, and interactions could have been manipulated for ZINB
                
                transit_tools.log("reading genes")
                genes = tnseq_tools.read_genes(self.inputs.annotation_path)
            
            # 
            # process data
            # 
            if True:
                transit_tools.log("processing data")
                TASiteindexMap = {ta: i for i, ta in enumerate(sites)}
                rv_site_indexes_map = tnseq_tools.rv_siteindexes_map(
                    genes, TASiteindexMap, n_terminus=self.inputs.n_terminus, c_terminus=self.inputs.c_terminus
                )
                means_by_rv = self.means_by_rv(data, rv_site_indexes_map, genes, conditions)

                transit_tools.log("Running Anova")
                msrs, mses, f_stats, pvals, qvals, run_status = self.calculate_anova(
                    data, genes, means_by_rv, rv_site_indexes_map, conditions
                )
            
            # 
            # write output
            # 
            if True:
                transit_tools.log(f"Adding File: {self.inputs.output_path}")
                
                # 
                # generate rows
                # 
                rows = []
                for gene in genes:
                    each_rv = gene["rv"]
                    if each_rv in means_by_rv:
                        means = [ means_by_rv[each_rv][condition_name] for condition_name in conditions_list]
                        refs  = [ means_by_rv[each_rv][ref_condition ] for ref_condition in self.inputs.refs]
                        lfcs = self.calc_lfcs(means, refs, self.inputs.pseudocount)
                        rows.append(
                            [
                                each_rv,
                                gene["gene"],
                                str(len(rv_site_indexes_map[each_rv])),
                            ] + [
                                "%0.2f" % x for x in means
                            ] + [
                                "%0.3f" % x for x in lfcs
                            ] +  [
                                "%f" % x for x in [msrs[each_rv], mses[each_rv], f_stats[each_rv], pvals[each_rv], qvals[each_rv]]
                            ] + [
                                run_status[each_rv]
                            ]
                        )
                
                # 
                # write to file
                # 
                transit_tools.write_result(
                    path=self.inputs.output_path,
                    file_kind=Analysis.identifier,
                    rows=rows,
                    column_names=[
                        "Rv",
                        "Gene",
                        "TAs",
                        *[ f"Mean_{condition_name}" for condition_name in conditions_list ],
                        *[  f"LFC_{condition_name}" for condition_name in conditions_list ],
                        "MSR",
                        "MSE+alpha",
                        "Fstat",
                        "Pval",
                        "Padj",
                        "status"
                    ],
                    extra_info=dict(
                        parameters=dict(
                            normalization=self.inputs.normalization,
                            trimming=f"{self.inputs.n_terminus}/{self.inputs.c_terminus} % (N/C)",
                            pseudocounts=self.inputs.pseudocount,
                            alpha=self.inputs.alpha,
                        ),
                    ),
                )
                transit_tools.log("Finished Anova analysis")
                transit_tools.log(f"Time: {time.time() - start_time:0.1f}s\n")
            results_area.add(self.inputs.output_path)

@transit_tools.ResultsFile
class File(Analysis):
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
                "Display Table": lambda *args: SpreadSheet(title="Anova",heading="",column_names=self.column_names,rows=self.rows, sort_by=["Padj", "Pval"]).Show(),
                "Display Heatmap": lambda *args: self.create_heatmap(infile=self.path, output_path=self.path+".heatmap.png"),
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
            f"Gene Count": len(self.rows),
            f"Padj<{Analysis.significance_threshold}": len([
                1 for each in self.rows
                    if each.get("Padj", 0) < Analysis.significance_threshold 
            ]),
        })
    
    def __str__(self):
        return f"""
            File for {self.short_name}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()
    
    def display_histogram(self, display_frame, event):
        pass
        # gene = display_frame.grid.GetCellValue(display_frame.row, 0)
        # filepath = os.path.join(
        #     ntpath.dirname(display_frame.path),
        #     transit_tools.fetch_name(display_frame.path),
        # )
        # filename = os.path.join(filepath, gene + ".png")
        # if os.path.exists(filename):
        #     imgWindow = pytransit.components.file_display.ImgFrame(None, filename)
        #     imgWindow.Show()
        # else:
        #     transit_tools.show_error_dialog("Error Displaying File. Histogram image not found. Make sure results were obtained with the histogram option turned on.")
        #     print("Error Displaying File. Histogram image does not exist.")

    def create_heatmap(self, infile, output_path, topk=-1, qval=0.05, low_mean_filter=5):
        with gui_tools.nice_error_log:
            if not HAS_R:
                raise Exception(f'''Error: R and rpy2 (~= 3.0) required to run Heatmap''')
            headers = None
            data, hits = [], []
            number_of_conditions = -1

            with open(infile) as file:
                for line in file:
                    w = line.rstrip().split("\t")
                    if line[0] == "#" or (
                        "pval" in line and "padj" in line
                    ):  # check for 'pval' for backwards compatibility
                        headers = w
                        continue  # keep last comment line as headers
                    # assume first non-comment line is header
                    if number_of_conditions == -1:
                        # ANOVA header line has names of conditions, organized as 3+2*number_of_conditions+3 (2 groups (means, lfc_s) X number_of_conditions conditions)
                        number_of_conditions = int((len(w) - 6) / 2)
                        headers = headers[3 : 3 + number_of_conditions]
                        headers = [x.replace("Mean_", "") for x in headers]
                    else:
                        means = [
                            float(x) for x in w[3 : 3 + number_of_conditions]
                        ]  # take just the columns of means
                        lfcs = [
                            float(x) for x in w[3 + number_of_conditions : 3 + number_of_conditions + number_of_conditions]
                        ]  # take just the columns of lfc_s
                        each_qval = float(w[-2])
                        data.append((w, means, lfcs, each_qval))
            
            data.sort(key=lambda x: x[-1])
            hits, lfc_s = [], []
            for k, (w, means, lfcs, each_qval) in enumerate(data):
                if (topk == -1 and each_qval < qval) or (
                    topk != -1 and k < topk
                ):
                    mm = round(numpy.mean(means), 1)
                    if mm < low_mean_filter:
                        print("excluding %s/%s, mean(means)=%s" % (w[0], w[1], mm))
                    else:
                        hits.append(w)
                        lfc_s.append(lfcs)

            print("heatmap based on %s genes" % len(hits))
            gene_names = ["%s/%s" % (w[0], w[1]) for w in hits]
            hash = {}
            headers = [h.replace("Mean_", "") for h in headers]
            for i, col in enumerate(headers):
                hash[col] = FloatVector([x[i] for x in lfc_s])
            df = DataFrame(hash)
            transit_tools.r_heatmap_func(df, StrVector(gene_names), output_path)
            
            # add it as a result
            results_area.add(output_path)
            gui_tools.show_image(output_path)
    
Method = GUI = Analysis
Analysis() # make sure there's one instance
