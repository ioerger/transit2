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

from pytransit.generic_tools import csv, misc, informative_iterator
from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components import samples_area, results_area, parameter_panel, file_display

from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.specific_tools.transit_tools import wx, basename, HAS_R, FloatVector, DataFrame, StrVector, EOL, SEPARATOR, rpackages
from pytransit.specific_tools.tnseq_tools import ProtTable
from pytransit.components.spreadsheet import SpreadSheet

cli_args = LazyDict(
    gene=None, # for debugging
)

@misc.singleton
class Method:
    name = "Zinb"
    identifier  = name
    cli_name    = name.lower()
    menu_name   = f"{name} - Perform {name} analysis"
    description   = f"""Perform {name} analysis"""
    
    valid_cli_flags = [
        "-n",
        "--exclude-conditions",
        "--include-conditions",
        "--ref",
        "-iN",
        "-iC",
        "-winz",
        "-PC",
        "--group-by",
        "--condition",
        "--covars",
        "--interactions",
        "--prot_table_path",
        "--gene",
    ]
    
    usage_string = f"""{console_tools.subcommand_prefix} {cli_name} <combined wig file> <samples_metadata file> <annotation .prot_table_path> <output file> [Optional Arguments]
        Optional Arguments:
            -n <string>         :=  Normalization method. Default: -n TTR
            --exclude-conditions <cond1,cond2> :=  Comma separated list of conditions to exclude, for the analysis.
            --include-conditions <cond1,cond2> :=  Comma separated list of conditions to include, for the analysis. Conditions not in this list, will be excluded.
            --ref <cond>    := which condition(s) to use as a reference for calculating log_fold_changes (comma-separated if multiple conditions)
            -iN <float>     := Ignore TAs occuring within given percentage (as integer) of the N terminus. Default: -iN 5
            -iC <float>     := Ignore TAs occuring within given percentage (as integer) of the C terminus. Default: -iC 5
            -winz           := winsorize insertion counts for each gene in each condition (replace max cnt with 2nd highest; helps mitigate effect of outliers)
            -PC <N>         := pseudocounts to use for calculating log_fold_changes. Default: -PC 5
            --group-by      := columnname (in samples_metadata) to use as the Condition. Default: "Condition"
            --condition     := alias for --group-by
            --covars <covar1,covar2...>       := Comma separated list of covariates (in metadata file) to include, for the analysis.
            --interactions <covar1,covar2...> := Comma separated list of covariates to include, that interact with the condition for the analysis. Must be factors
            --prot_table_path <filename>           := for appending annotations of genes
            --gene <RV number or Gene name>   := Run method for one gene and print model output.
    """.replace("\n        ", "\n")

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=4)
        
        cli_args.gene = kwargs["-gene"]
        
        # save the data
        Method.output(
            combined_wig_path = args[0],
            metadata_path = args[1],
            annotation_path = args[2],
            output_path = args[3],
            prot_table_path=kwargs["prot_table"],
            group_by=kwargs.get("-group-by", kwargs["-condition"]),
            covars=console_tools.string_arg_to_list(kwargs["-covars"]),
            interactions=console_tools.string_arg_to_list(kwargs["-interactions"]),
            refs=console_tools.string_arg_to_list(kwargs["-ref"]),
            excluded_conditions=console_tools.string_arg_to_list(kwargs["-exclude-conditions"]),
            included_conditions=console_tools.string_arg_to_list(kwargs["-include-conditions"]),
            winz="winz" in kwargs,
            normalization=kwargs["n"],
            n_terminus=kwargs["iN"],
            c_terminus=kwargs["iC"],
            pseudocount=kwargs["PC"],
            disable_logging=False,
        )

    
    @gui.add_menu("Method", "himar1", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    @gui.add_menu("Method", "tn5", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    def define_panel(self, _):
        from pytransit.components import panel_helpers
        with panel_helpers.NewPanel() as (panel, main_sizer):
            parameter_panel.set_instructions(
                method_short_text=self.name,
                method_long_text="",
                method_specific_instructions="""
                    The ZINB (Zero-Inflated Negative Binomial) method is used to determine which genes exhibit statistically significant variability across multiple conditions, in either the magnitude of insertion counts or local saturation, agnostically (in any one condition compared to the others). Like ANOVA, the ZINB method takes a combined_wig file (which combines multiple datasets in one file) and a samples_metadata file (which describes which samples/replicates belong to which experimental conditions).
                    
                    ZINB can be applied to two or more conditions at a time. Thus it subsumes resampling. Our testing suggests that ZINB typically identifies 10-20% more varying genes than resampling (and vastly out-performs ANOVA for detecting significant variability across conditions). Furthermore, because of how ZINB treats magnitude of read counts separately from local saturation in a gene, it occasionally identifies genes with variability not detectable by resampling analysis.
                    
                    1. Add an annotation file for the organism corresponding to the desired datasets
                    
                    2. FIXME
                    
                    3. FIXME
                    
                    4. FIXME
                """.replace("\n                    ","\n"),
            )
            self.value_getters = LazyDict()
            try:
                metadata_headers = gui.combined_wigs[-1].metadata.headers
            except Exception as error:
                pass
            
            self.value_getters = LazyDict(
                included_conditions= panel_helpers.create_selected_condition_names_input(panel, main_sizer),
                excluded_conditions= (lambda *args: []), # never needed, but exists to comply with CLI interface
                group_by=            panel_helpers.create_choice_input(panel, main_sizer, label="Group By", options=["Condition", *metadata_headers], default_option="Condition", tooltip_text="Column name (in samples_metadata) to use as the primary Condition being evaluated (to test for significant variability of insertions among conditions)."),
                # FIXME: covars needs to be an empty list by default
                covars=              panel_helpers.create_choice_input(panel, main_sizer, label="Covars", options=["Condition", *metadata_headers], default_option="Condition", tooltip_text="Comma separated list of covariates (columns in metadata file) to include, for the analysis. If additional covariates distinguishing the samples are available, such as library, timepoint, or genotype, they may be “factored out” of the test of the primary condition. (variation due to covars is accounted for in the model, but not considered in evaluating the effect on variability due to the primary condition)"), 
                # FIXME: interactions needs to be an empty list by default
                interactions=        panel_helpers.create_choice_input(panel, main_sizer, label="Interactions", options=["Condition", *metadata_headers], default_option="Condition", tooltip_text="Comma separated list of covariates (cols in metadata) to include, that interact with the condition for the analysis. (variation due to these variables is included in testing the effect of the main condition)"), 
                refs=                panel_helpers.create_reference_condition_input(panel, main_sizer), # FIXME: currently returns a string, should return a list
                n_terminus=          panel_helpers.create_n_terminus_input(panel, main_sizer),
                c_terminus=          panel_helpers.create_c_terminus_input(panel, main_sizer),
                normalization=       panel_helpers.create_normalization_input(panel, main_sizer),
                pseudocount=         panel_helpers.create_pseudocount_input(panel, main_sizer),
                alpha=               panel_helpers.create_alpha_input(panel, main_sizer),
                winz=                panel_helpers.create_winsorize_input(panel, main_sizer),
                prot_table_path=     panel_helpers.create_file_input(panel, main_sizer, button_label="Add ProtTable (optional)", tooltip_text="FIXME", popup_title="ProtTable"),
            )
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function=self.from_gui)
            
    @staticmethod
    def from_gui(frame):
        arguments = LazyDict()
        
        # 
        # global data
        # 
        arguments.combined_wig_path = gui.combined_wigs[-1].main_path
        arguments.metadata_path     = gui.combined_wigs[-1].metadata_path
        arguments.annotation_path   = gui.combined_wigs[-1].annotation_path
        
        # 
        # call all GUI getters, puts results into respective arguments key-value
        # 
        for each_key, each_getter in Method.value_getters.items():
            try:
                arguments[each_key] = each_getter()
            except Exception as error:
                logging.error(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')
        
        # 
        # ask for output path(s)
        # 
        arguments.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Method.cli_name}_output.tsv",
            output_extensions=transit_tools.result_output_extensions,
        )
        # if user didn't select an output path
        if not arguments.output_path:
            return None

        Method.output(**arguments)

    @staticmethod
    def output(*,
        combined_wig_path,
        metadata_path,
        annotation_path,
        output_path,
        prot_table_path=None,
        group_by=None,
        covars=None,
        interactions=None,
        refs=None,
        excluded_conditions=None,
        included_conditions=None,
        winz=None,
        pseudocount=None,
        normalization=None,
        n_terminus=None,
        c_terminus=None,
        disable_logging=False,
    ): # output()
        # Defaults (even if argument directly provided as None)
        prot_table_path     = prot_table_path     if prot_table_path     is not None else None
        group_by            = group_by            if group_by            is not None else "Condition"
        covars              = covars              if covars              is not None else []
        interactions        = interactions        if interactions        is not None else []
        refs                = refs                if refs                is not None else []
        excluded_conditions = excluded_conditions if excluded_conditions is not None else []
        included_conditions = included_conditions if included_conditions is not None else []
        winz                = winz                if winz                is not None else False
        pseudocount         = pseudocount         if pseudocount         is not None else 5.0 # TODO: check later to make sure this is the correct default --Jeff
        normalization       = normalization       if normalization       is not None else "TTR"
        n_terminus          = n_terminus          if n_terminus          is not None else 5.0
        c_terminus          = c_terminus          if c_terminus          is not None else 5.0
        
        # transit_tools.require_r_to_be_installed(required_r_packages=[ "MASS", "pscl" ]) # FIXME: uncomment this
        with transit_tools.TimerAndOutputs(method_name=Method.identifier, output_paths=[output_path], disable=disable_logging) as timer:
            # 
            # process data
            # 
            if True:
                # 
                # create gene_name_to_description from prot_table
                # 
                gene_name_to_description = None
                if prot_table_path != None:
                    gene_name_to_description = {}
                    with open(prot_table_path) as file:
                        for line in file:
                            row = line.rstrip().split("\t")
                            gene_name        = row[ProtTable.gene_name_index]
                            gene_description = row[ProtTable.gene_name_index]
                            gene_name_to_description[gene_name] = gene_description
                
                # 
                # gather comwig data
                # 
                logging.log("Getting Data")
                (sites, data, filenames_in_comb_wig) = tnseq_tools.CombinedWigData.load(
                    combined_wig_path
                )
                
                # 
                # normalize if needed
                # 
                logging.log("Normalizing using: %s" % normalization)
                (data, factors) = norm_tools.normalize_data(data, normalization)
                
                
                # 
                # process the metadata
                # 
                if True:
                    (
                        conditions_by_wig_fingerprint,
                        covariates_by_wig_fingerprint_list,
                        interactions_by_wig_fingerprint_list,
                        ordering_metadata,
                    ) = tnseq_tools.CombinedWigMetadata.read_condition_data(
                        metadata_path, covars, interactions, column_name_for_condition=group_by,
                    )
                    
                    # 
                    # [Condition] in the order of files in combined wig
                    # 
                    conditions_from_comwig = [
                        conditions_by_wig_fingerprint.get(f, "FLAG-UNMAPPED-CONDITION-IN-WIG")
                            for f in filenames_in_comb_wig
                    ]
                    # 
                    # [Covariate] in the order of files in combined wig
                    # 
                    try:
                        covariates_from_comwig = [
                            [covars_by_file.get(f, "?") for f in filenames_in_comb_wig]
                            for covars_by_file in covariates_by_wig_fingerprint_list
                        ]
                    except KeyError:
                        for f in filenames_in_comb_wig:
                            if f not in covariates_by_wig_fingerprint_list[0]:
                                logging.log(f)
                        logging.error("Error: Covariates not found for sample:")
                    # 
                    # [Interaction] in the order of files in combined wig
                    # 
                    try:
                        interactions_from_wigs = [
                            [covars_by_file.get(f, "?") for f in filenames_in_comb_wig]
                            for covars_by_file in interactions_by_wig_fingerprint_list
                        ]
                    except KeyError:
                        for f in filenames_in_comb_wig:
                            if f not in interactions_by_wig_fingerprint_list[0]:
                                logging.log(f)
                        logging.error("Error: Interaction var not found for sample")
                    
                    metadata = transit_tools.get_samples_metadata(metadata_path)
                    condition_names = metadata["Condition"]  # original Condition names for each sample, as ordered list
                    wig_fingerprints = metadata["Filename"]

                    # this is the new way to filter samples, where --include/exclude-conditions refers to Condition column in metadata, regardless of whether --condition was specified
                    (
                        data,
                        wig_fingerprints,
                        condition_names,
                        filtered_conditions_from_comwig,
                        filtered_covariates_from_comwig,
                        interactions_from_wigs,
                    ) = transit_tools.filter_wigs_by_conditions2(
                        data,
                        wig_fingerprints=wig_fingerprints,
                        condition_names=condition_names,  # original Condition column in samples metadata file
                        included_cond=included_conditions,
                        excluded_cond=excluded_conditions,
                        conditions=conditions_from_comwig,  # user might have specified a column other than Condition
                        covariates=covariates_from_comwig,
                        interactions=interactions_from_wigs,
                    )
                    wig_fingerprint_set = set(wig_fingerprints)

                    # show the samples associated with each condition (and covariates or interactions, if defined), and count samples in each cross-product of vars
                    selected_metadata_headers = misc.no_duplicates([group_by] + covars + interactions)
                    headers_to_vals = {}
                    headers_to_vals[group_by] = misc.no_duplicates(filtered_conditions_from_comwig)
                    for index_of_covariate, each_covariate_header in enumerate(covars):
                        # TODO: this indexing looks suspicious to me at the moment, maybe test it
                        headers_to_vals[each_covariate_header] = misc.no_duplicates(filtered_covariates_from_comwig[index_of_covariate])
                    for index_of_interaction, each_interaction_header in enumerate(interactions):
                        headers_to_vals[each_interaction_header] = misc.no_duplicates(interactions_from_wigs[index_of_interaction])
                    
                    metadata_header_index_to_mapping_of_value_to_wig_fingerprints = [conditions_by_wig_fingerprint] + covariates_by_wig_fingerprint_list + interactions_by_wig_fingerprint_list
                    for index, each_selected_header in enumerate(selected_metadata_headers):
                        logging.log("\nsamples for Condition/Covariate/Interaction: %s" % selected_metadata_headers[index])
                        value_to_wig_fingerprints = misc.invert_dict(metadata_header_index_to_mapping_of_value_to_wig_fingerprints[index])
                        for each_metadata_value, each_wig_fingerprints in value_to_wig_fingerprints.items():
                            samples = list(wig_fingerprint_set.intersection(set(each_wig_fingerprints)))
                            if each_metadata_value in headers_to_vals.get(each_selected_header, []):
                                logging.log("        %s: %s" % (each_metadata_value, " ".join(samples)))
                    
                    logging.log("\nsamples in cross-product:")
                    any_empty = transit_tools.expand_var([], selected_metadata_headers, metadata_header_index_to_mapping_of_value_to_wig_fingerprints, headers_to_vals, set(wig_fingerprint_set), enable_logging=(not disable_logging))
                    if any_empty:
                        logging.warn("ZINB requires samples in all combinations of conditions; the fact that one is empty could result in Model Errors")
                
                # 
                # process gene data 
                # 
                if True:
                    genes = tnseq_tools.read_genes(annotation_path)
                    ta_site_index_map = {
                        each_ta_site: index
                            for index, each_ta_site in enumerate(sites) 
                    }
                    rv_site_indexes_map = tnseq_tools.rv_site_indexes_map(genes, ta_site_index_map, n_terminus=n_terminus, c_terminus=c_terminus)
                    stats_by_rv, stat_group_names = transit_tools.get_stats_by_rv(data, rv_site_indexes_map, genes, filtered_conditions_from_comwig, interactions_from_wigs, winz)
                    log_z_perc_by_replicate, non_zero_mean_by_replicate = Method.global_stats_for_rep(data)
            # 
            # Main computation
            # 
            if True:
                logging.log("Running ZINB")
                gene_name_to_p_value, gene_name_to_adj_p_value, run_status = Method.calculate(
                    group_by=group_by,
                    interactions=interactions,
                    covars=covars,
                    winz=winz,
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
                        group_by == "Condition"
                        and len(included_conditions) > 0
                        and len(excluded_conditions) == 0
                    ):
                        cond_diff = included_conditions.index(
                            c1
                        ) - included_conditions.index(c2)
                        ## Order by interaction, if stat belongs to same condition
                        if cond_diff == 0 and i1 is not None and i2 is not None:
                            return ordering_metadata["interaction"].index(i1) - ordering_metadata[
                                "interaction"
                            ].index(i2)
                        return cond_diff

                    ## Order by samples metadata, if include flag not provided.
                    cond_diff = ordering_metadata["condition"].index(c1) - ordering_metadata[
                        "condition"
                    ].index(c2)
                    if cond_diff == 0 and i1 is not None and i2 is not None:
                        return ordering_metadata["interaction"].index(i1) - ordering_metadata[
                            "interaction"
                        ].index(i2)
                    return cond_diff
                
                import functools
                ordered_stat_group_names = sorted(stat_group_names, key=functools.cmp_to_key(order_stats))
                headers_stat_group_names = [ x.replace(SEPARATOR, "_") for x in ordered_stat_group_names ]
            
            # 
            # create column_names, rows, extra_info
            # 
            if True:
                only_have_one_lfc_column = len(ordered_stat_group_names) == 2
                have_gene_name_column = gene_name_to_description != None
                # 
                # rows
                #
                rows = [] 
                for gene in genes:
                    gene_name = gene["rv"]
                    means_per_group = [
                        stats_by_rv[gene_name]["mean"][group]
                            for group in ordered_stat_group_names
                    ]
                    if only_have_one_lfc_column:
                        # TODO: still need to adapt this to use --ref if defined
                        log_fold_changes = [
                            numpy.math.log(
                                (means_per_group[1] + pseudocount) / (means_per_group[0] + pseudocount),
                                2
                            ) 
                        ]  
                    else:
                        if len(refs) == 0:
                            grand_means = numpy.mean(means_per_group)  # grand mean across all conditions
                        else:
                            grand_means = numpy.mean(
                                [stats_by_rv[gene_name]["mean"][group] for group in refs]
                            )
                        log_fold_changes = [numpy.math.log((each_mean + pseudocount) / (grand_means + pseudocount), 2) for each_mean in means_per_group]
                    
                    row = [
                        gene_name,
                        gene["gene"],
                        str(len(rv_site_indexes_map[gene_name])),
                        *[
                            "%0.1f" % stats_by_rv[gene_name]["mean"][each_group]
                                for each_group in ordered_stat_group_names
                        ],
                        *[
                            "%0.3f" % each_lfc
                                for each_lfc in log_fold_changes
                        ],
                        *[
                            "%0.1f" % stats_by_rv[gene_name]["nz_mean"][each_group]
                                for each_group in ordered_stat_group_names
                        ],
                        *[
                            "%0.2f" % stats_by_rv[gene_name]["nz_perc"][each_group]
                                for each_group in ordered_stat_group_names
                        ],
                        gene_name_to_p_value[gene_name],
                        gene_name_to_adj_p_value[gene_name],
                        run_status[gene_name],
                    ]
                    
                    if have_gene_name_column:
                        gene_name = gene_name_to_description.get(gene_name, "?")
                        row.append(gene_name)
                    
                    rows.append(row)
                # 
                # column_names
                # 
                if only_have_one_lfc_column:
                    lfc_names = [ "Log 2 FC" ]
                else:
                    lfc_names = [ "Log 2 FC "+each_name for each_name in headers_stat_group_names ]
                
                column_names = [
                    "Rv",
                    "Gene",
                    "TA Sites",
                    *[
                        "Mean "+each_name
                            for each_name in headers_stat_group_names
                    ],
                    *lfc_names,
                    *[
                        "Non Zero Mean "+each_name
                            for each_name in headers_stat_group_names
                    ],
                    *[
                        "Non Zero Percentage "+each_name
                            for each_name in headers_stat_group_names
                    ],
                    "P Value",
                    "Adj P Value",
                    "Status",
                ]
                if have_gene_name_column:
                    column_names.append("Gene Name")
                
                # 
                # extra_info
                # 
                extra_info = dict(
                    files=dict(
                        combined_wig=combined_wig_path,
                        annotation_path=annotation_path,
                        prot_table=prot_table_path,
                    ),
                    parameters=dict(
                        normalization=normalization,
                        n_terminus=n_terminus,
                        c_terminus=c_terminus,
                        pseudocount=pseudocount
                    )
                )
            
            # 
            # write output
            # 
            logging.log(f"Adding File: {output_path}")
            transit_tools.write_result(
                path=output_path, # path=None means write to STDOUT
                file_kind=Method.identifier,
                rows=rows,
                column_names=column_names,
                extra_info=extra_info,
            )
            
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
        from pytransit.specific_tools.transit_tools import DataFrame, IntVector, FloatVector, StrVector
        import statsmodels.stats.multitest

        count = 0
        
        p_values, gene_names, status = [], [], []
        r_zinb_signif = Method.def_r_zinb_signif()
        logging.log("Running analysis...")
        if winz:
            logging.log("Winsorizing insertion count data")

        logging.log(f"Grouping by: {group_by}")

        comp1a = "1+cond"
        comp1b = "1+cond"

        # include cond in mod0 only if testing interactions
        comp0a = "1" if len(interactions) == 0 else "1+cond" # >>> "1+cond"
        comp0b = "1" if len(interactions) == 0 else "1+cond" # >>> "1+cond"
        for each_interaction in interactions: # >>>  [ "KO_5849" ]
            comp1a += "*" + each_interaction
            comp1b += "*" + each_interaction
            comp0a += "+" + each_interaction
            comp0b += "+" + each_interaction
        for each_covariate in covars: # >>> [ ]
            comp1a += "+" + each_covariate
            comp1b += "+" + each_covariate
            comp0a += "+" + each_covariate
            comp0b += "+" + each_covariate
        zinb_mod1 = "cnt~%s+offset(log(non_zero_mean))|%s+offset(logit_z_perc)" % (comp1a, comp1b)
        zinb_mod0 = "cnt~%s+offset(log(non_zero_mean))|%s+offset(logit_z_perc)" % (comp0a, comp0b)

        nb_mod1 = "cnt~%s" % (comp1a)
        nb_mod0 = "cnt~%s" % (comp0a)
        to_r_float_or_str_vec = lambda xs: (
            FloatVector([float(x) for x in xs])   if misc.str_is_float(xs[0])   else StrVector(xs)
        )

        for gene in genes:
            count += 1
            gene_name = gene["rv"]
            
            # 
            # Single gene case for debugging
            # 
            if cli_args.gene:
                gene_name = None
                if cli_args.gene in rv_site_indexes_map:
                    gene_name = cli_args.gene
                else:
                    for each_gene in genes:
                        if each_gene["gene"] == cli_args.gene:
                            gene_name = each_gene["rv"]
                            break
                if not gene_name:
                    logging.error("Cannot find gene: {0}".format(cli_args.gene))
            if debugging_enabled:
                logging.log("======================================================================")
                logging.log(gene["rv"] + " " + gene["gene"])
            
            
            if len(rv_site_indexes_map[gene_name]) <= 1:
                status.append("TA sites <= 1, not analyzed")
                p_values.append(1)
            else:
                norm_data = list(
                    map(lambda wigData: wigData[rv_site_indexes_map[gene_name]], data)
                )
                
                if winz:
                    norm_data = tnseq_tools.winsorize(norm_data)
                (
                    read_counts,
                    condition,
                    covars_data,
                    interactions_data,
                    non_zero_mean,
                    logit_z_perc,
                ) = Method.melt_data(
                    read_counts_for_rv=norm_data,
                    conditions=filtered_conditions_from_comwig,
                    covariates=filtered_covariates_from_comwig,
                    interactions=interactions_from_wigs,
                    non_zero_mean_by_replicate=non_zero_mean_by_replicate,
                    log_z_perc_by_replicate=log_z_perc_by_replicate,
                )
                if numpy.sum(read_counts) == 0:
                    status.append(
                        "pan-essential (no counts in all conditions) - not analyzed"
                    )
                    p_values.append(1)
                else:
                    df_args = {
                        "cnt": IntVector(read_counts),
                        "cond": to_r_float_or_str_vec(condition),
                        "non_zero_mean": FloatVector(non_zero_mean),
                        "logit_z_perc": FloatVector(logit_z_perc),
                    }
                    ## Add columns for covariates and interactions if they exist.
                    df_args.update(
                        list(
                            map(
                                lambda t_ic: (
                                    t_ic[1],
                                    to_r_float_or_str_vec(covars_data[t_ic[0]]),
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
                                    to_r_float_or_str_vec(interactions_data[t_ic[0]]),
                                ),
                                enumerate(interactions),
                            )
                        )
                    )

                    melted = DataFrame(df_args)
                    # r_args = [IntVector(read_counts), StrVector(condition), melted, map(lambda x: StrVector(x), covars), FloatVector(non_zero_mean), FloatVector(logit_z_perc)] + [True]
                    debugging = debugging_enabled or cli_args.gene
                    pval, msg = r_zinb_signif(
                        melted, zinb_mod1, zinb_mod0, nb_mod1, nb_mod0, debugging
                    )
                    status.append(msg)
                    p_values.append(float(pval))
                if debugging_enabled or cli_args.gene:
                    logging.log(
                        "P Value for Gene {0}: {1}, status: {2}".format(
                            gene_name, p_values[-1], status[-1]
                        )
                    )
                if cli_args.gene:
                    logging.log("Ran for single gene. Exiting...")
                    sys.exit(0)
            gene_names.append(gene_name)
            # Update progress
            percentage = (100.0 * count / len(genes))
            text = "Running ZINB Method... %5.1f%%" % percentage
            parameter_panel.progress_update(text, percentage)

        p_values = numpy.array(p_values)
        mask = numpy.isfinite(p_values)
        adj_p_values = numpy.full(p_values.shape, numpy.nan)
        adj_p_values[mask] = statsmodels.stats.multitest.fdrcorrection(p_values)[1]  # BH, alpha=0.05

        gene_name_to_p_value, gene_name_to_adj_p_value, status_map = {}, {}, {}
        for gene_index, gene_name in enumerate(gene_names):
            gene_name_to_p_value[gene_name], gene_name_to_adj_p_value[gene_name], status_map[gene_name] = p_values[gene_index], adj_p_values[gene_index], status[gene_index]
        
        return gene_name_to_p_value, gene_name_to_adj_p_value, status_map

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
        read_counts_for_rv,
        conditions,
        covariates,
        interactions,
        non_zero_mean_by_replicate,
        log_z_perc_by_replicate,
    ):
        rv_sites_length = len(read_counts_for_rv[0])
        repeat_and_flatten = lambda xs: numpy.repeat(xs, rv_sites_length)

        return [
            numpy.concatenate(read_counts_for_rv).astype(int),
            repeat_and_flatten(conditions),
            list(map(repeat_and_flatten, covariates)),
            list(map(repeat_and_flatten, interactions)),
            repeat_and_flatten(non_zero_mean_by_replicate),
            repeat_and_flatten(log_z_perc_by_replicate),
        ]
    
    @staticmethod
    def def_r_zinb_signif():
        if not HAS_R: return lambda *args, **kwargs: Exception("You're getting this if debugging is on and you dont have R")
        from pytransit.specific_tools.transit_tools import r, globalenv
        r(
            """
            zinb_signif = function(
                df,
                zinbMod1,
                zinbMod0,
                nbMod1,
                nbMod0,
                DEBUG = F
            ) {
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
                  newvec$cnt = 1 # note: non_zero_mean and NZperc are copied from last dataset in condition
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
        
@transit_tools.ResultsFile
class ResultFileType1:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier)
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier,
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier,
                    heading=self.comments_string or misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[
                        # HANDLE_THIS
                    ],
                ).Show(),
            })
        )
        
        # 
        # read in data
        # 
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        self.values_for_result_table.update(self.extra_data.get("parameters", {}))
        
        # 
        # get summary stats
        #
        self.values_for_result_table.update({
            # HANDLE_THIS (additional summary_info for results table)
            # examples:
                # f"Gene Count": len(self.rows),
                # f"Adj P Value < {Method.significance_threshold}": len([
                #     1 for each in self.rows
                #         if each.get("Adj P Value", 0) < Method.significance_threshold 
                # ]),
        })
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()

