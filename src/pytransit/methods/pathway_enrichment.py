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
from statsmodels.stats import multitest

from pytransit.specific_tools import logging, gui_tools
from pytransit.specific_tools.transit_tools import wx
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.parameter_panel import panel,progress_update, set_instructions

from pytransit.specific_tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import csv, misc
from pytransit.specific_tools.transit_tools import wx, basename, HAS_R, FloatVector, DataFrame, StrVector
from pytransit.globals import gui, cli, root_folder, debugging_enabled
from pytransit.components import file_display, results_area, parameter_panel, panel_helpers

@misc.singleton
class Method:
    name = "Pathway Enrichment"
    identifier  = name.replace(" ", "")
    cli_name    = name.replace(" ", "_").lower()
    menu_name   = f"{identifier} - Perform {name} analysis"
    description = f"""Perform {name} analysis"""
    rows = []
    
    inputs = LazyDict(
        resampling_file = None,
        associations_file = None,
        pathways_file = None,
        output_path= None,
        organism_pathway = None,
        pval_col = -2,
        qval_col = -1,
        lfc_col = 6,
        num_permutations = 10000,
        enrichment_exponent = 0, # for GSEA
        pseudocount = 2
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

    #-Pval_col <int>    : indicate column with *raw* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for sorting) (default: -2)
    #-Qval_col <int>    : indicate column with *adjusted* P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for significant cutoff) (default: -1)
    #-LFC_col <int>     : indicate column with log2FC (starting with 0; can also be negative, i.e. -1 means last col) (used for ranking genes by SLPV or LFC) (default: 6)

    usage_string = f"""{console_tools.subcommand_prefix} pathway_enrichment <resampling_file> <associations> <pathways> <output_file> [-M <FET|GSEA|GO>] [-PC <int>] [-ranking SLPV|LFC] [-p <float>] [-Nperm <int>] [-Pval_col <int>] [-Qval_col <int>]  [-LFC_col <int>]

        Optional parameters:
        -M FET|GSEA|ONT:     method to use, FET for Fisher's Exact Test (default), GSEA for Gene Set Enrichment Method (Subramaniam et al, 2005), or ONT for Ontologizer (Grossman et al, 2007)

        for GSEA...
        -ranking SLPV|LFC  : SLPV is signed-log-p-value (default); LFC is log2-fold-change from resampling 
        -p <float>         : exponent to use in calculating enrichment score; recommend trying 0 or 1 (as in Subramaniam et al, 2005)
        -Nperm <int>       : number of permutations to simulate for null distribution to determine p-value (default=10000)
        for FET...
        -PC <int>          :  pseudo-counts to use in calculating p-value based on hypergeometric distribution (default=2)
    """.replace("\n        ", "\n")
    
    @gui.add_menu("Method", "himar1", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)
    
    @gui.add_menu("Method", "tn5", menu_name)
    def on_menu_click(event):
        Method.define_panel(event)

    def call_from_results_panel(self, results_file):
        self.inputs.resampling_file = results_file
        self.define_panel()

    def create_default_pathway_button(self,panel, sizer, *, button_label, tooltip_text=""):
        import csv
        COG_orgs = []
        with open(root_folder+"src/pytransit/data/cog-20.org.tsv") as file_obj:
            reader_obj = csv.reader(file_obj)
            for row in reader_obj:
                COG_orgs.append(row[1])
        row_sizer = wx.BoxSizer(wx.HORIZONTAL)
        if True:
            # 
            # tooltip
            # 
            if tooltip_text:
                from pytransit.components.icon import InfoIcon
                row_sizer.Add(
                    InfoIcon(panel, wx.ID_ANY, tooltip=tooltip_text),
                    0,
                    wx.ALIGN_CENTER_VERTICAL,
                    gui_tools.default_padding,
                )
            # 
            # button
            # 
            if True:
                popup_button = wx.Button(
                    panel,
                    wx.ID_ANY,
                    button_label,
                    wx.DefaultPosition,
                    wx.DefaultSize,
                    0,
                )
                # whenever the button is clicked, popup
                organism_pathway_text = None
                organism_pathway = None
                @gui_tools.bind_to(popup_button, wx.EVT_BUTTON)
                def when_button_clicked(*args,**kwargs):
                    nonlocal organism_pathway
                    win = wx.Dialog(panel,wx.FRAME_FLOAT_ON_PARENT)
                    popup_sizer = wx.BoxSizer(wx.VERTICAL)
                    win.SetSizer(popup_sizer)

                    pathway_label_text= wx.StaticText(win, wx.ID_ANY, label="Select A Pathway Type : ", style=wx.ALIGN_LEFT)
                    popup_sizer.Add(pathway_label_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                    pathway_type = wx.ComboBox(win,choices = ["Sanger", "COG" ,"GO", "KEGG"])
                    popup_sizer.Add(pathway_type,wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)

                    select_btn = wx.Button(win, wx.ID_OK, label = "Select", size = (50,20), pos = (75,50))
                    popup_sizer.Add(select_btn,wx.EXPAND, gui_tools.default_padding)

                    win.Layout()
                    popup_sizer.Fit(win)
                    selected_path = win.ShowModal()

                    if selected_path == wx.ID_OK:
                        pathway_type_selected = pathway_type.GetValue()
                        organism_label_text= wx.StaticText(win, wx.ID_ANY, label="Select An Organism : ", style=wx.ALIGN_LEFT)
                        popup_sizer.Add(organism_label_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                        if pathway_type_selected== "COG":                           
                            organism = wx.ComboBox(win,choices = sorted(COG_orgs))               
                        elif pathway_type_selected== "KEGG":
                            organism = wx.ComboBox(win,choices = ["H37Rv"])
                        elif pathway_type_selected== "Sanger":
                            organism = wx.ComboBox(win,choices = ["H37Rv"])
                        else:
                            organism = wx.ComboBox(win,choices = ["H37Rv", "Smeg"])

                        popup_sizer.Add(organism,wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                        ok_btn = wx.Button(win, wx.ID_OK, label = "Ok", size = (50,20), pos = (75,50))
                        popup_sizer.Add(ok_btn,wx.EXPAND, gui_tools.default_padding)

                        win.Layout()
                        popup_sizer.Fit(win)
                        res = win.ShowModal()
                        if res == wx.ID_OK:
                            organism_pathway = "-".join([organism.GetValue(),pathway_type_selected])
                            organism_pathway_text.SetLabel(basename(organism_pathway or ""))
                        win.Destroy()
                    

            row_sizer.Add(popup_button, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)

            organism_pathway_text= wx.StaticText(panel, wx.ID_ANY, label="", style=wx.ALIGN_LEFT)
            row_sizer.Add(organism_pathway_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
        
        sizer.Add(row_sizer, 0, wx.ALL | wx.ALIGN_CENTER_HORIZONTAL, gui_tools.default_padding)
        return lambda *args, **kwargs: organism_pathway

    def define_panel(self,_=None):
        from pytransit.components import panel_helpers 
        with panel_helpers.NewPanel() as (panel, main_sizer):
            set_instructions(
                method_short_text=self.name,
                method_long_text="",
                method_specific_instructions="""
                Pathway Enrichment Analysis provides a method to identify enrichment of functionally-related genes among those that are conditionally essential (i.e. significantly more or less essential between two conditions). The analysis is typically applied as post-processing step to the hits identified by a comparative analysis, such as resampling. Several analytical method are provided: Fisher’s exact test (FET, hypergeometric distribution), GSEA (Gene Set Enrichment Analysis) by Subramanian et al (2005), and Ontologizer. 


                    1. If you have selected this method from the menu bar, ensure you select a resampling file from the Select 
                       Input File button below

                    2. Add in Associaiton/Pathway Information
                        a. Choose from one of our provided associations and pathways using the "Select from Provided Files" OR
                        b. Select your own using the Select Custom Associations and Select Custom Pathways Buttons

                    3. Select Pathway Enrichement Method. 
                        * Ontologizer is best suited for GO Terms

                    4. [Optional] Adjust parameters
                    
                    5. Click Run
                """.replace("\n                    ","\n"),
            )
            self.value_getters = LazyDict()

            if Method.inputs.resampling_file == None:
                self.value_getters.resampling_file = panel_helpers.create_file_input(panel, main_sizer, 
                    button_label="Select Resampling File", 
                    tooltip_text="The resampling file is the one obtained after using the resampling method in Transit. (It is a tab separated file with 11 columns.) GSEA method makes usage of the last column (adjusted P-value)", 
                    popup_title="Select File with Hits",
                    allowed_extensions='All files (*.*)|*.*'
                )

            self.value_getters.associations_file = panel_helpers.create_file_input(panel, main_sizer, 
                button_label="Select Custom Associations File", 
                tooltip_text="This is a tab-separated text file with 2 columns: pathway id, and pathway name. If a gene is in multiple pathways, the associated ids should be listed on separate lines. It is OK if there are no associations listed for some genes. Important: if pathways are hierarchical, you should expand this file to explicitly include associations of each gene with all parent nodes. Files with GO term associations will have to be pre-processed this way too.", 
                popup_title="Select Associations File",
                allowed_extensions='All files (*.*)|*.*'
            )

            self.value_getters.pathways_file = panel_helpers.create_file_input(panel, main_sizer, 
                button_label="Select Custom Pathways File", 
                tooltip_text="This is a tab-separated text file with 2 columns: pathway id, and pathway name.", 
                popup_title="Select Pathways File",
                allowed_extensions='All files (*.*)|*.*'
            )

            self.value_getters.organism_pathway =  self.create_default_pathway_button(panel, main_sizer, 
                button_label="Select from Provided Files", 
                tooltip_text="We have a few Associaiton and Pathway files pre-loaded for for your use. When this button is clicked, a pop-up will appear that will allow you to select a Pathway type and organism", 
            )
    
            self.value_getters.method = panel_helpers.create_choice_input(panel, main_sizer,
                label = "Method",
                options= ["FET", "GSEA", "ONT"],
                tooltip_text = "Method to use, FET for Fisher's Exact Test (default), GSEA for Gene Set Enrichment Method (Subramaniam et al, 2005), or ONT for Ontologizer (Grossman et al, 2007)"
            )

            self.value_getters.ranking = panel_helpers.create_choice_input(panel, main_sizer,
                label = "Ranking",
                options= ["SLPV", "LFC"],
                tooltip_text="SLPV is signed-log-p-value (default); LFC is log2-fold-change from resampling"
            )

            self.value_getters.pval_col            = panel_helpers.create_int_getter(panel, main_sizer, label_text="P-Value Column Index", default_value=-2, tooltip_text="Column index with raw P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for sorting) (default: -2, i.e. second-to-last column)")
            self.value_getters.qval_col            = panel_helpers.create_int_getter(panel, main_sizer, label_text="Q-Value Column Index", default_value=-1, tooltip_text="Column index with adjusted P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for significant cutoff) (default: -1)")
            self.value_getters.lfc_col             = panel_helpers.create_int_getter(panel, main_sizer, label_text="LFC Column Index", default_value=6, tooltip_text="Column index with log2FC (starting with 0; can also be negative, i.e. -1 means last col) (used for ranking genes by SLPV or LFC) (default: 6)")
                
            self.value_getters.enrichment_exponent = panel_helpers.create_int_getter(  panel, main_sizer, label_text="Enrichment Exponent",    default_value=0,      tooltip_text="Exponent to use in calculating enrichment score; recommend trying 0 or 1 (as in Subramaniam et al, 2005)")
            self.value_getters.num_permutations    = panel_helpers.create_int_getter(  panel, main_sizer, label_text="Number of Permutations", default_value=10000,  tooltip_text="Number of permutations to simulate for null distribution to determine p-value")
            self.value_getters.pseudocount         = panel_helpers.create_pseudocount_input(panel, main_sizer, default_value=2)
            
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function = self.from_gui)


    @staticmethod
    def from_gui(frame):       
        # 
        # get annotation
        # 
        #Method.inputs.annotation_path = gui.annotation_path
        # 
        # setup custom inputs
        #        
        for each_key, each_getter in Method.value_getters.items():
            try:
                Method.inputs[each_key] = each_getter()
            except Exception as error:
                logging.error(f'''Failed to get value of "{each_key}" from GUI:\n{error}''')


        if Method.inputs.organism_pathway != None:
            organism,pathway = Method.inputs.organism_pathway.split("-")
            if pathway == "COG":
                try:
                    import requests
                    URL = "https://orca1.tamu.edu/essentiality/transit/COG2020/"+organism+"_COG_20_roles.associations.txt"
                    response = requests.get(URL)
                    open(root_folder+"src/pytransit/data/"+organism+"_COG_20_roles.associations.txt", "wb").write(response.content)
                except requests.exceptions.ConnectionError:
                    logging.error("Please Connect to the Internet to get this COG files for "+organism)
                    #sys.exit()
                
                Method.inputs.associations_file = root_folder+"src/pytransit/data/"+organism+"_COG_20_roles.associations.txt"
                Method.inputs.pathways_file = root_folder+"src/pytransit/data/COG_20_roles.txt"

            elif Method.inputs.organism_pathway =="H37Rv-Sanger":
                logging.log("Loading in H37Rv Associations for Sanger Pathways")
                Method.inputs.associations_file = root_folder+"src/pytransit/data/H37Rv_sanger_roles.dat"
                Method.inputs.pathways_file = root_folder+"src/pytransit/data/sanger_roles.dat"
            elif Method.inputs.organism_pathway =="H37Rv-KEGG":
                logging.log("Loading in H37Rv Associations for KEGG Pathways")
                Method.inputs.associations_file = root_folder+"src/pytransit/data/H37Rv_KEGG_roles.txt"
                Method.inputs.pathways_file = root_folder+"src/pytransit/data/KEGG_roles.txt"
            elif Method.inputs.organism_pathway =="H37Rv-GO" and Method.inputs.method == "FET":
                logging.log("Loading in H37Rv Associations for GO Pathways")
                Method.inputs.associations_file = root_folder+"src/pytransit/data/H37Rv_GO_terms.txt"
                Method.inputs.pathways_file = root_folder+"src/pytransit/data/GO_term_names.dat"
            elif Method.inputs.organism_pathway =="H37Rv-GO" and Method.inputs.method == "ONT":
                logging.log("Loading in H37Rv Associations for GO Pathways")
                Method.inputs.associations_file = root_folder+"src/pytransit/data/H37Rv_GO_terms.txt"
                Method.inputs.pathways_file = root_folder+"src/pytransit/data/gene_ontology.1_2.3-11-18.obo"
            elif Method.inputs.organism_pathway =="Smeg-GO":
                logging.log("Loading in Smeg Associations for GO Pathways")
                Method.inputs.associations_file = root_folder+"src/pytransit/data/smeg_GO_terms.txt"
                Method.inputs.pathways_file = root_folder+"src/pytransit/data/GO_term_names.dat"  

        Method.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Method.cli_name}_output.txt",
            output_extensions='Common output extensions (*.tsv,*.dat,*.txt,*.out)|*.tsv;*.dat;*.txt;*.out;|\nAll files (*.*)|*.*',
        )


        # 
        # validate input files have been selected
        # 
        assert Method.inputs.resampling_file != None, "Please select a resampling file"
        assert Method.inputs.associations_file != None, "Please select an associations file from the pre-provided options or upload your own"
        assert Method.inputs.pathways_file != None, "Please select a pathways file from the pre-provided options or upload your own"

        return Method

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, exactly=4)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)

        # save the data
        Method.inputs.update(dict(
            resampling_file = args[0],
            associations_file = args[1],
            pathways_file = args[2],
            output_path=args[3],
            method = kwargs.get("M", "FET"),
            pval_col = int(kwargs.get("Pval_col", Method.inputs.pval_col)),
            qval_col = int(kwargs.get("Qval_col", Method.inputs.qval_col)),
            ranking = kwargs.get("ranking", "SLPV"),
            lfc_col = int(kwargs.get("LFC_col", Method.inputs.lfc_col)),
            enrichment_exponent = int(kwargs.get("p", "1")),
            num_permutations = int(kwargs.get("Nperm", Method.inputs.num_permutations)),
            pseudocount = int(kwargs.get("PC", "2")),
        ))
        
        Method.Run()
        
    def Run(self):
        self.rows = []
        with gui_tools.nice_error_log:
            from pytransit.specific_tools import stat_tools
            logging.log(f"Starting {Method.identifier} analysis")
            start_time = time.time()
            
            #checking validation of inputs
            if self.inputs.method == "FET":
                self.hit_summary = self.fisher_exact_test()
                file_output_type = Method.identifier+"FET"
                file_columns = [
                        "Pathway",
                        "Total Genes", 
                        "Genes In Path",
                        "Significant Genes",
                        "Significent Genes In Path",
                        "Expected", 
                        "K Plus PC",
                        "Number Adjusted By PC",
                        "Enrichment" , 
                        "P Value", 
                        "Adj P Value", 
                        "Description", 
                        "Genes"
                    ]
            elif self.inputs.method == "GSEA":
                up,down = self.GSEA()
                self.hit_summary = str(up)+str(" Siginificant Pathways for Conditional Essential Genes, ") + str(down) + str(" Siginificant Pathways for Conditional Non-Essential Genes, ")
                file_output_type = Method.identifier+"GSEA"
                file_columns = [
                        "Pathway",
                        "Pathway Description"
                        "Genes in Path", 
                        "Mean Rank",
                        "Enrichment" , 
                        "P Value", 
                        "Adj P Value", 
                        "Genes"
                    ]
            elif self.inputs.method == "ONT":
                self.hit_summary = self.Ontologizer()
                file_output_type = Method.identifier+"ONT"
                file_columns = [
                        "Pathway",
                        "Total Genes", 
                        "Genes In Path",
                        "Significant Genes",
                        "Significent Genes In Path",
                        "Expected", 
                        "K Plus PC",
                        "Number Adjusted By PC",
                        "Enrichment" , 
                        "P Value", 
                        "Adj P Value", 
                        "Description", 
                        "Genes"
                    ]
            else:
                self.inputs.method = "Not a valid method"
                progress_update("Not a valid method", 100)
        
        # 
        # write output
        # 
        logging.log(f"Adding File: {self.inputs.output_path}")
        # 
        # write to file
        # 
        transit_tools.write_result(
            path=self.inputs.output_path, # path=None means write to STDOUT
            file_kind=file_output_type,
            rows=self.rows,
            column_names=file_columns,
            extra_info=dict(
                parameters=dict(
                    method = self.inputs.method,
                    ranking = self.inputs.ranking,
                    enrichment_exponent = self.inputs.enrichment_exponent,
                    num_permutations = self.inputs.num_permutations,
                    pseudocount = self.inputs.pseudocount,
                    hit_summary = self.hit_summary
                ),
            ),
        )
        logging.log(f"Finished {Method.identifier} analysis in {time.time() - start_time:0.1f}sec")
        results_area.add(self.inputs.output_path)
        
    def read_resampling_file(self, filename):
        logging.log("Reading in Resampling File", filename)
        genes, hits, standardized_headers= [], [], []
        with open(filename) as file:
            for line in file:
                if line[0] == "#":
                    headers = line.split("\t")
                    if len (headers)>2:
                        standardized_headers = [misc.pascal_case_with_spaces(col) for col in headers]
                    continue
                w = line.rstrip().split("\t")
                genes.append(w)
                qval = float(w[self.inputs.qval_col])
                if qval < 0.05:
                    hits.append(w[0])
        
        return genes, hits, standardized_headers
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
        logging.log("Reading in Pathways File")
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
        powers = [math.pow(abs(float(x)), float(p)) for x in Ascores]
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

    # during initialization, self.inputs.resampling_file etc have been set, and self.inputs.output_path has been opened
    
    def GSEA(self):  
        logging.log("Running GSEA")
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

        # rank by SLPV=sign(LFC)*log10(pval)
        # note: genes with lowest p-val AND negative LFC have highest scores (like positive correlation)
        # there could be lots of ties with pval=0 or 1, but so be it (there are probably fewer such ties than with Qvals) (randomize order below)
        pairs = []  # pair are: rv and score (SLPV)
        for w in data:
            orf = w[0]
            if self.inputs.ranking == "SLPV":
                p_value = float(w[self.inputs.pval_col])
                LFC = float(w[self.inputs.lfc_col])
                SLPV = (-1 if LFC < 0 else 1) * math.log(p_value + 0.000001, 10)
                pairs.append((orf, SLPV))
            elif self.inputs.ranking == "LFC":
                # lfc_col = headers.index("log2FC")
                LFC = float(w[self.inputs.lfc_col])
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
        num_permutations = self.inputs.num_permutations
        results, Total = [], len(terms)
        for i, term in enumerate(terms):
            sys.stdout.flush()
            orfs = terms2orfs.get(term, [])
            num_genes_in_pathway = len(orfs)
            if num_genes_in_pathway < 2:
                continue  # skip pathways with less than 2 genes
            mr = self.mean_rank(orfs, orfs2rank)
            es = self.enrichment_score(
                orfs, orfs2rank, orfs2score, p=self.inputs.enrichment_exponent,
            )  # always positive, even if negative deviation, since I take abs
            higher = 0
            for n in range(num_permutations):
                perm = random.sample(
                    allgenes, num_genes_in_pathway
                )  # compare to enrichment score for random sets of genes of same size
                e2 = self.enrichment_score(perm, orfs2rank, orfs2score, p=self.inputs.enrichment_exponent)
                if e2 > es:
                    higher += 1
                if n > 100 and higher > 10:
                    break  # adaptive: can stop after seeing 10 events (permutations with higher ES)
            pval = higher / float(n)
            # vals = [
            #     "#",
            #     term,
            #     #num_genes_in_pathway,
            #     mr,
            #     es,
            #     pval,
            #     ontology.get(term, "?"),
            # ]
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

        for term,mr,es,pval,qval in results:
            if qval<0.05 and mr<n2: 
                self.rows.append("#   %s %s (mean_rank=%s)" % (term,ontology.get(term,"?"),mr)
                )
    
        for term,mr,es,pval,qval in results:
            if qval<0.05 and mr>n2: 
                self.rows.append("#   %s %s (mean_rank=%s)" % (term,ontology.get(term,"?"),mr)
                )

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
            self.rows.append(vals)

        return up, down
        
    
    # ########## Fisher Exact Test ###############

    # HYPERGEOMETRIC
    # scipy.stats.hypergeom.sf() is survival function (1-cdf), so only enriched genes will be significant
    # M = all genes
    # n = category members overall
    # N = sample size (resampling hits)
    # k = number of hits in category (intersection)

    def fisher_exact_test(self):
        logging.log("Running FET")
        import scipy.stats
        
        
        genes, hits, headers = self.read_resampling_file(
            self.inputs.resampling_file
        )  # use self.inputs.qval_col to determine hits
        if len(hits) > 1:

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

            results.sort(key=lambda x: x[-2])  # pvals
            for res in results:
                vals = res
                term = res[0]
                vals.append(pathways[term])
                intersection = list(filter(lambda x: x in associations[term], hits))
                intersection = ["%s/%s" % (x, genenames[x]) for x in intersection]
                vals.append(" ".join(intersection))
                self.rows.append(vals)

            logging.log("Finishing up FET")
            return len([i for i in qvals if i<0.05])

        else:
            logging.log("The file you passed in has no hits to run Pathway Enrichment")

    # ########## Ontologizer ###############

    # this method is restricted to GO terms

    # this implements the union method of:
    #  Grossman et al (2007). Improved Detection of overrepresentation of Gene-Ontology
    #  annotation with parent-child analysis. Bioinformatics, 23(22):3024-3031.

    def Ontologizer(self):
        import scipy.stats
        
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
        logging.warn(
            "GO terms with at least one ORF: %s" % len(go2rvs.keys())
        )  # what about between MIN and MAX?
        for go in go2rvs.keys():
            if go not in ontology:
                logging.warn("not found: %s" % go)  # also indicate which gene?

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
            self.rows.append(vals)
        return len([i for i in qvals if i<0.05])


@transit_tools.ResultsFile
class FETResultsFile:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier+"FET")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier+"FET",
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier+"FET",
                    heading=misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[
                        "Adj P Value"
                    ],
                ).Show(),
            })
        )
        
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        self.values_for_result_table.update(self.extra_data.get("parameters", {}))
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()

@transit_tools.ResultsFile
class GSEAResultsFile:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier+"GSEA")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier+"GSEA",
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier+"GSEA",
                    heading=misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[
                        "Adj P Value"
                    ],
                ).Show(),
            })
        )
        
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        self.values_for_result_table.update(self.extra_data.get("parameters", {}))
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()


@transit_tools.ResultsFile
class ONTResultsFile:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier+"ONT")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier+"ONT",
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier+"ONT",
                    heading=misc.human_readable_data(self.extra_data),
                    column_names=self.column_names,
                    rows=self.rows,
                    sort_by=[
                        "Adj P Value"
                    ],
                ).Show(),
            })
        )
        
        self.column_names, self.rows, self.extra_data, self.comments_string = tnseq_tools.read_results_file(self.path)
        self.values_for_result_table.update(self.extra_data.get("parameters", {}))
    
    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()