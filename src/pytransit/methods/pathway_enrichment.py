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

from pytransit.specific_tools.transit_tools import wx
from pytransit.components.spreadsheet import SpreadSheet
from pytransit.components.parameter_panel import panel, progress_update, set_instructions

from pytransit.specific_tools import gui_tools, transit_tools, tnseq_tools, norm_tools, console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools import csv, misc
from pytransit.specific_tools.transit_tools import wx, basename
from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

@misc.singleton
class Method:
    name = "Pathway Enrichment"
    identifier  = name.replace(" ", "")
    cli_name    = name.replace(" ", "_").lower()
    menu_name   = f"{identifier} - Perform {name} analysis"
    description = f"""Perform {name} analysis"""
    rows = []
    
    inputs = LazyDict(
        input_file = None,
        input_type = "",
        associations_file = None,
        pathways_file = None,
        output_path= None,
        organism_pathway = None,
        pval_col = -2,
        qval_col = -1,
        lfc_col = 6,
        num_permutations = 10000,
        enrichment_exponent = 0, # for GSEA
        pseudocount = 2,
        focusLFC = "all",
        minLFC = 0,
        topk = -1,
        qvalCutoff=0.05,
    )
    
    default_significance = 0.05
    
    valid_cli_flags = [
        "--M", 
        #"--p-val-col",
        #"--q-val-col",
        "--ranking",
        "--LFC-col",
        "--p",
        "--n-perm",
        "--PC",
        "--focusLFC",
        "--minLFC",
        "--qval",
        "--topk"
    ]


    usage_string = f"""
        Usage 1: # --M FET for Fisher's Exact Test (default)
            {console_tools.subcommand_prefix} pathway_enrichment <input_file> <associations> <pathways> <output_file> --M FET [Optional Arguments]
            
            Optional Arguments:
                --PC <int>        := pseudo-counts to use in calculating p-value based on hypergeometric distribution. Default: --PC 2
                --LFC-col   <int> := column index (starting at 0) for LFC's
                --focusLFC pos|neg  :=  filter the output to focus on results with positive (pos) or negative (neg) LFCs (default: "all", no filtering)
                --minLFC <float>    :=  filter the output to include only genes that have a magnitude of LFC greater than the specified value (default: 0) (e.g. '--minLFC 1' means analyze only genes with 2-fold change or greater)
                --qval <float>      :=  filter the output to include only genes that have Qval less than to the value specified (default: 0.05)
                --topk <int>        :=  calculate enrichment among top k genes ranked by significance (Qval) regardless of cutoff (can combine with --focusLFC)
        
        Usage 2:
            # GSEA for Gene Set Enrichment Method (Subramaniam et al, 2005)
            {console_tools.subcommand_prefix} pathway_enrichment <input_file> <associations> <pathways> <output_file> --M GSEA [Optional Arguments]
            
            Optional Arguments:
                --ranking <SLPV or LFC> := SLPV is signed-log-p-value, LFC is log2-fold-change from input. Default --ranking SLPV
                --p         <float>     := exponent to use in calculating enrichment score; recommend trying 0 or 1 (as in Subramaniam et al, 2005)
                --n-perm    <int>       := number of permutations to simulate for null distribution to determine p-value. Default --n-perm 10000
                --LFC-col   <int> := column index (starting at 0) for LFC's
            
        Usage 3:
            # ONT for Ontologizer (Grossman et al, 2007)
            {console_tools.subcommand_prefix} pathway_enrichment <input_file> <associations> <pathways> <output_file> --M ONT [Optional Arguments]

            Optional Arguments:
                --LFC-col   <int>       := column index (starting at 0) for LFC's
    """.replace("\n        ", "\n")
    
    @gui.add_menu("Post-Processing", menu_name)
    def on_menu_click(event):
        Method.inputs.input_file = None
        Method.define_panel(event)
    
    def call_from_results_panel(self, results_file):
        Method.inputs.input_file = results_file
        Method.define_panel()


    def create_default_pathway_button(self,panel, sizer, *, button_label, tooltip_text="Click this button to select from TRANSIT provided files"):
        import csv
        COG_orgs = []
        with open(root_folder+"src/pytransit/data/cog-20.org.csv") as file_obj:
            reader_obj = csv.reader(file_obj)
            for row in reader_obj:
                COG_orgs.append(row[1])

        path_to_assoc_dict={"Sanger":["H37Rv"], "COG": COG_orgs, "KEGG":["H37Rv"], "GO":["H37Rv", "Smeg"]}
        assoc_to_path_dict={"H37Rv":["Sanger", "KEGG", "GO"]}
        for org in COG_orgs:
            assoc_to_path_dict[org]=["COG"]

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
                    popup_sizer = wx.GridSizer(rows=3, cols=3, hgap=2, vgap=2) 
                    win.SetSizer(popup_sizer,  wx.EXPAND|wx.ALL)

                    
                    pathway_label_text= wx.StaticText(win, wx.ALL | wx.ALIGN_CENTER, label="Select A Pathway Type : ", style=wx.EXPAND)
                    popup_sizer.Add(pathway_label_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                    pathway_type = wx.ComboBox(win, size = (250, 20))
                    pathway_type.Clear()
                    pathway_type.SetItems(list(path_to_assoc_dict.keys())+["Upload my Own Pathway file"])
                    popup_sizer.Add(pathway_type,wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                    pathway_text = wx.StaticText(win, wx.ALL | wx.ALIGN_CENTER, label="", style=wx.EXPAND, size = (250, 20))
                    popup_sizer.Add(pathway_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)


                    associations_label_text= wx.StaticText(win, wx.ALL | wx.ALIGN_CENTER, label="Select An Association : ", style=wx.EXPAND)
                    popup_sizer.Add(associations_label_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                    association_type = wx.ComboBox(win, size = (250, 20))
                    association_type.Clear()
                    association_type.SetItems(list(assoc_to_path_dict.keys())+["Upload my Own Associations file"])
                    popup_sizer.Add(association_type,wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
                    association_text = wx.StaticText(win, wx.ALL | wx.ALIGN_CENTER, label="", style=wx.EXPAND, size = (250, 20))
                    popup_sizer.Add(association_text, 0, wx.ALL | wx.ALIGN_CENTER, gui_tools.default_padding)
          

                    reset_btn = wx.Button(win, wx.ID_OK, label = "Reset Choices")
                    popup_sizer.Add(reset_btn,wx.EXPAND, gui_tools.default_padding)

                    
                    @gui_tools.bind_to(pathway_type, wx.EVT_COMBOBOX)
                    def onPathwaySelect(*args,**kwargs):
                        with gui_tools.nice_error_log:
                            selected_path = pathway_type.GetStringSelection()
                            logging.log("You selected "+selected_path + " Pathway")                               
                            if "Upload" in selected_path: #upload your own pathway file
                                pathway_text.SetLabel("Custom")
                                association_type.SetItems(["Upload my Own Associations file"])

                                pathway_file_path = gui_tools.ask_for_file(
                                    message="Select Pathways File",
                                    default_folder=None,
                                    default_file_name="",
                                    allowed_extensions='All files (*.*)|*.*',
                                )
                                Method.inputs.pathways_file = pathway_file_path
                                pathway_text.SetLabel(pathway_file_path.split("/")[-1])
                                association_type.SetItems(["Upload my Own Associations file"])

                            else:
                                pathway_text.SetLabel(selected_path)
                                association_type.SetItems(path_to_assoc_dict[selected_path]+["Upload my Own Associations file"])
                                if selected_path =="COG": 
                                    association_type.SetValue('Mycobacterium_tuberculosis_H37Rv')
                                    association_text.SetLabel('Mycobacterium_tuberculosis_H37Rv')


                    @gui_tools.bind_to(association_type, wx.EVT_COMBOBOX)
                    def onAssociationSelect(*args,**kwargs):

                        with gui_tools.nice_error_log:
                            selected_org = association_type.GetStringSelection()
                            logging.log("You selected "+selected_org + " Associations")
                            if "Upload" in selected_org:
                                # set the file path variable
                                associations_file_path = gui_tools.ask_for_file(
                                    message="Select Associations File",
                                    default_folder=None,
                                    default_file_name="",
                                    allowed_extensions='All files (*.*)|*.*',
                                )
                                Method.inputs.associations_file = associations_file_path
                                association_text.SetLabel(associations_file_path.split("/")[-1])
                                pathway_type.SetItems(["Upload my Own Associations file"])
                            else:
                                association_text.SetLabel(selected_org)
                                pathway_type.SetItems(assoc_to_path_dict[selected_org]+["Upload my Own Pathways file"])
                                    

                    @gui_tools.bind_to(reset_btn, wx.EVT_BUTTON)
                    def when_reset_button_clicked(*args,**kwargs):
                        pathway_text.SetLabel("")
                        pathway_type.Clear()
                        pathway_type.SetItems(list(path_to_assoc_dict.keys())+["Upload my Own Pathway file"])
                        association_text.SetLabel("")
                        association_type.Clear()
                        association_type.SetItems(list(assoc_to_path_dict.keys())+["Upload my Own Associations file"])                    
                    
                    ok_btn = wx.Button(win, wx.ID_OK, label = "Done")
                    popup_sizer.Add(ok_btn,wx.EXPAND, gui_tools.default_padding)
                    win.Layout()
                    popup_sizer.Fit(win)
                    res = win.ShowModal()
                    if res == wx.ID_OK :
                        organism_pathway = association_text.GetLabel() + "-" + pathway_text.GetLabel()
                        display_text= association_text.GetLabel().split("_")[-1]+ "-" + pathway_text.GetLabel()
                        organism_pathway_text.SetLabel(display_text)
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
                title_text=self.name,
                sub_text="",
                method_specific_instructions="""
                Pathway Enrichment Analysis provides a method to identify enrichment of functionally-related genes among those that are conditionally essential (i.e. significantly more or less essential between two conditions). The analysis is typically applied as post-processing step to the hits identified by a comparative analysis, such as input. Several analytical method are provided: Fisher’s exact test (FET, hypergeometric distribution), GSEA (Gene Set Enrichment Analysis) by Subramanian et al (2005), and Ontologizer. 


                    1. If you have selected this method from the menu bar, ensure you select a input file from the Select 
                       Input File button below

                    2. Add in Association/Pathway Information
                        a. Choose from one of our provided associations and pathways using the "Select from Provided Files" OR
                        b. Select your own using the Select Custom Associations and Select Custom Pathways Buttons

                    3. Select Pathway Enrichement Method. 
                        * Ontologizer is best suited for GO Terms

                    4. [Optional] Adjust parameters
                    
                    5. Click Run
                """.replace("\n                    ","\n"),
            )
            panel_helpers.create_run_button(panel, main_sizer, from_gui_function = self.from_gui)
            self.value_getters = LazyDict()

            if Method.inputs.input_file == None:
                self.value_getters.input_file = panel_helpers.create_file_input(panel, main_sizer, 
                    button_label="Select Input File", 
                    tooltip_text="The input file is the one obtained after using an comparitive analysis method in Transit that provides significance. GSEA method makes usage of adjusted P-value column", 
                    popup_title="Select File with Hits",
                    allowed_extensions='All files (*.*)|*.*'
                )

            self.value_getters.organism_pathway =  self.create_default_pathway_button(panel, main_sizer, 
                button_label="Select Pathway system", 
                tooltip_text="We have a few Association and Pathway files pre-loaded for for your use. When this button is clicked, a pop-up will appear that will allow you to select a Pathway type and organism", 
            )

   
            self.value_getters.method = panel_helpers.create_choice_input(panel, main_sizer,
                label = "Method",
                options= ["FET", "GSEA", "ONT"],
                tooltip_text = "Method to use, FET for Fisher's Exact Test (default), GSEA for Gene Set Enrichment Method (Subramaniam et al, 2005), or ONT for Ontologizer (Grossman et al, 2007)"
            )

            self.value_getters.ranking = panel_helpers.create_choice_input(panel, main_sizer,
                label = "Ranking",
                options= ["SLPV", "LFC"],
                tooltip_text="SLPV is signed-log-p-value (default); LFC is log2-fold-change from input. This parameter is only used in GSEA"
                #PV confoundes depleted AND enriched hits in the top. GSEA prefers ordered ranking by depletion (neg to pos). SPLV allows for positive LFC with significance at the top and negative LFC with signifincae is at the bottom
            )

            #self.value_getters.pval_col            = panel_helpers.create_int_getter(panel, main_sizer, label_text="P-Value Column Index", default_value=-2, tooltip_text="Column index with raw P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for sorting) (default: -2, i.e. second-to-last column)")
            #self.value_getters.qval_col            = panel_helpers.create_int_getter(panel, main_sizer, label_text="Q-Value Column Index", default_value=-1, tooltip_text="Column index with adjusted P-values (starting with 0; can also be negative, i.e. -1 means last col) (used for significant cutoff) (default: -1)")
            self.value_getters.lfc_col             = panel_helpers.create_int_getter(panel, main_sizer, label_text="LFC Column Index", default_value=6, tooltip_text="Column index with log2FC (starting with 0; can also be negative, i.e. -1 means last col) (used for ranking genes by SLPV or LFC) (default: 6)")
                
            self.value_getters.enrichment_exponent = panel_helpers.create_int_getter(  panel, main_sizer, label_text="Enrichment Exponent (p)",    default_value=0,      tooltip_text="Exponent to use in calculating enrichment score; recommend trying 0 or 1 (as in Subramaniam et al, 2005). FIXME  By anecdotal testing we found 0 is better ")
            self.value_getters.num_permutations    = panel_helpers.create_int_getter(  panel, main_sizer, label_text="Number of Permutations", default_value=10000,  tooltip_text="Number of permutations to simulate for null distribution to determine p-value")
            self.value_getters.pseudocount         = panel_helpers.create_pseudocount_input(panel, main_sizer, default_value=2, tooltip="Pseudo-counts used in calculating pathway enrichment. Useful to dampen the effects of small counts which may lead to deceptively high enrichment scores.")
            self.value_getters.focusLFC = panel_helpers.create_choice_input(panel, main_sizer,
                label = "FocusLFC",
                # default_value="all",
                options= ["all", "pos", "neg"],
                tooltip_text="Filter the output to focus on results with positive (pos) or negative (neg) LFCs (default: \"all\", no filtering)")
            self.value_getters.minLFC           = panel_helpers.create_int_getter(panel, main_sizer, label_text="MinLFC", default_value=0, tooltip_text="Filter the output to include only genes that have a magnitude of LFC greater than the specified value (default: 0) (e.g. '-minLFC 1' means analyze only genes with 2-fold change or greater)")
            self.value_getters.qvalCutoff       = panel_helpers.create_float_getter(panel, main_sizer, label_text="Qval cutoff", default_value=0.05, tooltip_text="Filter the output to include only genes that have Qval less than to the value specified (default: 0.05)")
            self.value_getters.topk             = panel_helpers.create_int_getter(panel, main_sizer, label_text="Top-k", default_value=-1, tooltip_text="Calculate enrichment among top k genes ranked by significance (Qval) regardless of cutoff (can combine with -focusLFC)")
            


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
        
        if Method.inputs.organism_pathway == "-":
        # if Method.inputs.associations_file == None or Method.inputs.pathways_file == None :
            # print(Method.inputs.organism_pathway)
            # print(Method.inputs.associations_file)
            # print(Method.inputs.pathways_file)
            logging.error("Select Association and Pathways File")
        else:
            organism,pathway = Method.inputs.organism_pathway.split("-")
            if pathway == "COG":
                try:
                    import requests
                    URL = "https://orca1.tamu.edu/essentiality/transit/COG2020/"+organism+"_COG_20_roles.associations.txt"
                    response = requests.get(URL)
                    with open(root_folder+"src/pytransit/data/"+organism+"_COG_20_roles.associations.txt", "wb") as file:
                        file.write(response.content)
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
            elif Method.inputs.organism_pathway =="H37Rv-GO": 
                if Method.inputs.method == "ONT":
                    logging.log("Loading in H37Rv Associations for GO Pathways")
                    Method.inputs.associations_file = root_folder+"src/pytransit/data/H37Rv_GO_terms.txt"
                    Method.inputs.pathways_file = root_folder+"src/pytransit/data/gene_ontology.1_2.3-11-18.obo"
                else: 
                    logging.log("Loading in H37Rv Associations for GO Pathways")
                    Method.inputs.associations_file = root_folder+"src/pytransit/data/H37Rv_GO_terms.txt"
                    Method.inputs.pathways_file = root_folder+"src/pytransit/data/GO_term_names.dat"
            elif Method.inputs.organism_pathway =="Smeg-GO":
                if Method.inputs.method == "ONT":
                    logging.log("Loading in H37Rv Associations for GO Pathways")
                    Method.inputs.associations_file = root_folder+"src/pytransit/data/smeg_GO_terms.txt"
                    Method.inputs.pathways_file = root_folder+"src/pytransit/data/gene_ontology.1_2.3-11-18.obo"
                else:
                    logging.log("Loading in Smeg Associations for GO Pathways")
                    Method.inputs.associations_file = root_folder+"src/pytransit/data/smeg_GO_terms.txt"
                    Method.inputs.pathways_file = root_folder+"src/pytransit/data/GO_term_names.dat"  
            
        Method.inputs.output_path = gui_tools.ask_for_output_file_path(
            default_file_name=f"{Method.cli_name}_output.tsv",
            output_extensions=transit_tools.result_output_extensions,
        )


        # 
        # validate input files have been selected
        # 
        assert Method.inputs.input_file != None, "Please select a input file"
        assert Method.inputs.associations_file != None, "Please select an associations file from the pre-provided options or upload your own"
        assert Method.inputs.pathways_file != None, "Please select a pathways file from the pre-provided options or upload your own"

        return Method

    @staticmethod
    @cli.add_command(cli_name)
    def from_args(args, kwargs):
        console_tools.handle_help_flag(kwargs, Method.usage_string)
        console_tools.enforce_number_of_args(args, Method.usage_string, at_least=4)
        console_tools.handle_unrecognized_flags(Method.valid_cli_flags, kwargs, Method.usage_string)

        # save the data
        Method.inputs.update(dict(
            input_file = args[0],
            associations_file = args[1],
            pathways_file = args[2],
            output_path=args[3],
            method = kwargs.get("M", "FET"),
            #pval_col = int(kwargs.get("p-val-col", Method.inputs.pval_col)),
            #qval_col = int(kwargs.get("q-val-col", Method.inputs.qval_col)),
            ranking = kwargs.get("ranking", "SLPV"),
            lfc_col = int(kwargs.get("LFC-col", Method.inputs.lfc_col)),
            enrichment_exponent = int(kwargs.get("p", "0")),
            num_permutations = int(kwargs.get("n-perm", Method.inputs.num_permutations)),
            pseudocount = int(kwargs.get("PC", "2")),

            focusLFC = kwargs.get("focusLFC", "all"),# for FET
            minLFC = float(kwargs.get("minLFC", "0")),# for FET
            topk = int(kwargs.get("topk", "-1")),# for FET
            qvalCutoff = float(kwargs.get("qval", "0.05")),# for FET
        ))
        
        Method.Run()
        
    def Run(self):
        #self.inputs.input_type = self.inputs.method
        self.rows = []
        with gui_tools.nice_error_log:
            from pytransit.components import results_area
            from pytransit.specific_tools import stat_tools
            logging.log(f"Starting {Method.identifier} analysis")
            start_time = time.time()
            
            #checking validation of inputs
            if self.inputs.method == "FET":
                self.hit_summary = {
                    "Sig Pathways":self.fisher_exact_test()
                }
                file_output_type = Method.identifier+"_FET"
                file_columns = [
                        "Pathway",
                        "Total Genes", 
                        "Number of Genes in Path",
                        "Significant Genes",
                        "Significent Genes In Path",
                        "Expected", 
                        "K Plus PC",
                        "Number Adjusted By PC",
                        "Enrichment Score" , 
                        "P Value", 
                        "Adj P Value", 
                        "Pathway Description",
                        "Relevant Genes"
                    ]
            elif self.inputs.method == "GSEA":
                up,down = self.GSEA()
                #hit summary shows # up Siginificant Pathways for Conditional Essential Genes and # down Siginificant Pathways for Conditional Non-Essential Genes
                self.hit_summary = {
                    "Sig Pathways": str(up) + " enriched;"+str(down) + " depleted",
                }
                file_output_type = Method.identifier+"_GSEA"
                file_columns = [
                        "Pathway",
                        "Pathway Description",
                        "Number of Genes in Path", 
                        "Mean Rank",
                        "Enrichment Score" , 
                        "P Value", 
                        "Adj P Value", 
                        "Relevant Genes"
                    ]
            elif self.inputs.method == "ONT":
                self.hit_summary = {
                    "Sig Pathways":self.Ontologizer()
                }
                file_output_type = Method.identifier+"_ONT"
                file_columns = [
                        "Pathway",
                        "Pathway Description",
                        "Total Genes", 
                        "Number of Genes in Path",
                        "Significant Genes",
                        "Significent Genes In Path",
                        "Expected", 
                        "K Plus PC",
                        "Number Adjusted By PC",
                        "Enrichment Score" ,  
                        "P Value", 
                        "Adj P Value",  
                        "Relevant Genes"
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
        input_column_names, input_rows, input_extra_data, input_comments_string = tnseq_tools.read_results_file(Method.inputs.input_file)
        if len(input_extra_data)==0:
            conditions = ""
        else: conditions = input_extra_data["conditions"]
        transit_tools.write_result(
            path=self.inputs.output_path, # path=None means write to STDOUT
            file_kind=file_output_type,
            rows=self.rows,
            column_names=file_columns,
            extra_info=dict(
                calculation_time=f"{(time.time() - start_time):0.1f}seconds",
                analysis_type=Method.identifier,
                files=dict(
                    input_file=Method.inputs.input_file,
                    association_path=Method.inputs.associations_file,
                    pathways_path=Method.inputs.pathways_file,
                ),
                parameters=dict(
                    conditions_tested_in_input = conditions,
                    input_type = Method.inputs.input_type,
                    method = self.inputs.method,
                    ranking = self.inputs.ranking,
                    enrichment_exponent = self.inputs.enrichment_exponent,
                    num_permutations = self.inputs.num_permutations,
                    pseudocount = self.inputs.pseudocount,
                    focusLFC = self.inputs.focusLFC,
                    minLFC = self.inputs.minLFC,
                    topk = self.inputs.topk,
                    qvalCutoff = self.inputs.qvalCutoff,
                ),
                summary_info = self.hit_summary
            ),
        )
        logging.log(f"Finished {Method.identifier} analysis in {time.time() - start_time:0.1f}sec")
        results_area.add(self.inputs.output_path)
        
        
    def read_input_file(self, filename):
        import pandas as pd
        logging.log("Reading in Input File", filename)
        lines = []
        with open(self.inputs.input_file, 'r') as file:
            lines = file.readlines()
        comments = [line for line in lines if line.startswith("#")]
        
        genes, hits, standardized_headers= [], [], []
        self.inputs.input_type = comments[0][1:]
        if "GI" in self.inputs.input_type:
            
            pval_list, qval_list=[],[]
            headers = comments[-1].split("\t")
            if len (headers)>2:
                standardized_headers = [misc.pascal_case_with_spaces(col) for col in headers]
                interactions_col = standardized_headers.index("Type Of Interaction")

            
        else:
            with open(filename) as file:
                for line in file:
                    if line[0] == "#": continue
                    headers = comments[-1].split("\t")
                    if len (headers)>2:
                        standardized_headers = [misc.pascal_case_with_spaces(col) for col in headers]
                        if "GI" in self.inputs.input_type: interactions_col = standardized_headers.index("Type Of Interaction")
                        else: 
                            self.inputs.pval_col = standardized_headers.index("P Value")
                            self.inputs.qval_col = standardized_headers.index("Adj P Value")
                    w = line.rstrip().split("\t")
                    if len(w)<len(headers): continue
                    genes.append(w)
                    if "GI" in self.inputs.input_type:
                        interaction = float(w[interactions_col])
                        if "No Interaction" in interaction:
                           hits.append(w[0]) 
                    else:
                        qval = float(w[self.inputs.qval_col])
                        if qval < Method.default_significance:
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
                    continue  # skip genes in association file that are not relevant (i.e. not in input file)
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

    # during initialization, self.inputs.input_file etc have been set, and self.inputs.output_path has been opened
    
    def GSEA(self):  
        logging.log("Running GSEA")
        data, hits, headers = self.read_input_file(
            self.inputs.input_file
        )  # hits are not used in GSEA()
        if "GI" in self.inputs.input_type : 
            logging.error("GSEA cannot be run with a GI output")
        orfs_in_input_file = [w[0] for w in data]
        headers = headers[-1].rstrip().split("\t")  # last line prefixed by '#'
        associations = self.read_associations(
            self.inputs.associations_file, filter=orfs_in_input_file
        )  # bidirectional map; includes term->genelist and gene->termlist
        # filter: project associations (of orfs to pathways) onto only those orfs appearing in the input file

        ontology = self.read_pathways(self.inputs.pathways_file)
        gene_names = {}
        for gene in data:
            gene_names[gene[0]] = gene[1]
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
            if qval < Method.default_significance:
                if mr < n2:
                    up += 1
                else:
                    down += 1

        for term, mr, es, pval, qval in results:
            rvs = terms2orfs[term]
            rvinfo = [(x, gene_names.get(x, "?"), orfs2rank.get(x, n2)) for x in rvs]
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
    # N = sample size (input hits)
    # k = number of hits in category (intersection)

    def fisher_exact_test(self):
        logging.log("Running FET")
        import scipy.stats
        
        
        genes, hits, headers = self.read_input_file(
            self.inputs.input_file
        )  # use self.inputs.qval_col to determine hits
        if len(hits) > 1:

            associations = self.read_associations(self.inputs.associations_file)
            pathways = self.read_pathways(self.inputs.pathways_file)

            # how many genes are there, and how many have associations?
            # how many genes in associations are not listed in input file?
            # do all associations have a definition in pathways?
            # how many pathways have >1 gene? (out of total?) what is max?

            # Hiearchy of flags:     
            #  
            #  START: qval (OPTIONAL) / topk [Mutually Exclusive | topk removes qval] -> focusLFC -> minLFC [END]
            #
            focus_genes = genes

            # Filter by only returning the top k genes (by q-value)
            if self.inputs.topk != -1:
            # should this account if there are mulitple genes with qval == 0 ??

                k_list = [(w[0], w[-1]) for w in focus_genes] # get a list of tuples, where it's just (orf, q-value)
                k_list = sorted(k_list, key=lambda tup: tup[1]) # sort
                k_list = k_list[:self.inputs.topk] # get top k genes

                k_list = [k[0] for k in k_list] # remove the q-values, getting just the orfs

                focus_genes = list(filter(lambda w: w[0] in k_list, focus_genes)) # then get all data points that are in our top-k subset
                hits = list(set([w[0] for w in focus_genes]) & set(hits)) 

            # Q-value filtering
            if self.inputs.qvalCutoff != 0.05 and self.inputs.topk == -1:# don't run the qvalCutoff filter if it's the default value and if topk is default (not being used)
                focus_genes = list(filter(lambda w: float(w[self.inputs.qval_col]) <= self.inputs.qvalCutoff, focus_genes))
                hits = list(set([w[0] for w in focus_genes]) & set(hits)) 

            # Sign-based log-fold-change filtering
            if self.inputs.focusLFC == "pos":
                focus_genes = list(filter(lambda w: float(w[self.inputs.lfc_col]) > 0, focus_genes))
                hits = list(set([w[0] for w in focus_genes]) & set(hits)) # filter the hits to only include positive LFCs by doing an intersection between the newly filtered orfs and the hits (that include all LFCs)
                                                                # by turning both lists into sets and intersecting (&) them, seemed to be the fastest way without adding too much more to MEM space
            elif self.inputs.focusLFC == "neg":
                focus_genes = list(filter(lambda w: float(w[self.inputs.lfc_col]) < 0, focus_genes))
                hits = list(set([w[0] for w in focus_genes]) & set(hits))

            # Minimum log-fold change filtering
            if self.inputs.minLFC != 0: # don't run the minLFC filter if it's the default value
                focus_genes = list(filter(lambda w: abs(float(w[self.inputs.lfc_col])) >= self.inputs.minLFC, focus_genes)) # we only want to keep values that are greater or equal to than the flag value
                                                                                                        # This is done intentionally after the focusLFC filter and uses its results
                hits = list(set([w[0] for w in focus_genes]) & set(hits))


            genes_with_associations = 0
            for gene in focus_genes:
                orf = gene[0]
                if orf in associations:
                    genes_with_associations += 1
            # self.rows.append("# method=FET, PC=%s" % self.inputs.pseudocount)
            # self.rows.append(
            #     "# genes with associations=%s out of %s total"
            #     % (genes_with_associations, len(genes))
            # )
            # self.rows.append("# significant genes (qval<Method.default_significance): %s" % (len(hits)))

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
                N = len(hits)  # number of input hits
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

            gene_names = {}
            for gene in genes:
                gene_names[gene[0]] = gene[1]

            results.sort(key=lambda x: x[-2])  # pvals
            for res in results:
                vals = res
                term = res[0]
                vals.append(pathways[term])
                intersection = list(filter(lambda x: x in associations[term], hits))
                intersection = ["%s/%s" % (x, gene_names[x]) for x in intersection]
                vals.append(" ".join(intersection))
                self.rows.append(vals)

            logging.log("Finishing up FET")
            return len([i for i in qvals if i<Method.default_significance])

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

        genes, pvals = {}, []
        all_orfs, studyset = [], []
        comments, headers, rows = tnseq_tools.read_result(self.inputs.input_file)
        self.inputs.input_type = comments[0]
        if "GI" in self.inputs.input_type:
            logging.error("ONT cannot be run with a GI output")
        try:
            self.inputs.qval_col = headers.index("Adj P Value")
            self.inputs.pval_col = headers.index("P Value")
        except Exception as error:
            logging.error(f'''missing header 'Adj P Value' or 'P Value' in file: {self.inputs.input_file}''')
        
        for each_row in rows:
            genes[each_row["ORF"]] = each_row
            all_orfs.append(each_row["ORF"])
                
            if each_row["Adj P Value"] < Method.default_significance:
                studyset.append(each_row["ORF"])
            pvals.append((each_row["ORF"], each_row["P Value"]))
        
        pvals.sort(key=lambda x: x[1])
        ranks = {}
        for i, (rv, pval) in enumerate(pvals):
            ranks[rv] = i + 1
        #self.rows.append("# number of input hits (qval<Method.default_significance): %s" % len(studyset))
        counts = []
        n, a = len(all_orfs), len(studyset)
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
        return len([i for i in qvals if i<Method.default_significance])


@transit_tools.ResultsFile
class FETResultsFile:
    @staticmethod
    def can_load(path):
        return transit_tools.file_starts_with(path, '#'+Method.identifier+"_FET")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier+"_FET",
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier+"_FET",
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
        #self.values_for_result_table.update(self.extra_data.get("parameters", {}))
        summary = self.extra_data.get("summary_info", {})
        summary_str = [str(summary[key])+" "+str(key) for key in sorted(summary.keys())] 
        self.values_for_result_table.update({"summary": "; ".join(summary_str) })
    

        parameters = self.extra_data.get("parameters",{})
        parameters_str = [str(key)+" : "+str(parameters[key]) for key in ["method", "conditions_tested_in_input", "input_type"]]
        self.values_for_result_table.update({"parameters": "; ".join(parameters_str) })

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
        return transit_tools.file_starts_with(path, '#'+Method.identifier+"_GSEA")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier+"_GSEA",
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier+"_GSEA",
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
        #self.values_for_result_table.update(self.extra_data.get("parameters", {}))
        summary = self.extra_data.get("summary_info", {})
        summary_str = [str(summary[key])+" "+str(key) for key in sorted(summary.keys())] 
        self.values_for_result_table.update({"summary": "; ".join(summary_str) })

        parameters = self.extra_data.get("parameters",{})
        parameters_str = [str(key)+" : "+str(parameters[key]) for key in ["method", "ranking","conditions_tested_in_input", "input_type"]]
        self.values_for_result_table.update({"parameters": "; ".join(parameters_str) })
    
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
        return transit_tools.file_starts_with(path, '#'+Method.identifier+"_ONT")
    
    def __init__(self, path=None):
        self.wxobj = None
        self.path  = path
        self.values_for_result_table = LazyDict(
            name=basename(self.path),
            type=Method.identifier+"_ONT",
            path=self.path,
            # anything with __ is not shown in the table
            __dropdown_options=LazyDict({
                "Display Table": lambda *args: SpreadSheet(
                    title=Method.identifier+"_ONT",
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
    
        summary = self.extra_data.get("summary_info", {})
        summary_str = [str(summary[key])+" "+str(key) for key in sorted(summary.keys())] 
        self.values_for_result_table.update({"summary": "; ".join(summary_str) })

        parameters = self.extra_data.get("parameters",{})
        parameters_str = [str(key)+" : "+str(parameters[key]) for key in ["method", "conditions_tested_in_input", "input_type"]]
        self.values_for_result_table.update({"parameters": "; ".join(parameters_str) })

    def __str__(self):
        return f"""
            File for {Method.identifier}
                path: {self.path}
                column_names: {self.column_names}
        """.replace('\n            ','\n').strip()