import os

from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.universal_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename, subscribe, working_directory, load_known_transit_file

import pytransit.gui_tools as gui_tools
import pytransit.file_display as file_display

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table

results = LazyDict(
    header=None,
    table=None,
    file_action_choice_element=None,
)
def create_results_area(frame):
    
    results_sizer = wx.BoxSizer(wx.VERTICAL)
    
    # 
    # Box
    # 
    if True:
        results_sizer.Add(
            wx.StaticText(
                universal.frame,
                label="Results Files",
                size=(330, -1),
            ),
            proportion=0,
            flag=wx.LEFT,
            border=gui_tools.default_padding,
        )
        
        results.header = wx.BoxSizer(wx.HORIZONTAL)
            
        # 
        # add_file_button
        # 
        if True:
            add_file_button = GenBitmapTextButton(
                frame,
                1,
                gui_tools.bit_map,
                "Add Results File",
                size=wx.Size(150, 30)
            )
            results.header.Add(add_file_button, 0, wx.ALL, gui_tools.default_padding)
            
            @gui_tools.bind_to(add_file_button, wx.EVT_BUTTON)
            def _(event):
                with gui_tools.nice_error_log:
                    dlg = wx.FileDialog(
                        frame,
                        message="Choose a file",
                        defaultDir=working_directory,
                        defaultFile="",
                        wildcard="Results Files (*.dat)|*.dat;|\nResults Files (*.txt)|*.txt;|\nAll files (*.*)|*.*",
                        style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
                    )
                    paths = []
                    if dlg.ShowModal() == wx.ID_OK:
                        paths = dlg.GetPaths()
                    dlg.Destroy()
                    # ^close the menu before processing the files to be more responsive
                    
                    print("You chose the following Results file(s):")
                    for fullpath in paths:
                        add(fullpath)
            
        # 
        # fileActionChoice
        # 
        if True:
            # by default no options are available
            change_file_action_choices({
                "[None]": lambda event: None
            })

            
        results_sizer.Add(results.header, 0, 0, gui_tools.default_padding)
    
    # 
    # results.table
    # 
    with Table(
        initial_columns=[ "name", "type", "path"],
        max_size=(-1, 200)
    ) as results.table:
    
        @results.table.events.on_select
        def _(event):
            row = results.table.rows[event.GetIndex()]
            dropdown_options_for_row = row.get("__dropdown_options", None)
            if isinstance(dropdown_options_for_row, dict):
                # attach all their callbacks to the dropdown
                change_file_action_choices(dropdown_options_for_row)
            
        results_sizer.Add(
            results.table.wx_object,
            1,
            wx.ALL | wx.EXPAND,
            gui_tools.default_padding,
        )
        
    return results_sizer


def change_file_action_choices(new_choices):
    new_choices = {
        "[None]": lambda event: None, # always have a none option, and always make it the first option
        **new_choices,
    }
    
    # hide the old one before showing the new one
    if results.file_action_choice_element != None:
        results.file_action_choice_element.Hide()
    
    results.file_action_choice_element = wx.Choice(
        universal.frame,
        wx.ID_ANY,
        wx.DefaultPosition,
        wx.DefaultSize,
        list(new_choices.keys()),
        0,
    )
    results.file_action_choice_element.SetSelection(0)
    results.header.Add(results.file_action_choice_element, proportion=1, flag=wx.ALL, border=gui_tools.default_padding)
    
    @gui_tools.bind_to(results.file_action_choice_element, wx.EVT_CHOICE)
    def _(event):
        choice = results.file_action_choice_element.GetString(results.file_action_choice_element.GetCurrentSelection())
        # run the callback that corrisponds the the choice
        new_choices[choice](event)
    

def add(path):
    # if not a recognized file type
    values_for_result_table = dict(
        name=basename(path),
        type="Unknown",
        path=path,
    )
    # if recognized
    result_object = load_known_transit_file(path)
    if result_object:
        values_for_result_table = result_object.values_for_result_table
    
    print(f'''result_object = {result_object}''')
    print(f'''values_for_result_table = {values_for_result_table}''')
    if HAS_WX:
        results.table.add(values_for_result_table)
# 
# 
# helpers
# 
# 
def file_action_func(event):
    # 0 - nothing
    # 1 - Volcano
    # 2 - Hist gene counts ratio
    plot_choice = results.file_action_choice_element.GetCurrentSelection()
    plot_name = results.file_action_choice_element.GetString(plot_choice)
    if plot_name == "[Choose Action]":
        return
    next = results.table.GetNextSelected(-1)
    if next > -1:
        dataset_path = results.table.GetItem(next, 3).GetText()
        dataset_name = results.table.GetItem(next, 0).GetText()
        dataset_type = results.table.GetItem(next, 1).GetText()

        if frame.verbose:
            transit_tools.log(
                "Performing the '%s' action on dataset '%s'"
                % (plot_name, dataset_name)
            )

        if plot_name == "Create a Volcano Plot":
            graph_volcano_plot(dataset_name, dataset_type, dataset_path)
        elif plot_name == "Plot Histogram of logFC of Gene Counts":
            graph_gene_counts(dataset_name, dataset_type, dataset_path)
        elif plot_name == "Plot Ranked Probability of Essentiality":
            graph_ranked_zbar(dataset_name, dataset_type, dataset_path)
        else:
            return

        results.file_action_choice_element.SetSelection(0)
    else:
        transit_tools.show_error_dialog("Please select a results file to plot!")

def graph_gene_counts(dataset_name, dataset_type, dataset_path):
    try:
        if dataset_type == "Resampling":
            X = []
            for line in open(dataset_path):
                if line.startswith("#"):
                    continue
                tmp = line.strip().split("\t")
                try:
                    log2FC = float(tmp[-3])
                except:
                    log2FC = 0
                X.append(log2FC)

            n, bins, patches = plt.hist(
                X, density=1, facecolor="c", alpha=0.75, bins=100
            )
            plt.xlabel("log2 FC - Total Gene Counts")
            plt.ylabel("Probability")
            plt.title(
                "Histogram of log2 Fold Change for Total Normalized Counts within Genes"
            )
            plt.axvline(0, color="r", linestyle="dashed", linewidth=3)
            plt.grid(True)
            plt.show()
        else:
            transit_tools.show_error_dialog("Need to select a 'Resampling' results file for this type of plot.")

    except Exception as e:
        transit_tools.log("Error occurred creating plot: %s" % str(e))
        traceback.print_exc()


def graph_volcano_plot(dataset_name, dataset_type, dataset_path):
    try:
        if dataset_type == "Resampling":
            X = []
            Y = []
            header = []
            qval_list = []
            bad = []
            col_logFC = -6
            col_pval = -2
            col_qval = -1
            ii = 0
            for line in open(dataset_path):
                if line.startswith("#"):
                    tmp = line.split("\t")
                    temp_col_logfc = [
                        i
                        for (i, x) in enumerate(tmp)
                        if "logfc" in x.lower()
                        or "log-fc" in x.lower()
                        or "log2fc" in x.lower()
                    ]
                    temp_col_pval = [
                        i
                        for (i, x) in enumerate(tmp)
                        if ("pval" in x.lower() or "p-val" in x.lower())
                        and "adj" not in x.lower()
                    ]
                    if temp_col_logfc:
                        col_logFC = temp_col_logfc[-1]
                    if temp_col_pval:
                        col_pval = temp_col_pval[-1]
                    continue

                tmp = line.strip().split("\t")
                try:
                    log10qval = -math.log(float(tmp[col_pval].strip()), 10)
                except ValueError as e:
                    bad.append(ii)
                    log10qval = 0

                log2FC = float(tmp[col_logFC])

                qval_list.append(
                    (float(tmp[col_qval]), float(tmp[col_pval].strip()))
                )
                X.append(log2FC)
                Y.append(log10qval)
                ii += 1
            count = 0
            threshold = 0.00001
            backup_thresh = 0.00001
            qval_list.sort()
            for (q, p) in qval_list:
                backup_thresh = p
                if q > 0.05:
                    break
                threshold = p
                count += 1

            if threshold == 0:
                threshold = backup_thresh
            for ii in bad:
                Y[ii] = max(Y)
            plt.plot(X, Y, "bo")
            plt.axhline(
                -math.log(threshold, 10), color="r", linestyle="dashed", linewidth=3
            )
            plt.xlabel("Log Fold Change (base 2)")
            plt.ylabel("-Log p-value (base 10)")
            plt.suptitle("Resampling - Volcano plot")
            plt.title("Adjusted threshold (red line): %1.8f" % threshold)
            plt.show()
        else:
            transit_tools.show_error_dialog("Need to select a 'Resampling' results file for this type of plot.")

    except Exception as e:
        print("Error occurred creating plot:", str(e))


def graph_ranked_zbar(dataset_name, dataset_type, dataset_path):
    try:
        X = []
        Y = []
        for line in open(dataset_path):
            if line.startswith("#"):
                continue
            tmp = line.strip().split("\t")
            try:
                # log2FC = math.log(float(tmp[6])/float(tmp[5]),2)
                zbar = float(tmp[-2])
            except:
                zbar = 0
            if zbar >= 0:
                Y.append(zbar)

        Y.sort()
        index = range(1, len(Y) + 1)
        plt.plot(index, Y, "bo")
        plt.xlabel("Rank of Gene According to Probability of Essentiality")
        plt.ylabel("Probability of Essentiality")
        plt.title("Ranked Probability of Essentiality")
        plt.show()

    except Exception as e:
        print("Error occurred creating plot:", str(e))
