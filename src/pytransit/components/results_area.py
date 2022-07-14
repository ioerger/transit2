import os

from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list
from pytransit.core_data import universal
from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename, subscribe, working_directory
import pytransit.gui_tools as gui_tools

from pytransit.components.generic.box import Column, Row
from pytransit.components.generic.text import Text
from pytransit.components.generic.button import Button
from pytransit.components.generic.table import Table

results = LazyDict(
    table=None,
    file_action_choice_element=None,
)
def create_results_area(frame):
    
    results_sizer = wx.StaticBoxSizer(
        wx.StaticBox(
            frame,
            wx.ID_ANY,
            "Results Files"
        ),
        wx.VERTICAL,
    )
    
    # 
    # Box
    # 
    if True:
        row = wx.BoxSizer(wx.HORIZONTAL)
        
        # 
        # displayButton
        # 
        if True:
            display_button = wx.Button(
                frame,
                wx.ID_ANY,
                "Display Table",
                wx.DefaultPosition,
                wx.DefaultSize,
                0,
            )
            row.Add(display_button, 0, wx.ALL, 5)
            
            @gui_tools.bind_to(display_button, wx.EVT_BUTTON)
            def display_file_func(event):
                next = results.table.GetNextSelected(-1)
                if next > -1:
                    dataset = results.table.GetItem(next, 3).GetText()
                    if dataset.verbose:
                        transit_tools.transit_message(
                            "Displaying results: %s"
                            % results.table.GetItem(next, 0).GetText()
                        )

                    try:
                        fileWindow = file_display.TransitGridFrame(frame, dataset)
                        fileWindow.Show()
                    except Exception as e:
                        transit_tools.transit_message(
                            "Error occurred displaying file: %s" % str(e)
                        )
                        traceback.print_exc()

                else:
                    if dataset.verbose:
                        transit_tools.transit_message("No results selected to display!")

            
        # 
        # fileActionButton
        # 
        if True:
            file_action_button = wx.Button(
                frame,
                wx.ID_ANY,
                "Display Graph",
                wx.DefaultPosition,
                wx.DefaultSize,
                0,
            )
            file_action_button.Hide()
            @gui_tools.bind_to(file_action_button, wx.EVT_BUTTON)
            def _(event):
                file_action_func(event)

            
            row.Add(file_action_button, 0, wx.ALL, 5)
            
            
            
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
            row.Add(add_file_button, 0, wx.ALL, 5)
            
            @subscribe("file")
            def add_file(data):
                results.table.add(data)
                # BOOKMARK: previous method:
                    # fullpath = data["path"]
                    # type     = data["type"]
                    # date     = data["date"]
                    # name = transit_tools.basename(fullpath)
                    # results.table.InsertItem(self.index_file, name)
                    # results.table.SetItem(self.index_file, 1, f"{type}")
                    # results.table.SetItem(self.index_file, 2, f"{date}")
                    # results.table.SetItem(self.index_file, 3, f"{fullpath}")
                    # self.index_file += 1
            
            @gui_tools.bind_to(add_file_button, wx.EVT_BUTTON)
            def _(event):
                try:
                    dlg = wx.FileDialog(
                        frame,
                        message="Choose a file",
                        defaultDir=working_directory,
                        defaultFile="",
                        wildcard="Results Files (*.dat)|*.dat;|\nResults Files (*.txt)|*.txt;|\nAll files (*.*)|*.*",
                        style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
                    )
                    if dlg.ShowModal() == wx.ID_OK:
                        paths = dlg.GetPaths()
                        print("You chose the following Results file(s):")
                        for fullpath in paths:
                            print("\t%s" % fullpath)
                            name = transit_tools.basename(fullpath)
                            line = open(fullpath).readline()
                            if line.startswith("#Gumbel"):
                                type = "Gumbel"
                            elif line.startswith("#Binomial"):
                                type = "Binomial"
                            elif line.startswith("#HMM - Sites"):
                                type = "HMM - Sites"
                            elif line.startswith("#HMM - Genes"):
                                type = "HMM - Genes"
                            elif line.startswith("#Resampling"):
                                type = "Resampling"
                            elif line.startswith("#DE-HMM - Sites"):
                                type = "DE-HMM - Sites"
                            elif line.startswith("#DE-HMM - Segments"):
                                type = "DE-HMM - Segments"
                            elif line.startswith("#GI"):
                                type = "GI"
                            else:
                                type = "Unknown"
                            data = {
                                "path": fullpath,
                                "type": type,
                                "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p"),
                            }
                            wx.CallAfter(pub.sendMessage, "file", data=data)
                    dlg.Destroy()
                except Exception as e:
                    transit_tools.transit_message("Error: %s" % e)
                    print("PATH", fullpath)
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno)
            
        # 
        # fileActionChoice
        # 
        if True:
            results.file_action_choice_element = wx.Choice(
                frame,
                wx.ID_ANY,
                wx.DefaultPosition,
                wx.DefaultSize,
                [
                    "[Choose Action]",  # list of available choices
                    # BOOKMARK
                ],
                0,
            )
            results.file_action_choice_element.SetSelection(0)
            row.Add(results.file_action_choice_element, 0, wx.ALL, 5)
            
            @gui_tools.bind_to(results.file_action_choice_element, wx.EVT_CHOICE)
            def _(event):
                file_action_func(event)

            
        results_sizer.Add(row, 0, 0, 5)
    
    # 
    # results.table
    # 
    with Table(
        initial_columns=[ "Name", "Type", "Date", "Full Path"],
        max_size=(-1, 200)
    ) as results.table:
        
        results_sizer.Add(
            results.table.wx_object,
            1,
            wx.ALL | wx.EXPAND,
            5,
        )
        
    return results_sizer

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
            transit_tools.transit_message(
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
        transit_tools.ShowError(MSG="Please select a results file to plot!")

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
            transit_tools.ShowError(
                MSG="Need to select a 'Resampling' results file for this type of plot."
            )

    except Exception as e:
        transit_tools.transit_message("Error occurred creating plot: %s" % str(e))
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
            transit_tools.ShowError(
                MSG="Need to select a 'Resampling' results file for this type of plot."
            )

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
