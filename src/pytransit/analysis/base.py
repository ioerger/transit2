# __all__ = []
import sys

from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, WX_VERSION

import traceback
import datetime
import numpy
import pytransit.transit_tools as transit_tools
from pytransit.components.icon import InfoIcon
from pytransit.transit_tools import InvalidArgumentException

file_prefix = "[FileDisplay]"

class TransitFile:
    def __init__(self, identifier="#Unknown", colnames=[]):
        self.wxobj      = None
        self.short_name = "TRANSIT"
        self.long_name  = "TRANSIT"
        self.short_desc = "TRANSIT - Short Description"
        self.long_desc  = "TRANSIT - Long Description"
        self.identifier = identifier
        self.colnames   = colnames

    def console_message(self, text):
        
        sys.stdout.write("[%s] %s\n" % (self.short_name, text))

    def console_message_inplace(self, text):
        
        sys.stdout.write("[%s] %s   \r" % (self.short_name, text))
        sys.stdout.flush()

    def transit_message_inplace(self, text):
        
        self.console_message_inplace(text)

    def transit_error(self, text):
        transit_tools.log(text)
        if self.wxobj:
            transit_tools.show_error_dialog(text)

    def transit_warning(self, text):
        transit_tools.log(text)
        if self.wxobj:
            transit_tools.ShowWarning(text)

    def getData(self, path, colnames):
        # TODO write docstring
        row = 0
        data = []
        shownError = False
        for line in open(path):
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            tmp[-1] = tmp[-1].strip()
            # print(colnames)
            # print( len(colnames), len(tmp))
            try:
                rowdict = dict([(colnames[i], tmp[i]) for i in range(len(colnames))])
            except Exception as e:
                if not shownError:
                    self.transit_warning(
                        "Error reading data! This may be caused by trying to load a old results file, when the format has changed."
                    )
                    shownError = True
                rowdict = dict(
                    [(colnames[i], tmp[i]) for i in range(min(len(colnames), len(tmp)))]
                )
            data.append((row, rowdict))
            row += 1
        return data

    def getHeader(self, path):
        # TODO write docstring
        return "Generic Transit File Type."

    def getMenus(self):
        menus = [("Display in Track View", self.displayInTrackView)]
        return menus

    def displayInTrackView(self, displayFrame, event):

        # print("Self:", self)
        # print("Frame:", displayFrame)
        # print("Event:", event)
        # print("Frame parent:", displayFrame.parent)
        try:
            gene = displayFrame.grid.GetCellValue(displayFrame.row, 0)
            displayFrame.parent.allViewFunc(displayFrame, gene)
        except Exception as e:
            print(file_prefix, "Error occurred: %s" % e)


class AnalysisGUI:
    def __init__(self):
        self.wxobj = None
        self.panel = None
        self.LABELSIZE = (100, -1)
        self.WIDGETSIZE = (100, -1)

    def Hide(self):
        if self.panel: self.panel.Hide()

    def Show(self):
        if self.panel: self.panel.Show()

    def Enable(self):
        if self.panel: self.panel.Enable()

    def defineTextBox(
        self,
        panel,
        labelText="",
        widgetText="",
        tooltipText="",
        labSize=None,
        widgetSize=None,
    ):
        if not labSize:
            labSize = self.LABELSIZE
        if not widgetSize:
            widgetSize = self.WIDGETSIZE

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(
            panel, wx.ID_ANY, labelText, wx.DefaultPosition, labSize, 0
        )
        label.Wrap(-1)
        textBox = wx.TextCtrl(
            panel, wx.ID_ANY, widgetText, wx.DefaultPosition, widgetSize, 0
        )
        sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(textBox, 0, wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(
            InfoIcon(panel, wx.ID_ANY, tooltip=tooltipText),
            0,
            wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        return (label, textBox, sizer)

    def defineChoiceBox(
        self,
        panel,
        labelText="",
        widgetChoice=[""],
        tooltipText="",
        labSize=None,
        widgetSize=None,
    ):
        if not labSize:
            labSize = self.LABELSIZE
        if not widgetSize:
            widgetSize = self.WIDGETSIZE

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(
            panel, wx.ID_ANY, labelText, wx.DefaultPosition, labSize, 0
        )
        label.Wrap(-1)
        choiceBox = wx.Choice(
            panel, wx.ID_ANY, wx.DefaultPosition, widgetSize, widgetChoice, 0
        )
        choiceBox.SetSelection(0)
        sizer.Add(label, 0, wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(choiceBox, 0, wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(
            InfoIcon(panel, wx.ID_ANY, tooltip=tooltipText),
            0,
            wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        return (label, choiceBox, sizer)

    def defineCheckBox(
        self, panel, labelText="", widgetCheck=False, tooltipText="", widgetSize=None
    ):
        if not widgetSize:
            widgetSize = self.WIDGETSIZE
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        checkBox = wx.CheckBox(panel, label=labelText, size=widgetSize)
        checkBox.SetValue(widgetCheck)
        sizer.Add(checkBox, 0, wx.ALIGN_CENTER_VERTICAL, 5)
        sizer.Add(
            InfoIcon(panel, wx.ID_ANY, tooltip=tooltipText),
            0,
            wx.ALIGN_CENTER_VERTICAL,
            5,
        )
        return (checkBox, sizer)


class AnalysisMethod:
    """
    Basic class for analysis methods. Inherited by SingleMethod and ComparisonMethod.
    """

    def __init__(
        self,
        short_name,
        long_name,
        short_desc,
        long_desc,
        output,
        annotation_path,
        wxobj=None,
    ):
        self.short_name = short_name
        self.long_name = long_name
        self.short_desc = short_desc
        self.long_desc = long_desc
        self.output = output
        self.annotation_path = annotation_path

        self.WX_VERSION = WX_VERSION
        self.wxobj = wxobj

    @classmethod
    def from_gui(self, wxobj):
        
        raise NotImplementedError

    @classmethod
    def fromargs(self, rawargs):
        
        raise NotImplementedError

    @classmethod
    def fromconsole(self):
        
        try:
            return self.fromargs(sys.argv[2:])
        except InvalidArgumentException as e:
            print("Error: %s" % str(e))
            print(self.usage_string())
        except IndexError as e:
            print("Error: %s" % str(e))
            print(self.usage_string())
        except TypeError as e:
            print("Error: %s" % str(e))
            traceback.print_exc()
            print(self.usage_string())
        except ValueError as e:
            print("Error: %s" % str(e))
            traceback.print_exc()
            print(self.usage_string())
        except Exception as e:
            print("Error: %s" % str(e))
            traceback.print_exc()
            print(self.usage_string())
        sys.exit()

    @classmethod
    def usage_string(self):
        
        raise NotImplementedError

    def Run(self):
        # TODO write docstring
        raise NotImplementedError

    def print_members(self):
        
        members = sorted(
            [
                attr
                for attr in dir(self)
                if not callable(getattr(self, attr)) and not attr.startswith("__")
            ]
        )
        for m in members:
            print("%s = %s" % (m, getattr(self, m)))

    def add_file(self, path=None, filetype=None):

        
        if not path:
            path = self.output.name
        if not filetype:
            filetype = self.short_name

        data = {
            "path": path,
            "type": filetype,
            "date": datetime.datetime.today().strftime("%B %d, %Y %I:%M%p"),
        }

        if self.wxobj:
            wx.CallAfter(pub.sendMessage, "file", data=data)

    def finish(self):
        
        if self.wxobj:
            wx.CallAfter(pub.sendMessage, "finish", msg=self.short_name.lower())

    def console_message(self, text):
        
        sys.stdout.write("[%s] %s\n" % (self.short_name, text))

    def console_message_inplace(self, text):
        
        sys.stdout.write("[%s] %s   \r" % (self.short_name, text))
        sys.stdout.flush()

    def transit_message_inplace(self, text):
        
        self.console_message_inplace(text)

    def transit_error(self, text):
        transit_tools.log(text)
        if self.wxobj:
            transit_tools.show_error_dialog(text)

    def transit_warning(self, text):
        transit_tools.log(text)
        if self.wxobj:
            transit_tools.ShowWarning(text)


class SingleConditionMethod(AnalysisMethod):
    """
    Class to be inherited by analysis methods that determine essentiality in a single condition (e.g. Gumbel, Binomial, HMM).
    """

    def __init__(
        self,
        short_name,
        long_name,
        short_desc,
        long_desc,
        ctrldata,
        annotation_path,
        output,
        replicates="Sum",
        normalization=None,
        LOESS=False,
        ignoreCodon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        wxobj=None,
    ):
        AnalysisMethod.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            output,
            annotation_path,
            wxobj,
        )
        self.ctrldata = ctrldata
        self.replicates = replicates
        self.normalization = normalization
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.n_terminus = n_terminus
        self.c_terminus = c_terminus


class DualConditionMethod(AnalysisMethod):
    """
    Class to be inherited by analysis methods that determine changes in essentiality between two conditions (e.g. Resampling, DEHMM).
    """

    def __init__(
        self,
        short_name,
        long_name,
        short_desc,
        long_desc,
        ctrldata,
        expdata,
        annotation_path,
        output,
        normalization,
        replicates="Sum",
        LOESS=False,
        ignoreCodon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        wxobj=None,
    ):
        AnalysisMethod.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            output,
            annotation_path,
            wxobj,
        )
        self.ctrldata = ctrldata
        self.expdata = expdata
        self.normalization = normalization
        self.replicates = replicates
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.n_terminus = n_terminus
        self.c_terminus = c_terminus


class QuadConditionMethod(AnalysisMethod):
    """
    Class to be inherited by analysis methods that determine changes in essentiality between four conditions (e.g. GI).
    """

    def __init__(
        self,
        short_name,
        long_name,
        short_desc,
        long_desc,
        ctrldataA,
        ctrldataB,
        expdataA,
        expdataB,
        annotation_path,
        output,
        normalization,
        replicates="Sum",
        LOESS=False,
        ignoreCodon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        wxobj=None,
    ):
        AnalysisMethod.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            output,
            annotation_path,
            wxobj,
        )
        self.ctrldataA = ctrldataA
        self.ctrldataB = ctrldataB
        self.expdataA = expdataA
        self.expdataB = expdataB
        self.normalization = normalization
        self.replicates = replicates
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.n_terminus = n_terminus
        self.c_terminus = c_terminus


class MultiConditionMethod(AnalysisMethod):
    """
    Class to be inherited by analysis methods that compare essentiality between multiple conditions (e.g Anova).
    """

    def __init__(
        self,
        short_name,
        long_name,
        short_desc,
        long_desc,
        combined_wig,
        metadata,
        annotation_path,
        output,
        normalization=None,
        LOESS=False,
        ignoreCodon=True,
        wxobj=None,
        excluded_conditions=None,
        included_conditions=None,
        n_terminus=0.0,
        c_terminus=0.0,
    ):
        AnalysisMethod.__init__(
            self,
            short_name,
            long_name,
            short_desc,
            long_desc,
            output,
            annotation_path,
            wxobj,
        )
        self.combined_wig        = combined_wig
        self.metadata            = metadata
        self.normalization       = normalization
        self.n_terminus           = n_terminus
        self.c_terminus           = c_terminus
        self.unknown_cond_flag   = "FLAG-UNMAPPED-CONDITION-IN-WIG"
        self.excluded_conditions = excluded_conditions or []
        self.included_conditions = included_conditions or []

    def filter_wigs_by_conditions(
        self,
        data,
        conditions,
        covariates=[],
        interactions=[],
        excluded_conditions=[],
        included_conditions=[],
    ):
        """
            Filters conditions that are excluded/included.
            ([[Wigdata]], [Condition], [[Covar]], [Condition], [Condition]) -> Tuple([[Wigdata]], [Condition])
        """
        excluded_conditions, included_conditions = (
            set(excluded_conditions),
            set(included_conditions),
        )
        d_filtered, cond_filtered, filtered_indexes = [], [], []

        if len(excluded_conditions) > 0 and len(included_conditions) > 0:
            self.transit_error("Both excluded and included conditions have len > 0")
            sys.exit(0)
        elif len(excluded_conditions) > 0:
            transit_tools.log("conditions excluded: {0}".format(excluded_conditions))
            for i, c in enumerate(conditions):
                if (c != self.unknown_cond_flag) and (c not in excluded_conditions):
                    d_filtered.append(data[i])
                    cond_filtered.append(conditions[i])
                    filtered_indexes.append(i)
        elif len(included_conditions) > 0:
            transit_tools.log("conditions included: {0}".format(included_conditions))
            for i, c in enumerate(conditions):
                if (c != self.unknown_cond_flag) and (c in included_conditions):
                    d_filtered.append(data[i])
                    cond_filtered.append(conditions[i])
                    filtered_indexes.append(i)
        else:
            for i, c in enumerate(conditions):
                if c != self.unknown_cond_flag:
                    d_filtered.append(data[i])
                    cond_filtered.append(conditions[i])
                    filtered_indexes.append(i)

        covariates_filtered = [[c[i] for i in filtered_indexes] for c in covariates]
        interactions_filtered = [[c[i] for i in filtered_indexes] for c in interactions]

        return (
            numpy.array(d_filtered),
            numpy.array(cond_filtered),
            numpy.array(covariates_filtered),
            numpy.array(interactions_filtered),
        )

    # input: conditions are per wig; orderingMetdata comes from tnseq_tools.read_samples_metadata()
    # output: conditionsList is selected subset of conditions (unique, in preferred order)

    def select_conditions(
        self, conditions, included_conditions, excluded_conditions, orderingMetadata
    ):
        if len(included_conditions) > 0:
            conditionsList = included_conditions
        else:
            conditionsList = []
            for c in orderingMetadata[
                "condition"
            ]:  # the order conds appear in metadata file, duplicated for each sample
                if c not in conditionsList:
                    conditionsList.append(c)
        for c in excluded_conditions:
            if c in conditionsList:
                conditionsList.remove(c)
        return conditionsList

    def filter_wigs_by_conditions2(
        self, data, conditions, conditionsList, covariates=[], interactions=[]
    ):
        """
            Filters conditions that are excluded/included.
            ([[Wigdata]], [Condition], [[Covar]], [Condition], [Condition]) -> Tuple([[Wigdata]], [Condition])
        """
        d_filtered, cond_filtered, filtered_indexes = [], [], []

        for i, c in enumerate(conditions):
            if (c != self.unknown_cond_flag) and (c in conditionsList):
                d_filtered.append(data[i])
                cond_filtered.append(conditions[i])
                filtered_indexes.append(i)

        covariates_filtered = [[c[i] for i in filtered_indexes] for c in covariates]
        interactions_filtered = [[c[i] for i in filtered_indexes] for c in interactions]

        return (
            numpy.array(d_filtered),
            numpy.array(cond_filtered),
            numpy.array(covariates_filtered),
            numpy.array(interactions_filtered),
        )

    def filter_wigs_by_conditions3(
        self,
        data,
        fileNames,
        conditionNames,
        included_cond,
        excluded_cond,
        conditions,
        covariates=[],
        interactions=[],
    ):
        """
            Filters conditions that are excluded/included; also extract cond, covar, and interaction labels
            conditionNames: based on original Conditions column in metadata
            conditions: user might have specified an alternative column to analyze (list of labels parallel to wigs)
        """
        (
            fileNames_filtered,
            condNames_filtered,
            d_filtered,
            cond_filtered,
            filtered_indexes,
        ) = ([], [], [], [], [])

        for i in range(len(data)):
            if (
                len(included_cond) == 0 or conditionNames[i] in included_cond
            ) and conditionNames[i] not in excluded_cond:
                d_filtered.append(data[i])
                fileNames_filtered.append(fileNames[i])
                condNames_filtered.append(conditionNames[i])
                cond_filtered.append(conditions[i])
                filtered_indexes.append(i)

        covariates_filtered = [[c[i] for i in filtered_indexes] for c in covariates]
        interactions_filtered = [[c[i] for i in filtered_indexes] for c in interactions]

        return (
            numpy.array(d_filtered),
            numpy.array(fileNames_filtered),
            numpy.array(condNames_filtered),
            numpy.array(cond_filtered),
            numpy.array(covariates_filtered),
            numpy.array(interactions_filtered),
        )

    # return a hash table of parallel lists, indexed by column header
    def get_samples_metadata(self):
        data = {}
        header = None
        for line in open(self.metadata):
            if line[0] == "#":
                continue
            w = line.rstrip().split("\t")
            if header == None:
                header = w
                for col in header:
                    data[col] = []
            else:
                for i in range(len(header)):
                    data[header[i]].append(w[i])
        return data


class TransitAnalysis:
    def __init__(
        self,
        sn,
        ln,
        short_desc,
        long_desc,
        tn,
        method_class=AnalysisMethod,
        gui_class=AnalysisGUI,
        filetypes=[TransitFile],
    ):
        self.short_name  = sn
        self.long_name   = ln
        self.short_desc  = short_desc
        self.long_desc   = long_desc
        self.full_name   = f"[{self.short_name}]  -  {self.short_desc}"
        self.transposons = tn
        self.method      = method_class
        self.gui         = gui_class()
        self.filetypes   = filetypes
        self.transposons_text = transit_tools.get_transposons_text(self.transposons)

    def __str__(self):
        return f"""
            Analysis Method:
                Short Name:  {self.short_name}
                Long Name:   {self.long_name}
                Short Desc:  {self.short_desc}
                Long Desc:   {self.long_desc}
                Method:      {self.method}
                GUI:         {self.gui}
        """.replace('\n            ','\n').strip()