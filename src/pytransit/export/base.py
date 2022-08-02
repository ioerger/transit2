# __all__ = []
import sys

from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub

import traceback
import pytransit.transit_tools as transit_tools
from pytransit.components.menu import selected_export_menu_item
from pytransit.components.icon import InfoIcon
from pytransit.transit_tools import InvalidArgumentException

prefix = "[Export]"

class ExportGUI:
    def __init__(self):
        self.wxobj = None
        self.menuitem = None
        self.LABELSIZE = (100, -1)
        self.WIDGETSIZE = (100, -1)

    def defineMenuItem(self, wxobj, label):
        # TODO: write docstring

        self.wxobj = wxobj

        self.menuitem = wx.MenuItem(
            selected_export_menu_item,
            wx.ID_ANY,
            label,
            wx.EmptyString,
            wx.ITEM_NORMAL,
        )


class ExportMethod:
    """
    Basic class for analysis methods. Inherited by SingleMethod and ComparisonMethod.
    """

    def __init__(
        self,
        short_name,
        long_name,
        description,
        label,
        output,
        annotation_path,
        wxobj=None,
    ):
        self.short_name = short_name
        self.long_name = long_name
        self.description = description
        self.label = label
        self.output = output
        self.annotation_path = annotation_path

        self.WX_VERSION = WX_VERSION
        self.wxobj = wxobj

    @classmethod
    def from_gui(self, wxobj):
        # TODO: write docstring
        raise NotImplementedError

    @classmethod
    def fromargs(self, rawargs):
        # TODO: write docstring
        raise NotImplementedError

    @classmethod
    def fromconsole(self):
        # TODO: write docstring
        try:
            return self.fromargs(sys.argv[3:])
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
        # TODO: write docstring
        raise NotImplementedError

    def Run(self):
        # TODO write docstring
        raise NotImplementedError

    def finish(self):
        # TODO: write docstring
        if self.wxobj:
            wx.CallAfter(pub.sendMessage, "finish", msg=self.short_name.lower())

    def console_message(self, text):
        # TODO: write docstring
        sys.stdout.write("[%s] %s\n" % (self.short_name, text))

    def console_message_inplace(self, text):
        # TODO: write docstring
        sys.stdout.write("[%s] %s   \r" % (self.short_name, text))
        sys.stdout.flush()

    def transit_message_inplace(self, text):
        # TODO: write docstring
        self.console_message_inplace(text)

    def transit_error(self, text):
        transit_tools.log(text)
        if self.wxobj:
            transit_tools.show_error_dialog(text)

    def transit_warning(self, text):
        transit_tools.log(text)
        if self.wxobj:
            transit_tools.ShowWarning(text)


class SingleConditionMethod(ExportMethod):
    """
    Class to be inherited by analysis methods that determine essentiality in a single condition (e.g. Gumbel, Binomial, HMM).
    """

    def __init__(
        self,
        short_name,
        long_name,
        description,
        label,
        ctrldata,
        annotation_path,
        output,
        normalization=None,
        LOESS=False,
        ignoreCodon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        wxobj=None,
    ):
        ExportMethod.__init__(
            self,
            short_name,
            long_name,
            description,
            label,
            output,
            annotation_path,
            wxobj,
        )
        self.ctrldata = ctrldata
        self.normalization = normalization
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.n_terminus = n_terminus
        self.c_terminus = c_terminus


class DualConditionMethod(ExportMethod):
    """
    Class to be inherited by analysis methods that determine changes in essentiality between two conditions (e.g. Resampling, DEHMM).
    """

    def __init__(
        self,
        short_name,
        long_name,
        description,
        label,
        ctrldata,
        expdata,
        annotation_path,
        output,
        normalization,
        LOESS=False,
        ignoreCodon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        wxobj=None,
    ):
        ExportMethod.__init__(
            self,
            short_name,
            long_name,
            description,
            label,
            output,
            annotation_path,
            wxobj,
        )
        self.ctrldata = ctrldata
        self.expdata = expdata
        self.normalization = normalization
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.n_terminus = n_terminus
        self.c_terminus = c_terminus


class QuadConditionMethod(ExportMethod):
    """
    Class to be inherited by analysis methods that determine changes in essentiality between four conditions (e.g. GI).
    """

    def __init__(
        self,
        short_name,
        long_name,
        description,
        label,
        ctrldataA,
        ctrldataB,
        expdataA,
        expdataB,
        annotation_path,
        output,
        normalization,
        LOESS=False,
        ignoreCodon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        wxobj=None,
    ):
        ExportMethod.__init__(
            self,
            short_name,
            long_name,
            description,
            label,
            output,
            annotation_path,
            wxobj,
        )
        self.ctrldataA = ctrldataA
        self.ctrldataB = ctrldataB
        self.expdataA = expdataA
        self.expdataB = expdataB
        self.normalization = normalization
        self.LOESS = LOESS
        self.ignoreCodon = ignoreCodon
        self.n_terminus = n_terminus
        self.c_terminus = c_terminus


class TransitExport:
    def __init__(
        self, sn, ln, desc, lab, tn, method_class=ExportMethod, gui_class=ExportGUI
    ):
        self.short_name = sn
        self.long_name = ln
        self.description = desc
        self.label = lab
        self.transposons = tn
        self.method = method_class
        self.gui = gui_class()

    def __str__(self):
        return """Export Method:
    Short Name:  %s
    Long Name:   %s
    Description: %s
    Method:      %s
    GUI:         %s""" % (
            self.short_name,
            self.long_name,
            self.description,
            self.method,
            self.gui,
        )

    def fullname(self):
        return "[%s]  -  %s" % (self.short_name, self.long_name)

    def getInstructionsText(self):
        return ""

    def getDescriptionText(self):
        return self.description

    def getTransposonsText(self):
        if len(self.transposons) == 0:
            return "Tn attribute missing!"
        elif len(self.transposons) == 1:
            return "Intended for %s only" % self.transposons[0]
        elif len(self.transposons) == 2:
            return "Intended for %s or %s" % tuple(self.transposons)
        else:
            return (
                "Intended for "
                + ", ".join(self.transposons[:-1])
                + ", and "
                + self.transposons[-1]
            )
