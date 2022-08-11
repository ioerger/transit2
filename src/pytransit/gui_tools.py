import traceback
import os
from typing import NamedTuple

from pytransit.transit_tools import wx
from pytransit.core_data import universal



bit_map = None
default_padding = 5


def rgba(*args):
    return tuple(args)

def bind_to(wxPythonObj, event):
    """
        Usage:
            @bind_to(gui_object, wx.EVENT_THING)
            def callback(event):
                pass
        
    """
    def wrapper2(function_to_attach):
        wxPythonObj.Bind(event, function_to_attach)
        def wrapper1(*args, **kwargs):
            try:
               return function_to_attach(*args, **kwargs)
            except Exception as error:
                print(error)
                handle_error(error)
        return wrapper1
    return wrapper2

def handle_error(error_obj):
    """
        Summary:
            logs error message in bottom corner of GUI
            prints it to console
            and avoids crashing the whole runtime
        Usage:
            try:
                pass
            except Exception as error:
                handle_error
    """
    traceback.print_exc()
    frame = universal.frame
    if frame and hasattr(frame, "statusBar") and hasattr(frame.statusBar, "SetStatusText"):
        frame.statusBar.SetStatusText("Error: "+str(error_obj.args))

def handle_traceback(traceback_obj):
    import traceback
    frame = universal.frame
    print(''.join(traceback.format_tb(traceback_obj)))
    if frame and hasattr(frame, "statusBar") and hasattr(frame.statusBar, "SetStatusText"):
        frame.statusBar.SetStatusText("Error: "+str(error.args))

def show_message(MSG=""):
    # TODO: Write docstring
    wx.MessageBox(MSG, "Info", wx.OK | wx.ICON_INFORMATION)

def set_status(message):
    frame = universal.frame
    if frame:
        frame.statusBar.SetStatusText(message)
        wx.Yield()

def ask_for_files(message):
    import os
    frame = universal.frame
    file_dialog = wx.FileDialog(
        frame,
        message=message,
        defaultDir=os.getcwd(),
        defaultFile="",
        style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
    )
    output = None
    if file_dialog.ShowModal() == wx.ID_OK:
        output = list(file_dialog.GetPaths())
    file_dialog.Destroy()
    return output

def ask_for_output_file_path(
        default_folder=None,
        default_file_name="",
        output_extensions=u'Common output extensions (*.txt,*.dat,*.out)|*.txt;*.dat;*.out;|\nAll files (*.*)|*.*"',
    ):
        path = None
        if not default_folder:
            default_folder = os.getcwd()

        
        file_dialog = wx.FileDialog(
            universal.frame,
            message="Save file as ...",
            defaultDir=default_folder,
            defaultFile=default_file_name,
            wildcard=output_extensions,
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
        )
        if file_dialog.ShowModal() == wx.ID_OK:
            path = file_dialog.GetPath()
        file_dialog.Destroy()
        return path

class color(NamedTuple):
    black          =  rgba(0, 0, 0)
    white          =  rgba(255, 255, 255)
    light_gray     =  rgba(199, 203, 205)
    gray           =  rgba(84, 110, 122)
    rust           =  rgba(193, 126, 112)
    orange         =  rgba(247, 140, 108)
    yellow         =  rgba(254, 195, 85)
    bananna_yellow =  rgba(221, 215, 144)
    lime           =  rgba(195, 232, 141)
    green           = rgba(187, 237, 181, 255)
    bold_green     =  rgba(78, 201, 176, 209)
    vibrant_green  =  rgba(4, 216, 149)
    dim_green      =  rgba(128, 203, 171)
    dark_slate     =  rgba(63, 132, 141)
    light_slate    =  rgba(100, 186, 197)
    light_blue     =  rgba(137, 221, 255)
    blue           =  rgba(130, 170, 255)
    electric_blue  =  rgba(0, 174, 255, 232)
    light_purple   =  rgba(199, 146, 234)
    pink           =  rgba(229, 126, 179)
    red            =  rgba(255, 85, 114)
    soft_red       =  rgba(240, 113, 120)

class NiceErrorLog(object):
    """
        Example:
            with nice_error_log:
                print("stuff")
                raise Exception("example error")
                
                # ((error magically gets logged and caught))
    """
    def __init__(*args, **kwargs):
        pass
    
    def __enter__(self):
        frame = universal.frame
        return frame
    
    def __exit__(self, _, error, traceback_obj):
        if error is not None:
            print(''.join(traceback.format_tb(traceback_obj)))
            frame = universal.frame
            if frame:
                frame.statusBar.SetStatusText("Error: "+str(error.args))

nice_error_log = NiceErrorLog()

import wx.grid
class TransitTable(wx.grid.GridTableBase):
    """
    A custom wx.Grid Table using user supplied data
    """

    def __init__(self, data, colnames):
        """data is a list of the form
        [(rowname, dictionary),
        dictionary.get(colname, None) returns the data for column
        colname
        """
        # The base class must be initialized *first*
        wx.grid.GridTableBase.__init__(self)
        self.data = data
        self.colnames = colnames
        # XXX
        # we need to store the row length and column length to
        # see if the table has changed size
        self._rows = self.GetNumberRows()
        self._cols = self.GetNumberCols()
        self.sorted_col = None
        self.sorted_dir = None

    def GetNumberCols(self):
        return len(self.colnames)

    def GetNumberRows(self):
        return len(self.data)

    def GetColLabelValue(self, col):
        return self.colnames[col]

    def GetRowLabelValue(self, row):
        return "%d" % int(self.data[row][0])

    def GetValue(self, row, col):
        return str(self.data[row][1].get(self.GetColLabelValue(col), ""))

    def GetRawValue(self, row, col):
        return self.data[row][1].get(self.GetColLabelValue(col), "")

    def SetValue(self, row, col, value):
        self.data[row][1][self.GetColLabelValue(col)] = value

    def SortColumn(self, col):
        if self.sorted_col == col:
            self.sorted_dir = not self.sorted_dir
        else:
            self.sorted_col = col
            self.sorted_dir = False

        name = self.colnames[col]
        tempdata = []

        for row in self.data:
            rowname, entry = row
            try:
                tempval = float(entry.get(name, None))
            except:
                tempval = entry.get(name, None)
            tempdata.append((tempval, row))

        tempdata.sort(reverse=self.sorted_dir)
        self.data = []

        for sortvalue, row in tempdata:
            self.data.append(row)

class SpreadSheet(wx.Frame):
    max_column_size = 200
    min_column_size = 100
    max_width = 1500
    max_height = 800
    
    def __init__(self, title, heading, column_names, rows):
        wx.Frame.__init__(self, universal.frame, size=(-1,-1))
        self.parent = universal.frame
        self.col = 0
        self.row = 0
        
        self.title        = title
        self.heading      = heading
        self.column_names = column_names
        self.rows         = list(enumerate(rows))

        self.SetTitle(title)
        if True:
            # 
            # inner box
            # 
            outer_box_sizer = wx.BoxSizer(wx.VERTICAL)
            if True:
                inner_box_sizer = wx.StaticBoxSizer(
                    wx.StaticBox(self, wx.ID_ANY, u"Information"), wx.HORIZONTAL
                )
                header_wxobj = wx.StaticText(self, wx.ID_ANY, self.heading, wx.DefaultPosition, wx.DefaultSize, 0)
                header_wxobj.Wrap(-1)
                inner_box_sizer.Add(header_wxobj, 0, wx.ALL, default_padding)
            
            # 
            # grid
            # 
            self.grid = wx.grid.Grid(self, -1)

            outer_box_sizer.Add(inner_box_sizer, 0, wx.EXPAND, default_padding)
            outer_box_sizer.Add(self.grid, 1, wx.EXPAND, default_padding)
        
        self.SetSizer(outer_box_sizer)
        self.Centre(wx.BOTH)

        self.Bind(wx.grid.EVT_GRID_LABEL_LEFT_DCLICK, self.OnLabelDoubleClicked)
        self.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.OnCellRightClicked)

        mytable = TransitTable(self.rows, self.column_names)
        self.grid.SetTable(mytable, True)

        self.grid.EnableEditing(False)
        self.grid.AdjustScrollbars()
        self.grid.SetColLabelSize(wx.grid.GRID_AUTOSIZE)

        self.grid.AutoSizeColumns()
        self.AutoResizeCols()
        self.grid.ForceRefresh()

        (width, height) = outer_box_sizer.GetMinSize()
        width = min(width + 50, self.max_width) # I'm not sure what the +50 is for
        height = min(height, self.max_height)

        self.SetMinSize((width, height))

        self.Layout()

    def AutoResizeCols(self):
        self.grid.AutoSizeColumns(False)
        for i, label in enumerate(self.column_names):
            size = self.grid.GetColSize(i)
            if size > self.max_column_size:
                self.grid.SetColSize(i, self.max_column_size)
            elif size < self.min_column_size:
                self.grid.SetColSize(i, self.min_column_size)

    def OnLabelDoubleClicked(self, evt):
        col = evt.GetCol()
        if col != -1:
            self.grid.GetTable().SortColumn(col)
            self.grid.ForceRefresh()

    def OnCellRightClicked(self, evt):

        menu = wx.Menu()
        id1 = wx.NewId()
        sortID = wx.NewId()

        xo, yo = evt.GetPosition()
        self.row = self.grid.YToRow(yo) - 1
        self.col = 0
        val = self.grid.GetCellValue(self.row, 0)

        self.Refresh()
        for (menuname, menufunc) in self.filetype.getMenus():
            newid = wx.NewId()
            menu.Append(newid, menuname)
            newmenufunc = partial(menufunc, self)
            self.Bind(wx.EVT_MENU, newmenufunc, id=newid)

        self.PopupMenu(menu)
        menu.Destroy()
