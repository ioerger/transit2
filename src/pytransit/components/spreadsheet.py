from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled

SpreadSheet = None
TransitTable = None
if gui.is_active:

    from pytransit.specific_tools import  gui_tools
    from pytransit.specific_tools.transit_tools import wx
    
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
        
        def __init__(self, title, heading, column_names, rows, sort_by=[]):
            """
            Arguments:
                title: string
                heading: string
                column_names: list of strings
                rows: list of dictionaries, keys=column names
                sort_by: list of strings
            """
            
            wx.Frame.__init__(self, gui.frame, size=(-1,-1))
            self.parent = gui.frame
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
                        wx.StaticBox(self, wx.ID_ANY, "Information"), wx.HORIZONTAL
                    )
                    header_wxobj = wx.StaticText(self, wx.ID_ANY, self.heading, wx.DefaultPosition, wx.DefaultSize, 0)
                    header_wxobj.Wrap(-1)
                    inner_box_sizer.Add(header_wxobj, 0, wx.ALL, gui_tools.default_padding)
                
                # 
                # grid
                # 
                self.grid = wx.grid.Grid(self, -1)

                outer_box_sizer.Add(inner_box_sizer, 0, wx.EXPAND, gui_tools.default_padding)
                outer_box_sizer.Add(self.grid, 1, wx.EXPAND, gui_tools.default_padding)
            
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
            
            for each_name in reversed(sort_by):
                if each_name not in self.column_names:
                    print(f"Warning: sort_by included this: {each_name}, but it wasnt a column_name: {self.column_names}")
                    continue
                else:
                    self.grid.GetTable().SortColumn(self.column_names.index(each_name))
            self.grid.ForceRefresh()
                

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
            print(f'''col = {col}''')
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
            for (menuname, menufunc) in self.filetype.get_menus():
                newid = wx.NewId()
                menu.Append(newid, menuname)
                newmenufunc = partial(menufunc, self)
                self.Bind(wx.EVT_MENU, newmenufunc, id=newid)

            self.PopupMenu(menu)
            menu.Destroy()
