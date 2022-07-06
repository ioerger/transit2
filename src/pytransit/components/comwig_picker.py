from pytransit.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub, basename
from pytransit.core_data import universal
import pytransit.gui_tools as gui_tools

def create_comwig_picker(window):
    # 
    # component
    # 
    combined_wig_file_picker = GenBitmapTextButton(
        window.mainWindow,
        1,
        gui_tools.bit_map,
        "Add Files",
        size=wx.Size(250, -1),
    )
    combined_wig_file_picker.SetBackgroundColour(gui_tools.color.green)
    
    # 
    # callback
    # 
    @gui_tools.bind_to(combined_wig_file_picker, wx.EVT_BUTTON)
    def load_combined_wig_file_func(event): # BOOKMARK: cwig_callback
        file_dialog = wx.FileDialog(
            window,
            message="Choose a cwig file",
            defaultDir=window.workdir,
            defaultFile="",
            wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
            style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
        )
        if file_dialog.ShowModal() == wx.ID_OK:
            cwig_paths = list(file_dialog.GetPaths())
            metadata_paths = []
            for fullpath in cwig_paths:
                metadata_dialog = wx.FileDialog(
                    window,
                    message=f"\n\nPick the sample metadata\nfor {basename(fullpath)}\n\n",
                    defaultDir=window.workdir,
                    defaultFile="",
                    wildcard=u"Read Files (*.wig)|*.wig;|\nRead Files (*.txt)|*.txt;|\nRead Files (*.dat)|*.dat;|\nAll files (*.*)|*.*",
                    style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR,
                )
                if metadata_dialog.ShowModal() == wx.ID_OK:
                    metadata_path = metadata_dialog.GetPaths()[0]
                    metadata_paths.append(
                        metadata_path
                    )
                
                metadata_dialog.Destroy()
            for each_cwig_path, each_metadata_path in zip(cwig_paths, metadata_paths):
                universal.session_data.add_cwig(
                    cwig_path=each_cwig_path,
                    metadata_path=each_metadata_path,
                )
        file_dialog.Destroy()
        
        # 
        # add graphical entries for each condition
        # 
        if True:
            for each_condition in universal.session_data.conditions:
                window.wig_table.actions.add(dict(
                    name=each_condition.name,
                    disabled=each_condition.is_disabled,
                ))
            
    return combined_wig_file_picker