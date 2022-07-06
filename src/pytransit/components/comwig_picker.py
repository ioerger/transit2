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
        with gui_tools.nice_error_log:
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
                print(f'''universal.session_data = {universal.session_data}''')
                for each_condition in universal.session_data.conditions:
                    # window.listConditions.InsertItem(window.condition_index, name)
                    window.listConditions.SetItem(window.condition_index, window.conditions_enum["Condition"], each_condition.name)
                    # window.listConditions.SetItem(window.condition_index, window.conditions_enum["Control"], each_condition.name)
                    # window.listConditions.SetItem(window.condition_index, window.conditions_enum["Experiment"], each_condition.name)
                    # window.listConditions.SetItem(window.condition_index, window.conditions_enum["Reference"], each_condition.name)
                    # window.listConditions.SetItem(window.condition_index, window.conditions_enum["Sample Count"], each_condition.name)
                    window.listConditions.Select(window.condition_index)
                    window.condition_index += 1
            
    return combined_wig_file_picker