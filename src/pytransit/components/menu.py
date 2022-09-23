from collections import defaultdict
from functools import partial

from pytransit.tools import logging, gui_tools, transit_tools, tnseq_tools, norm_tools, stat_tools
import pytransit.components.qc_display as qc_display
from pytransit.universal_data import SessionData, universal
import pytransit

selected_export_menu_item = None
convert_menu_item = None
documentation_url = "http://saclab.tamu.edu/essentiality/transit/transit.html"

def create_menu(frame):
    # must imported inside the function to avoid circular import
    import pytransit.components.parameter_panel as parameter_panel
    from pytransit.methods.analysis   import methods as analysis_methods
    from pytransit.methods.export     import methods as export_methods
    from pytransit.methods.convert    import methods as convert_methods
    from pytransit.tools.transit_tools import HAS_WX, wx, GenBitmapTextButton, pub
    
    global selected_export_menu_item
    global convert_menu_item
    menu_bar = wx.MenuBar(0)
    frame = universal.frame
    
    # 
    # File Menu
    # 
    if True:
        file_menu = wx.Menu()
    
        # 
        # Export Menu
        # 
        if True:
            export_menu_item = wx.Menu()
            
            # 
            # Selected Samples
            # 
            if True:
                selected_export_menu_item = wx.Menu()
                export_menu_item.AppendSubMenu(
                    selected_export_menu_item, "Selected Samples"
                )
                
                # 
                # find export options
                # 
                def when_export_clicked(selected_name, event=None):
                    with gui_tools.nice_error_log:
                        if frame.verbose: logging.log(f"Selected Export Method: {selected_name}")
                        gui_tools.run_method_by_label(method_options=export_methods, method_label=selected_name)
                
                for name in export_methods:
                    method = export_methods[name]
                    method.gui.define_menu_item(frame, method.label)
                    temp_menu_item = method.gui.menuitem
                    selected_export_menu_item.Append(temp_menu_item)
                    
                    frame.Bind(
                        wx.EVT_MENU,
                        partial(when_export_clicked, method.label),
                        temp_menu_item,
                    )

            file_menu.AppendSubMenu(export_menu_item, "Export")

        # 
        # Convert
        # 
        if True:
            convert_menu_item = wx.Menu()
            
            # 
            # prot_table to PTT
            # 
            annotation_convert_pt_to_ptt_menu = wx.MenuItem(
                convert_menu_item,
                wx.ID_ANY,
                "prot_table to PTT",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            convert_menu_item.Append(annotation_convert_pt_to_ptt_menu)
            def when_annotation_pt_to_ptt_clicked(event):
                with gui_tools.nice_error_log:
                    annotation_path = universal.session_data.annotation_path
                    default_file = transit_tools.fetch_name(annotation_path) + ".ptt.table"
                    # default_dir = os.path.dirname(os.path.realpath(__file__))
                    default_dir = os.getcwd()

                    if not annotation_path:
                        # NOTE: was a popup
                        logging.error("Error: No annotation file selected.")
                    else:

                        output_path = frame.SaveFile(default_dir, default_file)
                        if not output_path:
                            return
                        if frame.verbose:
                            logging.log(
                                "Converting annotation file from prot_table format to PTT format"
                            )
                        from pytransit.tools.transit_tools import gather_sample_data_for
                        data, position = gather_sample_data_for(selected_samples=True)
                        orf2info = transit_tools.get_gene_info(annotation_path)
                        hash = transit_tools.get_pos_hash(annotation_path)
                        (orf2reads, orf2pos) = tnseq_tools.get_gene_reads(
                            hash, data, position, orf2info
                        )

                        output = open(output_path, "w")
                        output.write("geneID\tstart\tend\tstrand\tTA coordinates\n")
                        with open(annotation_path) as file:
                            for line in file:
                                if line.startswith("#"):
                                    continue
                                tmp = line.strip().split("\t")
                                orf = tmp[8]
                                name = tmp[7]
                                desc = tmp[0]
                                start = int(tmp[1])
                                end = int(tmp[2])
                                strand = tmp[3]
                                ta_str = "no TAs"
                                if orf in orf2pos:
                                    ta_str = "\t".join([str(int(ta)) for ta in orf2pos[orf]])
                                output.write("%s\t%s\t%s\t%s\t%s\n" % (orf, start, end, strand, ta_str))
                        output.close()
                        if frame.verbose:
                            logging.log("Finished conversion")

            frame.Bind(wx.EVT_MENU, when_annotation_pt_to_ptt_clicked, id=annotation_convert_pt_to_ptt_menu.GetId(),  )
            
            # 
            # prot_table to GFF3
            # 
            annotation_convert_pt_to_gff3_menu = wx.MenuItem(
                convert_menu_item,
                wx.ID_ANY,
                "prot_table to GFF3",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            convert_menu_item.Append(annotation_convert_pt_to_gff3_menu)
            def when_annotation_pt_to_gff3_clicked(event):
                with gui_tools.nice_error_log:
                    annotation_path = universal.session_data.annotation_path
                    default_file = transit_tools.fetch_name(annotation_path) + ".gff3"
                    # default_dir = os.path.dirname(os.path.realpath(__file__))
                    default_dir = os.getcwd()
                    output_path = frame.SaveFile(default_dir, default_file)

                    ORGANISM = transit_tools.fetch_name(annotation_path)
                    if not annotation_path:
                        # NOTE: was a popup
                        logging.error("Error: No annotation file selected.")

                    elif output_path:
                        if frame.verbose:
                            logging.log(
                                "Converting annotation file from prot_table format to GFF3 format"
                            )
                        year = time.localtime().tm_year
                        month = time.localtime().tm_mon
                        day = time.localtime().tm_mday

                        output = open(output_path, "w")
                        output.write("##gff-version 3\n")
                        output.write("##converted to IGV with TRANSIT.\n")
                        output.write("##date %d-%d-%d\n" % (year, month, day))
                        output.write("##Type DNA %s\n" % ORGANISM)

                        with open(annotation_path) as file:
                            for line in file:
                                if line.startswith("#"):
                                    continue
                                tmp = line.strip().split("\t")
                                desc = tmp[0]
                                start = int(tmp[1])
                                end = int(tmp[2])
                                strand = tmp[3]
                                length = tmp[4]
                                name = tmp[7]
                                orf = tmp[8]
                                ID = name
                                desc.replace("%", "%25").replace(";", "%3B").replace(
                                    "=", "%3D"
                                ).replace(",", "%2C")
                                output.write(
                                    "%s\tRefSeq\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s;Alias=%s;locus_tag=%s;desc=%s\n"
                                    % (ORGANISM, start, end, strand, orf, ID, orf, orf, desc)
                                )

                        output.close()
                        if frame.verbose:
                            logging.log("Finished conversion")
                            
            frame.Bind(wx.EVT_MENU, when_annotation_pt_to_gff3_clicked, id=annotation_convert_pt_to_gff3_menu.GetId(), )

            # 
            # PTT to prot_table
            # 
            annotation_convert_ptt_to_pt = wx.MenuItem(
                convert_menu_item,
                wx.ID_ANY,
                "PTT to prot_table",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            convert_menu_item.Append(annotation_convert_ptt_to_pt)
            def when_annotation_ptt_to_pt_clicked(event):
                with gui_tools.nice_error_log:
                    
                    annotation_path = universal.session_data.annotation_path
                    default_file = transit_tools.fetch_name(annotation_path) + ".prot_table"
                    # default_dir = os.path.dirname(os.path.realpath(__file__))
                    default_dir = os.getcwd()

                    if not annotation_path:
                        # NOTE: was a popup
                        logging.error("Error: No annotation file selected.")
                    else:

                        output_path = frame.SaveFile(default_dir, default_file)
                        if not output_path:
                            return
                        if frame.verbose:
                            logging.log(
                                "Converting annotation file from PTT format to prot_table format"
                            )

                        output = open(output_path, "w")
                        with open(annotation_path) as file:
                            for line in file:
                                if line.startswith("#"):
                                    continue
                                if line.startswith("geneID"):
                                    continue
                                tmp = line.strip().split("\t")
                                orf = tmp[0]
                                if orf == "intergenic":
                                    continue
                                name = "-"
                                desc = "-"
                                start = int(tmp[1])
                                end = int(tmp[2])
                                length = ((end - start + 1) / 3) - 1
                                strand = tmp[3]
                                someID = "-"
                                someID2 = "-"
                                COG = "-"
                                output.write(
                                    "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n"
                                    % (
                                        desc,
                                        start,
                                        end,
                                        strand,
                                        length,
                                        someID,
                                        someID2,
                                        name,
                                        orf,
                                        COG,
                                    )
                                )
                        output.close()
                        if frame.verbose:
                            logging.log("Finished conversion")

            frame.Bind(wx.EVT_MENU, when_annotation_ptt_to_pt_clicked , id=annotation_convert_ptt_to_pt.GetId(),       )
            
            
            # 
            # find Convert options
            # 
            def when_convert_clicked(selected_name, event=None):
                with gui_tools.nice_error_log:
                    if frame.verbose: logging.log(f"Selected Convert Method: {selected_name}")
                    gui_tools.run_method_by_label(method_options=convert_methods, method_label=selected_name)

            for name in convert_methods:
                convert_methods[name].gui.define_menu_item(frame, convert_methods[name].label)
                temp_menu_item = convert_methods[name].gui.menuitem
                convert_menu_item.Append(temp_menu_item)

                frame.Bind(
                    wx.EVT_MENU,
                    partial(when_convert_clicked, convert_methods[name].label),
                    temp_menu_item,
                )
            
            
            file_menu.AppendSubMenu(convert_menu_item, "Convert")
        
        # 
        # Exit
        # 
        if True:
            exit_option = wx.MenuItem( file_menu, wx.ID_ANY, "&Exit", wx.EmptyString, wx.ITEM_NORMAL )
            file_menu.Append(exit_option)
            def when_exit_clicked(event):
                if frame.verbose: logging.log("Exiting Transit")
                frame.Close()
            frame.Bind(wx.EVT_MENU, when_exit_clicked, id=exit_option.GetId())
        
        menu_bar.Append(file_menu, "&File")
    
    # 
    # View Menu
    # 
    if True:
        view_menu_item = wx.Menu()
        
        # 
        # Scatter Plot
        # 
        if True:
            scatter_menu_item = wx.MenuItem(
                view_menu_item,
                wx.ID_ANY,
                "&Scatter Plot",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            view_menu_item.Append(scatter_menu_item)
            def when_scatter_plot_clicked(event):
                with gui_tools.nice_error_log:
                    import numpy
                    import matplotlib
                    import matplotlib.pyplot as plt
                    from pytransit.tools import stat_tools
                    selected_samples = universal.session_data.selected_samples
                    if len(selected_samples) == 2:
                        if frame.verbose: logging.log( f"Showing scatter plot for: {[ each_sample.id for each_sample in selected_samples ]}")
                        from pytransit.tools.transit_tools import gather_sample_data_for
                        data, position = gather_sample_data_for(selected_samples=True)
                        x = data[0, :]
                        y = data[1, :]

                        plt.plot(x, y, "bo")
                        plt.title("Scatter plot - Reads at TA sites")
                        plt.xlabel(selected_samples[0].id)
                        plt.ylabel(selected_samples[1].id)
                        plt.show()
                    else:
                        # NOTE: was a popup
                        logging.error("Please make sure only two samples are selected")
            frame.Bind(wx.EVT_MENU, when_scatter_plot_clicked, id=scatter_menu_item.GetId() )

        # 
        # Track View
        # 
        if True:
            track_view_option = wx.MenuItem(view_menu_item, wx.ID_ANY, "&Track View", wx.EmptyString, wx.ITEM_NORMAL)
            view_menu_item.Append(track_view_option)
            def when_track_view_clicked(event, gene=""):
                with gui_tools.nice_error_log:
                    import pytransit.components.trash as trash
                    annotation_path = universal.session_data.annotation_path
                    wig_ids = [ each_sample.id for each_sample in universal.session_data.selected_samples ]

                    if wig_ids and annotation_path:
                        if frame.verbose:
                            logging.log(
                                "Visualizing counts for: %s"
                                % ", ".join(wig_ids)
                            )
                        view_window = trash.TrashFrame(frame, wig_ids, annotation_path, gene=gene)
                        view_window.Show()
                    elif not wig_ids:
                        # NOTE: was a popup
                        logging.error("Error: No samples selected.")
                        return
                    else:
                        # NOTE: was a popup
                        logging.error("Error: No annotation file selected.")
                        return

            frame.Bind(wx.EVT_MENU, when_track_view_clicked, id=track_view_option.GetId())
        
        # 
        # Quality Control
        # 
        if True:
            quality_control_option = wx.MenuItem(view_menu_item, wx.ID_ANY, "&Quality Control", wx.EmptyString, wx.ITEM_NORMAL )
            view_menu_item.Append( quality_control_option )
            def when_quality_control_clicked(event):
                with gui_tools.nice_error_log:
                    wig_ids = [ each_sample.id for each_sample in universal.session_data.selected_samples ] 
                    number_of_files = len(wig_ids)

                    if number_of_files <= 0:
                        raise Exception(f'''No Datasets selected, unable to run''')
                    else:
                        logging.log(f"Displaying results: {wig_ids}")
                        try:
                            qc_window = qc_display.QualityControlFrame(frame, wig_ids)
                            qc_window.Show()
                        except Exception as error:
                            raise Exception(f"Error occured displaying file: {error}")
                        
            frame.Bind(wx.EVT_MENU, when_quality_control_clicked, id=quality_control_option.GetId())
        
        menu_bar.Append(view_menu_item, "&View")
    
    # 
    # Analysis Menu
    # 
    if True:
        analysis_menu = wx.Menu()
        
        # 
        # Himar1 & Tn5
        # 
        if True:
            himar1_menu = wx.Menu()
            tn5_menu = wx.Menu()

            # 
            # generate methods
            # 
            method_names = sorted(analysis_methods.keys())
            for name in method_names:
                method = analysis_methods[name]
                if hasattr(method, "define_panel"):
                    menu_callback = method.define_panel
                    
                    # 
                    # himar1 and tn5 menu children
                    # 
                    for transposon_name, parent_menu in [ ["himar1", himar1_menu], ["tn5", tn5_menu] ]:
                        if transposon_name in method.transposons:
                            temp_menu_item = wx.MenuItem(parent_menu, wx.ID_ANY, method.full_name, wx.EmptyString, wx.ITEM_NORMAL)
                            frame.Bind(wx.EVT_MENU, menu_callback, temp_menu_item)
                            parent_menu.Append(temp_menu_item)
            
            analysis_menu.AppendSubMenu(himar1_menu, "&Himar1 Methods")
            analysis_menu.AppendSubMenu(tn5_menu, "&Tn5 Methods")
        
        menu_bar.Append(analysis_menu, "&Analysis")

    # 
    # Help Menu
    # 
    if True:
        help_menu = wx.Menu()
        
        # 
        # Documentation
        # 
        if True:
            documentation_option = wx.MenuItem(
                help_menu,
                wx.ID_ANY,
                "&Documentation",
                wx.EmptyString,
                wx.ITEM_NORMAL,
            )
            help_menu.Append(documentation_option)
            def when_documentation_clicked(event):
                with gui_tools.nice_error_log:
                    from pytransit.basics.misc import open_url
                    open_url(documentation_url)
                    
            frame.Bind(wx.EVT_MENU, when_documentation_clicked, id=documentation_option.GetId())
        
        # 
        # About
        #
        if True: 
            about_option = wx.MenuItem(
                help_menu, wx.ID_ANY, "&About", wx.EmptyString, wx.ITEM_NORMAL
            )
            help_menu.Append(about_option)
            def when_about_option_clicked(event):
                with gui_tools.nice_error_log:
                    import pytransit.components.images as images
                    
                    description = """TRANSIT is a tool for analysing TnSeq data. It provides an easy to use graphical interface and access to several different analysis methods that allow the user to determine essentiality within a single condition as well as between two conditions.


                        If you need to cite this tool, please use the following reference:

                        DeJesus, M.A., Ambadipudi, C., Baker, R., Sassetti, C., and Ioerger, T.R. (2015). TRANSIT - a Software Tool for Himar1 TnSeq Analysis. PLOS Computational Biology, 11(10):e1004401


                    """.replace("\n                    ","\n")

                    licence = """
                        TRANSIT is free software: you can redistribute it and/or modify
                        it under the terms of the GNU General Public License as published by
                        the Free Software Foundation, either version 3 of the License.


                        TRANSIT is distributed in the hope that it will be useful,
                        but WITHOUT ANY WARRANTY; without even the implied warranty of
                        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
                        GNU General Public License for more details.

                        You should have received a copy of the GNU General Public License
                        along with TRANSIT.  If not, see <http://www.gnu.org/licenses/>.
                    """.replace("\n            ", "\n")

                    info = wx.adv.AboutDialogInfo()
                    info.SetIcon(images.transit_logo2.GetIcon())
                    # images.transit_logo2.GetImage().ConvertToBitmap()
                    info.SetName("TRANSIT")
                    info.SetVersion(pytransit.__version__)
                    info.SetDescription(description)
                    info.SetCopyright("(C) 2015\n Michael A. DeJesus\nThomas R. Ioerger")
                    info.SetWebSite("http://saclab.tamu.edu/essentiality/transit/")
                    info.SetLicence(licence)
                    info.AddDeveloper("Michael A. DeJesus")
                    info.AddDeveloper("Thomas R. Ioerger")
                    info.AddDeveloper("Chaitra Ambadipudi")
                    info.AddDeveloper("Richard Baker")
                    info.AddDeveloper("Christopher Sassetti")
                    info.AddDeveloper("Eric Nelson")
                    wx.adv.AboutBox(info)
            frame.Bind(wx.EVT_MENU, when_about_option_clicked, id=about_option.GetId())
            
            menu_bar.Append(help_menu, "&Help")
    
    frame.SetMenuBar(menu_bar)


# UNUSED 
def annotation_gff3_to_pt(event):
    with gui_tools.nice_error_log:
        annotation_path = universal.session_data.annotation_path
        default_file = transit_tools.fetch_name(annotation_path) + ".prot_table"
        # default_dir = os.path.dirname(os.path.realpath(__file__))
        default_dir = os.getcwd()

        if not annotation_path:
            # NOTE: was a popup
            logging.error("Error: No annotation file selected.")
        else:
            output_path = frame.SaveFile(default_dir, default_file)
            if not output_path:
                return
            if frame.verbose:
                logging.log(
                    "Converting annotation file from GFF3 format to prot_table format"
                )

            output = open(output_path, "w")
            with open(annotation_path) as file:
                for line in file:
                    if line.startswith("#"):
                        continue
                    tmp = line.strip().split("\t")
                    chr = tmp[0]
                    type = tmp[2]
                    start = int(tmp[3])
                    end = int(tmp[4])
                    length = ((end - start + 1) / 3) - 1
                    strand = tmp[6]
                    features = dict([tuple(f.split("=")) for f in tmp[8].split(";")])
                    if "ID" not in features:
                        continue
                    orf = features["ID"]
                    name = features.get("Name", "-")
                    if name == "-":
                        name = features.get("name", "-")

                    desc = features.get("Description", "-")
                    if desc == "-":
                        desc = features.get("description", "-")
                    if desc == "-":
                        desc = features.get("Desc", "-")
                    if desc == "-":
                        desc = features.get("desc", "-")

                    someID = "-"
                    someID2 = "-"
                    COG = "-"
                    output.write(
                        "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n"
                        % (
                            desc,
                            start,
                            end,
                            strand,
                            length,
                            someID,
                            someID2,
                            name,
                            orf,
                            COG,
                        )
                    )
            output.close()
            if frame.verbose:
                logging.log("Finished conversion")