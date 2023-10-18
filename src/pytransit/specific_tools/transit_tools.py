# Copyright 2015.
#   Michael A. DeJesus, Chaitra Ambadipudi, and  Thomas R. Ioerger.
#
#
#    This file is part of TRANSIT.
#
#    TRANSIT is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License.
#
#
#    TRANSIT is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with TRANSIT.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import math
import warnings
import ntpath
import heapq

import numpy
import matplotlib.pyplot as plt
from pytransit.globals import logging


SEPARATOR = "," # TODO remove this
EOL = "\n" # TODO remove this

# 
# optional import: wx
# 
try:
    import wx
    import wx.xrc
    import wx.adv
    import wx.lib.scrolledpanel
    import wx.lib.mixins.listctrl as listmix
    from wx.lib.buttons import GenBitmapTextButton
    #from pubsub import pub

    WX_VERSION = int(wx.version()[0])
    HAS_WX = True

except ModuleNotFoundError as e:
    HAS_WX = False
    WX_VERSION = 0
    wx                  = None
    GenBitmapTextButton = None
    pub                 = None
    listmix             = None

# 
# optional import: R
# 
try:
    import rpy2.robjects
    from rpy2.robjects import (
        r,
        DataFrame,
        globalenv,
        IntVector,
        FloatVector,
        StrVector,
        packages as rpackages,
    )
    HAS_R = True
except ModuleNotFoundError as e:
    HAS_R = False
    r = None
    DataFrame   = None
    globalenv   = None
    IntVector   = None
    FloatVector = None
    StrVector   = None
    rpackages   = None

import pytransit
from pytransit.specific_tools import tnseq_tools
from pytransit.specific_tools import norm_tools
import pytransit.generic_tools.csv as csv
from pytransit.specific_tools import console_tools
from pytransit.generic_tools.lazy_dict import LazyDict
from pytransit.generic_tools.named_list import named_list
from pytransit.specific_tools.console_tools import clean_args

if HAS_WX:
    class AssumeZerosDialog(wx.Dialog):
        def __init__(self, *args, **kw):

            self.ID_HIMAR1 = wx.NewId()
            self.ID_TN5 = wx.NewId()

            wx.Dialog.__init__(self, None, title="Dialog")

            self.ID_HIMAR1 = wx.NewId()
            self.ID_TN5 = wx.NewId()

            self.SetSize((500, 300))
            self.SetTitle("Warning:  Wig Files Do Not Include Empty Sites")

            mainSizer = wx.BoxSizer(wx.VERTICAL)
            self.SetSizer(mainSizer)

            warningText = """

                One or more of your .wig files does not include any empty sites (i.e. sites with zero read-counts). The analysis methods in TRANSIT require knowing ALL possible insertion sites, even those without reads.
                    
                    Please indicate how you want to proceed:

                    As Himar1: You will need to provide the DNA sequence (.fasta format) and TRANSIT will automatically determine empty TA sites.

                    As Tn5: TRANSIT will assume all nucleotides are possible insertion sites. Those not included in the .wig file are assumed to be zero.
            """.replace("\n                ", "\n")
            warningStaticBox = wx.StaticText(
                self, wx.ID_ANY, warningText, (-1, -1), (-1, -1), wx.ALL
            )
            warningStaticBox.Wrap(480)
            mainSizer.Add(warningStaticBox, flag=wx.CENTER, border=5)

            button_sizer = wx.BoxSizer(wx.HORIZONTAL)
            himar1Button = wx.Button(self, self.ID_HIMAR1, label="Proceed as Himar1")
            tn5Button = wx.Button(self, self.ID_TN5, label="Proceed as Tn5")
            cancelButton = wx.Button(self, wx.ID_CANCEL, label="Cancel")

            button_sizer.Add(himar1Button, flag=wx.LEFT, border=5)
            button_sizer.Add(tn5Button, flag=wx.LEFT, border=5)
            button_sizer.Add(cancelButton, flag=wx.LEFT, border=5)

            mainSizer.Add(
                button_sizer, flag=wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, border=10
            )

            himar1Button.Bind(wx.EVT_BUTTON, self.OnClose)
            tn5Button.Bind(wx.EVT_BUTTON, self.OnClose)
            cancelButton.Bind(wx.EVT_BUTTON, self.OnClose)

        def OnClose(self, event):

            if self.IsModal():
                self.EndModal(event.EventObject.Id)
            else:
                self.Close()

working_directory = os.getcwd()

# 
# Results read/write
# 
if True:
    result_output_extensions = 'Common output extensions (*.tsv,*.csv,*.dat,*.txt,*.out)|*.tsv;*.csv;*.dat;*.txt;*.out;|\nAll files (*.*)|*.*'
    result_file_classes = []
    def ResultsFile(a_class):
        """
        @ResultsFile
        class File:
            @staticmethod
            def can_load(args):
                return False
        """
        if not callable(getattr(a_class, "can_load", None)):
            raise Exception(f"""Everything that usese ResultsFile should have a can_load() static method, but {a_class} does not""")
        
        result_file_classes.append(a_class)
        return a_class
    
    def file_starts_with(path, identifier):
        with open(path) as in_file:
            for line in in_file:
                if line.startswith(identifier):
                    return True
                break
        return False
    
    def read_result(path):
        result_object = None
        for FileClass in result_file_classes:
            loadable = None
            try:
                loadable = FileClass.can_load(path)
            except Exception as error:
                print(error)
            if loadable:
                result_object = FileClass(path=path)
        
        return result_object

    def write_result(*, path, file_kind, extra_info, column_names, rows, comments=[]):
        assert file_kind.isidentifier(), f"The file_kind {file_kind} must not contain whitespace or anything else that makes it an invalid var name"
        
        import datetime
        import ez_yaml
        import pytransit
        import pytransit.generic_tools.csv as csv
        from pytransit.globals import logging, gui
        from pytransit.generic_tools.misc import indent, to_pure
        ez_yaml.yaml.version = None # disable the "%YAML 1.2\n" header
        ez_yaml.yaml.width = sys.maxint if hasattr(sys, "maxint") else sys.maxsize
        extra_info = extra_info or {}
        
        original_function = ez_yaml.yaml.representer.represent_sequence
        def my_represent_sequence(*args, **kwargs):
            if len(args) >= 2:
                tag, data, *_ = args
                if isinstance(data, (tuple, list)):
                    # if all primitives
                    if not any(isinstance(each, (tuple,list,dict)) for each in data):
                        # set the format
                        kwargs["flow_style"] = True
            return original_function(*args, **kwargs)

        ez_yaml.yaml.representer.represent_sequence = my_represent_sequence
        options = dict()

        yaml_string = ""
        try:
            yaml_string = ez_yaml.to_string(extra_info, options=options)
        except Exception as error:
            try:
                yaml_string = ez_yaml.to_string(to_pure(extra_info), options=options)
            except Exception as error:
                raise Exception(f'''There was an issue with turning this value (or its contents) into a yaml string: {extra_info}''')
        
        now = str(datetime.datetime.now())
        todays_date = now[: now.rfind(".")]
        
        # 
        # write to file
        # 
        csv.write(
            path=path,
            seperator="\t",
            comment_symbol="#",
            comments=[
                file_kind, # identifier always comes first
                *comments,
                f"yaml:",
                f"    date: {todays_date}",
                f"    transit_version: {pytransit.__version__}",
                f"    app_or_command_line: {'app' if gui.is_active else 'command_line'}",
                f"    console_command: {console_tools.full_commandline_command}",
                indent(yaml_string, by="    "),
                "\t".join(column_names) # column names always last
            ],
            rows=rows,
        )

import time
class TimerAndOutputs(object):
    def __init__(self, method_name, output_paths=[], disable=False):
        self.method_name = method_name
        self.output_paths = output_paths
        self.disable = disable
    
    def __enter__(self):
        if not self.disable:
            logging.log(f"Starting {self.method_name} analysis")
        self.start_time = time.time()
        return self
    
    def __exit__(self, _, error, traceback):
        from pytransit.globals import logging, gui
        from pytransit.components import results_area
        if error is None:
            if gui.is_active:
                for each in self.output_paths:
                    if not self.disable:
                        logging.log(f"Adding File: {each}")
                    results_area.add(each)
            if not self.disable:
                logging.log(f"Finished {self.method_name} analysis in {self.duration_in_seconds:0.1f}sec")
        else:
            raise error
    
    @property
    def duration_in_seconds(self):
        return time.time() - self.start_time


def require_r_to_be_installed(required_r_packages=[]):
    if not HAS_R:
        raise Exception(f'''Error: R and rpy2 (~= 3.0) required for this operation, see: https://www.r-project.org/''')
    
    if required_r_packages:
        missing_packages = [x for x in required_r_packages if not rpackages.isinstalled(x)]
        if len(missing_packages) > 0:
            logging.error(
                "Error: Following R packages are required: %(0)s. From R console, You can install them using install.packages(c(%(0)s))"
                % ({"0": '"{0}"'.format('", "'.join(missing_packages))})
            )

def fetch_name(filepath):
    return os.path.splitext(ntpath.basename(filepath))[0]

def basename(filepath):
    return ntpath.basename(filepath)

def dirname(filepath):
    return os.path.dirname(os.path.abspath(filepath))

def validate_wig_format(wig_list):
    # Check if the .wig files include zeros or not
    status = 0
    genome = ""
    includes_zeros = tnseq_tools.check_wig_includes_zeros(wig_list)

    if sum(includes_zeros) < len(includes_zeros):
        from pytransit.globals import logging, gui
        if not gui.is_active:
            warnings.warn(
                "\nOne or more of your .wig files does not include any empty sites (i.e. sites with zero read-counts). Proceeding as if data was Tn5 (all other sites assumed to be zero)!\n"
            )
            return (2, "")

        # Else check their decision
        dlg = AssumeZerosDialog()
        result = dlg.ShowModal()
        if result == dlg.ID_HIMAR1 and gui.is_active:
            status = 1
            # Get genome
            wc = "Known Sequence Extensions (*.fna,*.fasta)|*.fna;*.fasta;|\nAll files (*.*)|*.*"
            gen_dlg = wx.FileDialog(
                gui.frame,
                message="Save file as ...",
                defaultDir=os.getcwd(),
                defaultFile="",
                wildcard=wc,
                style=wx.FD_OPEN,
            )
            if gen_dlg.ShowModal() == wx.ID_OK:
                genome = gen_dlg.GetPath()
            else:
                genome = ""

        elif result == dlg.ID_TN5:
            status = 2
            genome = ""
        else:
            status = 3
            genome = ""
    return (status, genome)

def get_pos_hash(path):
    """Returns a dictionary that maps coordinates to a list of genes that occur at that coordinate.
    
    Arguments:
        path (str): Path to annotation in .prot_table or GFF3 format.
    
    Returns:
        dict: Dictionary of position to list of genes that share that position.
    """
    filename, file_extension = os.path.splitext(path)
    if file_extension.lower() in [".gff", ".gff3"]:
        return tnseq_tools.get_pos_hash_gff(path)
    else:
        return tnseq_tools.get_pos_hash_pt(path)

def expand_var(pairs, vars, vars_by_file_list, vars_to_vals, samples, enable_logging=False):
    '''
        pairs is a list of (var,val); samples is a set; varsByFileList is a list of dictionaries mapping values to samples for each var (parallel to vars)
        recursive: keep calling till vars reduced to empty
    '''
    from pytransit.generic_tools import misc
    
    if len(vars) == 0:
        s = "%s=%s" % (pairs[0][0], pairs[0][1])
        for i in range(1, len(pairs)):
            s += " & %s=%s" % (pairs[i][0], pairs[i][1])
        s += ": %s" % len(samples)
        enable_logging and logging.log(s)
        if len(samples) == 0:
            return True
    else:
        var, vals_dict = vars[0], vars_by_file_list[0]
        inv = misc.invert_dict(vals_dict)
        any_empty = False
        for val in vars_to_vals[var]:
            subset = samples.intersection(set(inv[val]))
            res = expand_var(
                pairs + [(var, val)],
                vars[1:],
                vars_by_file_list[1:],
                vars_to_vals,
                subset,
            )
            any_empty = any_empty or res
        return any_empty

def get_validated_data(wig_list):
    """ Returns a tuple of (data, position) containing a matrix of raw read-counts
        , and list of coordinates. 

    Arguments:
        wig_list (list): List of paths to wig files.

    Returns:
        tuple: Two lists containing data and positions of the wig files given.

    :Example:

        >>> from pytransit.specific_tools import tnseq_tools
        >>> (data, position) = tnseq_tools.get_validated_data(["data/cholesterol_glycerol.transit/glycerol_rep1.wig", "data/cholesterol_glycerol.transit/glycerol_rep2.wig"])
        >>> print(data)
        array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
               [ 0.,  0.,  0., ...,  0.,  0.,  0.]])

    .. seealso:: :class:`get_file_types` :class:`combine_replicates` :class:`get_data_zero_fill` :class:`pytransit.specific_tools.norm_tools.normalize_data`
    """
    (status, genome) = validate_wig_format(wig_list)

    # Regular file with empty sites
    if status == 0:
        return tnseq_tools.CombinedWig.gather_wig_data(wig_list)
    # No empty sites, decided to proceed as Himar1
    elif status == 1:
        return tnseq_tools.get_data_w_genome(wig_list, genome)
    # No empty sites, decided to proceed as Tn5
    elif status == 2:
        return tnseq_tools.get_data_zero_fill(wig_list)
    # Didn't choose either.... what!?
    else:
        return tnseq_tools.CombinedWig.gather_wig_data([])

# input: conditions are per wig; orderingmetadata comes from tnseq_tools.CombinedWigMetadata.read_condition_data()
# output: conditionsList is selected subset of conditions (unique, in preferred order)
def filter_wigs_by_conditions(
    data,
    conditions,
    covariates=[],
    interactions=[],
    excluded_conditions=[],
    included_conditions=[],
    unknown_cond_flag="FLAG-UNMAPPED-CONDITION-IN-WIG",
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
        raise Exception(f'''Both excluded and included conditions have len > 0''')
    elif len(excluded_conditions) > 0:
        logging.log("conditions excluded: {0}".format(excluded_conditions))
        for i, c in enumerate(conditions):
            if (c != unknown_cond_flag) and (c not in excluded_conditions):
                d_filtered.append(data[i])
                cond_filtered.append(conditions[i])
                filtered_indexes.append(i)
    elif len(included_conditions) > 0:
        logging.log("conditions included: {0}".format(included_conditions))
        for i, c in enumerate(conditions):
            if (c != unknown_cond_flag) and (c in included_conditions):
                d_filtered.append(data[i])
                cond_filtered.append(conditions[i])
                filtered_indexes.append(i)
    else:
        for i, c in enumerate(conditions):
            if c != unknown_cond_flag:
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

def select_conditions(conditions, included_conditions, excluded_conditions, ordering_metadata):
    if len(included_conditions) > 0:
        conditions_list = included_conditions
    else:
        conditions_list = []
        for each_condition in ordering_metadata[
            "condition"
        ]:  # the order conds appear in metadata file, duplicated for each sample
            if each_condition not in conditions_list:
                conditions_list.append(each_condition)
    for each_condition in excluded_conditions:
        if each_condition in conditions_list:
            conditions_list.remove(each_condition)
    
    return conditions_list

def filter_wigs_by_conditions2(
    data,
    wig_fingerprints,
    condition_names,
    included_cond,
    excluded_cond,
    conditions,
    covariates=[],
    interactions=[],
):
    """
        Filters conditions that are excluded/included; also extract cond, covar, and interaction labels
        condition_names: based on original Conditions column in metadata
        conditions: user might have specified an alternative column to analyze (list of labels parallel to wigs)
    """
    (
        file_names_filtered,
        cond_names_filtered,
        d_filtered,
        cond_filtered,
        filtered_indexes,
    ) = ([], [], [], [], [])

    for i in range(len(data)):
        if (
            len(included_cond) == 0 or condition_names[i] in included_cond
        ) and condition_names[i] not in excluded_cond:
            d_filtered.append(data[i])
            file_names_filtered.append(wig_fingerprints[i])
            cond_names_filtered.append(condition_names[i])
            cond_filtered.append(conditions[i])
            filtered_indexes.append(i)

    covariates_filtered = [[c[i] for i in filtered_indexes] for c in covariates]
    interactions_filtered = [[c[i] for i in filtered_indexes] for c in interactions]

    return (
        numpy.array(d_filtered),
        numpy.array(file_names_filtered),
        numpy.array(cond_names_filtered),
        numpy.array(cond_filtered),
        numpy.array(covariates_filtered),
        numpy.array(interactions_filtered),
    )

# return a hash table of parallel lists, indexed by column header
def get_samples_metadata(metadata_path):
    data = {}
    header = None
    with open(metadata_path) as file:
        for line in file:
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

def gather_sample_data_for(conditions=None, wig_ids=None, wig_fingerprints=None, selected_samples=False):
    from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
    from pytransit.specific_tools.tnseq_tools import Wig
    
    wig_objects = gui.samples
    # default to all samples unless selected_samples is true
    if selected_samples:
        wig_objects = gui.selected_samples
    
    # filter by conditions if needed
    if conditions:
        condition_names = [ (each if isinstance(each, str) else each.name) for each in conditions ]
        wig_objects = [ each for each in wig_objects if set(each.condition_names) & set(condition_names) ]
    
    # filter by wig_ids if needed
    if wig_ids:
        wig_objects = [ each for each in wig_objects if each.id in wig_ids ]
    
    # filter by wig_fingerprints if needed
    if wig_fingerprints:
        wig_objects = [ each for each in wig_objects if each.id in wig_fingerprints ]
    
    return Wig.selected_as_gathered_data(wig_objects)

##########################################
def winsorize(counts):
    # input is insertion counts for gene: list of lists: n_replicates (rows) X n_TA sites (cols) in gene
    unique_counts = numpy.unique(numpy.concatenate(counts))
    if len(unique_counts) < 2:
        return counts
    else:
        n, n_minus_1 = unique_counts[
            heapq.nlargest(2, range(len(unique_counts)), unique_counts.take)
        ]
        result = [
            [n_minus_1 if count == n else count for count in wig] for wig in counts
        ]
        return numpy.array(result)

def stats_for_gene(site_indexes, group_wig_index_map, data, winz):
    """
        Returns a dictionary of {Group: {Mean, NzMean, NzPerc}}
        ([SiteIndex], [Condition], [WigData]) -> [{Condition: Number}]
        SiteIndex :: Number
        WigData :: [Number]
        Group :: String (combination of '<interaction>_<condition>')
    """
    nonzero = lambda xs: xs[numpy.nonzero(xs)]
    nzperc = lambda xs: numpy.count_nonzero(xs) / float(xs.size)

    means = {}
    nz_means = {}
    nz_percs = {}

    for (group, wig_indexes) in group_wig_index_map.items():
        if len(site_indexes) == 0:  # If no TA sites, write 0
            means[group] = 0
            nz_means[group] = 0
            nz_percs[group] = 0
        else:
            arr = data[wig_indexes][:, site_indexes]
            if winz:
                arr = tnseq_tools.winsorize(arr)
            means[group] = numpy.mean(arr) if len(arr) > 0 else 0
            nonzero_arr = nonzero(arr)
            nz_means[group] = numpy.mean(nonzero_arr) if len(nonzero_arr) > 0 else 0
            nz_percs[group] = nzperc(arr)

    return {"mean": means, "nz_mean": nz_means, "nz_perc": nz_percs}

def get_stats_by_rv(data, rv_site_indexes_map, genes, conditions, interactions, winz):
    """
        Returns Dictionary of Stats by condition for each Rv
        ([[Wigdata]], {Rv: SiteIndex}, [Gene], [Condition], [Interaction]) -> {Rv: {Condition: Number}}
        Wigdata :: [Number]
        SiteIndex :: Number
        Gene :: {start, end, rv, gene, strand}
        Condition :: String
        Interaction :: String
    """
    from pytransit.specific_tools.transit_tools import SEPARATOR, stats_for_gene
    
    ## Group wigfiles by (interaction, condition) pair
    ## {'<interaction>_<condition>': [Wigindexes]}
    import collections
    group_wig_index_map = collections.defaultdict(lambda: [])
    for condition_index, condition_for_wig in enumerate(conditions):
        if len(interactions) > 0:
            for interaction in interactions:
                group_name = condition_for_wig + SEPARATOR + interaction[condition_index]
                group_wig_index_map[group_name].append(condition_index)
        else:
            group_name = condition_for_wig
            group_wig_index_map[group_name].append(condition_index)
    
    stats_by_rv = {}
    for gene in genes:
        gene_name = gene["rv"]
        stats_by_rv[gene_name] = stats_for_gene(
            rv_site_indexes_map[gene_name], group_wig_index_map, data, winz,
        )

    #TODO: Any ordering to follow?
    stat_group_names = group_wig_index_map.keys()
    return stats_by_rv, stat_group_names

def calc_gene_means(combined_wig_path=None, metadata_path=None, annotation_path=None, normalization="TTR", n_terminus=0, c_terminus=0, avg_by_conditions=False, combined_wig=None):
    if combined_wig==None: 
        combined_wig = tnseq_tools.CombinedWig.load(main_path=combined_wig_path, metadata_path=metadata_path, annotation_path=annotation_path)
    
    assert combined_wig.annotation_path != None, "When computing gene means, make sure the combined_wig.annotation_path is not None"
    
    logging.log(f"Normalizing using: {normalization}")
    combined_wig = combined_wig.normalized_with(kind=normalization)
    
    genes = combined_wig.get_genes(
        n_terminus=n_terminus,
        c_terminus=c_terminus,
    )

    # calculate gene means 
    means = []
    if combined_wig.metadata_path:
        labels = combined_wig.wig_ids
    else:
        labels = combined_wig.wig_fingerprints
    
    for gene in genes:
        if gene.n>=1:
            means.append(numpy.mean(gene.reads,axis=1)) # samples are in rows; columns are TA sites in gene
        else: means.append([0]*len(labels))
    means = numpy.vstack(means)
    
    if avg_by_conditions:
        labels = combined_wig.metadata.condition_names
        condition_per_wig_index = [
            combined_wig.metadata.condition_names_for(wig_fingerprint=each_fingerprint)[0]
                for each_fingerprint in combined_wig.wig_fingerprints
        ]
        conditions_array = numpy.array(condition_per_wig_index)
        
        # make a reduced numpy array by average over replicates of each condition
        count_lists = []
        for each_condition_name in combined_wig.metadata.condition_names:
            count_lists.append(
                # pick columns corresponding to condition; avg across rows (genes)
                numpy.mean(means[:,conditions_array==each_condition_name],axis=1) 
            )
        means = numpy.array(count_lists).transpose()
    
    return means, genes, labels
