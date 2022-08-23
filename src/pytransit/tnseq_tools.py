import sys
import os
import math
import warnings
from functools import total_ordering
from collections import namedtuple
from os.path import isabs, isfile, isdir, join, dirname, basename, exists, splitext, relpath

import numpy
import scipy.stats
import ez_yaml
from super_map import LazyDict, stringify
import pytransit.basics.csv as csv
from pytransit.basics.named_list import named_list
from pytransit.basics.misc import line_count_of, flatten_once, no_duplicates, indent

try:
    from pytransit import norm_tools

    noNorm = False
except ImportError:
    noNorm = True
    warnings.warn(
        "Problem importing the norm_tools.py module. Read-counts will not be normalized. Some functions may not work."
    )


def rv_siteindexes_map(genes, TASiteindexMap, n_terminus=0.0, c_terminus=0.0):
    """
    ([Gene], {TAsite: Siteindex}) -> {Rv: Siteindex}
    """
    RvSiteindexesMap = {}
    for g, gene in enumerate(genes):
        siteindexes = []
        start = gene["start"] if gene["strand"] == "+" else gene["start"] + 3
        end = gene["end"] - 3 if gene["strand"] == "+" else gene["end"]
        for i in range(start, end + 1):
            co = i
            if (co - start) / float(end - start) < (n_terminus / 100.0):
                continue
            if (co - start) / float(end - start) > ((100 - c_terminus) / 100.0):
                continue
            if co in TASiteindexMap:
                siteindexes.append(TASiteindexMap[co])
        RvSiteindexesMap[gene["rv"]] = siteindexes
    return RvSiteindexesMap


# format:
#   header lines (prefixed by '#'), followed by lines with counts
#   counts lines contain the following columns: TA coord, counts, other info like gene/annotation
#   for each column of counts, there must be a header line prefixed by "#File: " and then an id or filename

class Wig:
    def __init__(self, path, rows=None, extra_data=None):
        self.path            = path
        self.rows            = rows or []
        self.comments        = []
        
        from random import random
        self.extra_data = LazyDict(extra_data)
        if self.extra_data.get("condition", None) is None:
            self.extra_data["condition"] = basename(path)
        if self.extra_data.get("id", None) is None:
            self.extra_data["id"] = f"{self.extra_data.get('condition', None)}{random()}".replace(".", "")
        
    def load(self):
        pass # TODO
    
    def __repr__(self):
        return f"""Wig(
            path={self.path},
            rows_shape=({len(self.rows)}, {len(self.rows[0])}),
            extra_data={indent(self.extra_data, by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")

class Condition:
    def __init__(self, name, is_disabled=False, extra_data=None):
        self.name = name
        self.is_disabled = is_disabled
        self.extra_data = LazyDict(extra_data or {})
    
    def __repr__(self):
        return f"""Condition(
            name={self.name},
            is_disabled={self.is_disabled},
            extra_data={indent(self.extra_data, by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")
    
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, other):
        return self.name == other.name

class CombinedWigMetadata:
    _CacheClass = named_list([ 'conditions_by_file','covariates_by_file_list','interactions_by_file_list','ordering_metadata',])
    
    # can contain more data
    def __init__(self, path):
        self.path = path
        self.headers = []
        self.rows = []
        self.comments, self.headers, self.rows = csv.read(self.path, seperator="\t", first_row_is_headers=True)
        self.conditions = no_duplicates(
            Condition(
                name=each_row["Condition"],
            )
                for each_row in self.rows
        )
        self._cache_for_read = None
    
    def condition_for(self, wig_path=None, id=None):
        if wig_path:
            for each_row in self.rows:
                f = each_row["Filename"]
                if each_row["Filename"] == wig_path:
                    return each_row["Condition"]
        if id:
            for each_row in self.rows:
                if each_row["Id"] == id:
                    return each_row["Condition"]
    
    def id_for(self, wig_path=None):
        for each_row in self.rows:
            if each_row["Filename"] == wig_path:
                return each_row["Id"]
    
    def wigs_for(self, condition):
        return [
            each_row["Filename"]
                for each_row in self.rows
                    if each_row["Condition"] == condition
        ]
    
    def __repr__(self):
        return f"""CWigMetadata(
            path={self.path},
            rows_shape=({len(self.rows)}, {len(self.headers)}),
            conditions={indent(self.conditions, by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")
    
    # 
    # an "always filled" cache
    # 
    @property
    def cache_for_read(self):
        if not self._cache_for_read:
            self._cache_for_read = self.read_condition_data(self.path)
        return self._cache_for_read
    @cache_for_read.setter
    def cache_for_read(self, value):
        self._cache_for_read = value
    
    @classmethod
    def read_condition_data(self, path, covars_to_read=[], interactions_to_read=[], condition_name="Condition"):
        """
        Filename -> ConditionMap
        ConditionMap :: {Filename: Condition}, [{Filename: Covar}], [{Filename: Interaction}]
        Condition :: String
        Covar :: String
        Interaction :: String
        """
        _cache_for_read = getattr(self, "_cache_for_read", None)
        if _cache_for_read:
            return _cache_for_read
        
        conditions_by_file        = {}
        covariates_by_file_list   = [{} for i in range(len(covars_to_read))]
        interactions_by_file_list = [{} for i in range(len(interactions_to_read))]
        headers_to_read           = [condition_name.lower(), "filename"]
        ordering_metadata         = {"condition": [], "interaction": []}
        with open(path) as mfile:
            lines = mfile.readlines()
            head_indexes = [
                i
                for h in headers_to_read
                for i, c in enumerate(lines[0].split())
                if c.lower() == h.lower()
            ]
            covar_indexes = [
                i
                for h in covars_to_read
                for i, c in enumerate(lines[0].split())
                if c.lower() == h.lower()
            ]
            interaction_indexes = [
                i
                for h in interactions_to_read
                for i, c in enumerate(lines[0].split())
                if c.lower() == h.lower()
            ]

            for line in lines[1:]:
                if line[0] == "#":
                    continue
                vals = line.split()
                [condition, wfile] = vals[head_indexes[0]], vals[head_indexes[1]]
                conditions_by_file[wfile] = condition
                ordering_metadata["condition"].append(condition)
                for i, c in enumerate(covars_to_read):
                    covariates_by_file_list[i][wfile] = vals[covar_indexes[i]]
                for i, c in enumerate(interactions_to_read):
                    interactions_by_file_list[i][wfile] = vals[interaction_indexes[i]]

                    # This makes sense only if there is only 1 interaction variable
                    # For multiple interaction vars, may have to rethink ordering.
                    ordering_metadata["interaction"].append(vals[interaction_indexes[i]])

        return self._CacheClass([
            conditions_by_file,
            covariates_by_file_list,
            interactions_by_file_list,
            ordering_metadata,
        ])
    
    # 
    # conditions_by_file
    # 
    @property
    def conditions_by_file(self): return self.cache_for_read.conditions_by_file
    @conditions_by_file.setter
    def conditions_by_file(self, value):
        self.cache_for_read.conditions_by_file = value
    
    # 
    # covariates_by_file_list
    # 
    @property
    def covariates_by_file_list(self): return self.cache_for_read.covariates_by_file_list
    @covariates_by_file_list.setter
    def covariates_by_file_list(self, value):
        self.cache_for_read.covariates_by_file_list = value
    
    # 
    # interactions_by_file_list
    # 
    @property
    def interactions_by_file_list(self): return self.cache_for_read.interactions_by_file_list
    @interactions_by_file_list.setter
    def interactions_by_file_list(self, value):
        self.cache_for_read.interactions_by_file_list = value
    
    # 
    # ordering_metadata
    # 
    @property
    def ordering_metadata(self): return self.cache_for_read.ordering_metadata
    @ordering_metadata.setter
    def ordering_metadata(self, value):
        self.cache_for_read.ordering_metadata = value

class CombinedWigData(named_list(['sites','counts_by_wig','files',])):
    @classmethod
    def load(cls, file_path):
        """
            Read the combined wig-file generated by Transit
            :: Filename -> Tuple([Site], [WigData], [Filename])
            Site :: Integer
            WigData :: [Number]
            Filename :: String
        """
        import pytransit.transit_tools as transit_tools
        
        sites, counts_by_wig, files, extra_data = [], [], [], {}
        contained_yaml_data = False
        number_of_lines = line_count_of(file_path)
        with open(file_path) as f:
            yaml_mode_is_on = False
            yaml_string = "extra_data:\n"
            lines = f.readlines()
            comment_lines = []
            # 
            # handle header/comments
            # 
            for line in lines:
                if line.startswith("#"):
                    if line.startswith("#yaml:"):
                        yaml_mode_is_on = True
                        contained_yaml_data = True
                        continue
                    if yaml_mode_is_on and line.startswith("# "):
                        yaml_string += f"\n{line[-1:]}"
                        continue
                    else:
                        yaml_mode_is_on = False
                        # add to the extra_data dict when its done
                        if len(yaml_string) > 0:
                            an_object = ez_yaml.to_object(string=yaml_string)
                            extra_data.update(an_object["extra_data"] or {})
                            files += extra_data.get('files',[])
                    # 
                    # handle older file method
                    # 
                    if line.startswith("#File: "):
                        files.append(line.rstrip()[7:])  # allows for spaces in filenames
                        continue
            
            # 
            # handle body
            # 
            counts_by_wig = [ [] for _ in files ]
            for index, line in enumerate(lines):
                if index % 150 == 0: # 150 is arbitrary, bigger = slower visual update but faster read
                    percent_done = round(index/number_of_lines * 100, ndigits=2)
                    transit_tools.log(f"\rreading lines: {percent_done}%          ", end="")
                
                # lines to skip
                if line.startswith("#") or len(line) == 0:
                    continue
                
                #
                # actual parsing
                #
                cols = line.split("\t")[0 : 1+len(files)]
                cols = cols[: 1+len(files)]  # additional columns at end could contain gene info
                # Read in position as int, and readcounts as float
                cols = [
                    int(each_t_iv) if index == 0 else float(each_t_iv)
                        for index, each_t_iv in enumerate(cols)
                ]
                position, wig_counts = cols[0], cols[1:]
                sites.append(position)
                for index, count in enumerate(wig_counts):
                    counts_by_wig[index].append(count)
            transit_tools.log(f"\rreading lines: 100%          ")
        
        return CombinedWigData((numpy.array(sites), numpy.array(counts_by_wig), files))

class CombinedWig:
    def __init__(self, *, main_path, metadata_path, comments=None, extra_data=None):
        self.main_path     = main_path
        self.metadata_path = metadata_path
        self.metadata      = CombinedWigMetadata(self.metadata_path)
        self.data          = CombinedWigData.load(self.main_path) # for backwards compatibility (otherwise just used self.rows and helper methods)
        self.rows          = []
        self.comments      = comments or []
        self.extra_data    = LazyDict(extra_data or {})
        self.samples       = []
        self._load_main_path()
    
    def __repr__(self):
        return f"""CombinedWig(
            path={self.main_path},
            rows_shape=({len(self.rows)}, {len(self.rows[0])}),
            extra_data={indent(self.extra_data, by="            ", ignore_first=True)},
            samples={   indent(self.samples   , by="            ", ignore_first=True)},
            metadata={  indent(self.metadata  , by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")
    
    # 
    # sites
    # 
    @property
    def sites(self):
        return [ each.position for each in self.rows ]
    
    # 
    # files
    # 
    @property
    def files(self): return self.extra_data["files"]
    @files.setter
    def files(self, value):
        self.extra_data["files"] = value
    
    # 
    # conditions
    # 
    @property
    def conditions(self):
        return self.metadata.conditions
    
    # 
    # read_counts_array
    # 
    @property
    def read_counts_array(self):
        return [
            each_row[1:len(self.files)] 
                for each_row in self.rows 
        ]
    
    # 
    # read_counts_wig
    # 
    @property
    def read_counts_by_wig(self):
        counts_for_wig = { each_path: [] for each_path in self.files }
        for each_row in self.rows:
            for each_wig_path in self.files:
                counts_for_wig[each_wig_path].append(
                    each_row[each_wig_path]
                )
        return counts_for_wig
    
    def _load_main_path(self):
        comments, headers, rows = csv.read(self.main_path, seperator="\t", first_row_is_headers=False)
        comment_string = "\n".join(comments)
        
        sites, counts_by_wig, files, extra_data = [], [], [], {}
        yaml_switch_is_on = False
        yaml_string = "extra_data:\n"
        for line in comments:
            # 
            # handle yaml
            # 
            if line.startswith("yaml:"):
                yaml_switch_is_on = True
                continue
            if yaml_switch_is_on and line.startswith(" "):
                yaml_string += f"\n"+line
                continue
            else:
                yaml_switch_is_on = False
                
            # 
            # handle older file method
            # 
            if line.startswith("File: "):
                files.append(line.rstrip()[6:])  # allows for spaces in filenames
                continue
            
            # 
            # save other comments
            # 
            self.comments.append(line)
        
        # 
        # process yaml
        # 
        if len(yaml_string) > len("extra_data:\n"):
            extra = ez_yaml.to_object(string=yaml_string)["extra_data"]
            self.extra_data.update(
                extra or {}
            )
            files += self.extra_data.get('files',[])
            no_duplicates = []
            # remove any duplicate entries while preserving order
            for each_file in files:
                if each_file not in no_duplicates:
                    no_duplicates.append(each_file)
            files = no_duplicates
        
        # 
        # process data
        # 
        extra_data["files"] = files
        self.extra_data.update(extra_data)
        CWigRow = named_list([ "position", *files ])
        for row in rows:
            #
            # extract
            #
            position   = row[0]
            read_counts = row[ 1: 1+len(files) ]
            other      = row[ 1+len(files) :  ] # often empty
            
            # force types
            position   = int(position)
            read_counts = [ float(each) for each in read_counts ]
            
            # save
            self.rows.append(CWigRow([position]+read_counts+other))
            
            sites.append(position)
            for index, count in enumerate(read_counts):
                if len(counts_by_wig) < index+1:
                    counts_by_wig.append([])
                counts_by_wig[index].append(count)
        
        read_counts_by_wig = self.read_counts_by_wig
        for each_path in self.files:
            self.samples.append(
                Wig(
                    path=each_path,
                    rows=list(zip(self.sites, read_counts_by_wig[each_path])),
                    extra_data=LazyDict(
                        is_part_of_cwig=True,
                    ),
                )
            )
        
        return self
    
    @classmethod
    def gather_wig_data(cls, list_of_paths):
        """ Returns a tuple of (data, position) containing a matrix of raw read-counts
            , and list of coordinates.

        Arguments:
            wig_list (list): List of paths to wig files.

        Returns:
            tuple: Two lists containing data and positions of the wig files given.

        :Example:

            >>> from pytransit.tnseq_tools import CombinedWig
            >>> (data, position) = CombinedWig.gather_wig_data(["data/glycerol_H37Rv_rep1.wig", "data/glycerol_H37Rv_rep2.wig"])
            >>> print(data)
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                [ 0.,  0.,  0., ...,  0.,  0.,  0.]])

        .. seealso:: :class:`get_file_types` :class:`combine_replicates` :class:`get_data_zero_fill` :class:`pytransit.norm_tools.normalize_data`
        """
        # If empty just quickly return empty lists
        if not list_of_paths:
            return (numpy.zeros((1, 0)), numpy.zeros(0))

        # Check size of all wig file matches
        size_list = []
        for path in list_of_paths:
            line_count = 0
            for line in open(path):
                if line[0] not in "0123456789":
                    continue
                line_count += 1
            size_list.append(line_count)

        # If it doesn't match, report an error and quit
        if sum(size_list) != (line_count * len(size_list)):
            raise Exception(f'''
            
                Error:
                    Not all wig files have the same number of sites.
                    Make sure all .wig files come from the same strain.
                
                counts: {list(zip(list_of_paths, size_list))}
            ''')

        data_per_path = numpy.zeros((len(list_of_paths), line_count))
        position = numpy.zeros(line_count, dtype=int)
        for path_index, path in enumerate(list_of_paths):
            line_index = 0
            prev_pos   = 0
            for line in open(path):
                if line[0] not in "0123456789": # not sure why this is here -- Jeff
                    continue
                
                tmp      = line.split()
                pos      = int(tmp[0])
                rd       = float(tmp[1])
                prev_pos = pos

                try:
                    data_per_path[path_index, line_index] = rd
                except Exception as error:
                    raise Exception(f'''
                        
                        Make sure that all wig files have the same number of TA sites (i.e. same strain)
                        
                        Original Error:\n{error}
                        
                    ''')
                position[line_index] = pos
                line_index += 1
        
        return data_per_path, position
    
# backwards compatibility
read_samples_metadata = CombinedWigMetadata.read_condition_data
read_combined_wig = CombinedWigData.load

def read_genes(fname, descriptions=False):
    """
      (Filename, Options) -> [Gene]
      Gene :: {start, end, rv, gene, strand}
    """
    genes = []
    for line in open(fname):
        w = line.rstrip().split("\t")
        data = {
            "start": int(w[1]),
            "end": int(w[2]),
            "rv": w[8],
            "gene": w[7],
            "strand": w[3],
        }
        if descriptions == True:
            data.append(w[0])
        genes.append(data)
    return genes


@total_ordering
class Gene:
    """Class defining a gene with useful attributes for TnSeq analysis.

    This class helps define a "gene" with attributes that facilitate TnSeq
    analysis. Here "gene" can be defined to be any genomic region. The Genes
    class (with an s) can be used to define list of Gene objects with more
    useful operations on the "genome" level.

    Attributes:
        orf: A string defining the ID of the gene.
        name: A string with the human readable name of the gene.
        desc: A string with the description of the gene.
        reads: List of lists of read-counts in possible site replicate dataset.
        position: List of coordinates of the possible sites.
        start: An integer defining the start coordinate for the gene.
        end: An integer defining the end coordinate for the gene.
        strand: A string defining the strand of the gene.


    :Example:

        >>> import pytransit.tnseq_tools as tnseq_tools
        >>> G = tnseq_tools.Gene("Rv0001", "dnaA", "DNA Replication A", [[0,0,0,0,1,3,0,1]],  [1,21,32,37,45,58,66,130], strand="+" )
        >>> print(G)
        Rv0001  (dnaA)  k=3 n=8 r=4 theta=0.37500
        >>> print(G.phi())
        0.625
        >>> print(G.tosses)
        array([ 0.,  0.,  0.,  0.,  1.,  1.,  0.,  1.])

        .. seealso:: :class:`Genes`
        """

    def __init__(self, orf, name, desc, reads, position, start=0, end=0, strand=""):
        """Initializes the Gene object.

        Arguments:
            orf (str): A string defining the ID of the gene.
            name (str): A string with the human readable name of the gene.
            desc (str): A string with the description of the gene.
            reads (list): List of lists of read-counts in possible site replicate dataset.
            position (list): List of coordinates of the possible sites.
            start (int): An integer defining the start coordinate for the gene.
            end (int): An integer defining the end coordinate for the gene.
            strand (str): A string defining the strand of the gene.

        Returns:
            Gene: Object of the Gene class with the defined attributes.

        """

        self.orf = orf
        self.name = name
        self.desc = desc
        self.start = start
        self.end = end
        self.strand = strand
        self.reads = numpy.array(reads)
        self.position = numpy.array(position, dtype=int)
        self.tosses = tossify(self.reads)
        try:
            self.runs = runs(self.tosses)
        except Exception as e:
            print(orf, name, self.tosses)
            raise e

        self.k = int(numpy.sum(self.tosses))
        self.n = len(self.tosses)
        try:
            self.r = numpy.max(self.runs)
        except Exception as e:
            print(orf, name, self.tosses)
            print(self.runs)
            raise e

        self.s = self.get_gap_span()
        self.t = self.get_gene_span()

    #

    def __getitem__(self, i):
        """Return read-counts at position i.

        Arguments:
            i (int): integer of the index of the desired site.

        Returns:
            list: Reads at position i.
        """
        return self.reads[:, i]

    #

    def __str__(self):
        """Return a string representation of the object.

        Returns:
            str: Human readable string with some of the attributes.
        """
        return "%s\t(%s)\tk=%d\tn=%d\tr=%d\ttheta=%1.5f" % (
            self.orf,
            self.name,
            self.k,
            self.n,
            self.r,
            self.theta(),
        )

    #

    def __eq__(self, other):
        """Compares against other gene object.

        Returns:
            bool: True if the gene objects have same orf id.
        """
        return self.orf == other.orf

    #

    def __lt__(self, other):
        """Compares against other gene object.

        Returns:
            bool: True if the gene object id is less than the other.
        """
        return self.orf < other.orf

    #

    def get_gap_span(self):
        """Returns the span of the maxrun of the gene (i.e. number of nucleotides).

        Returns:
            int: Number of nucleotides spanned by the max run.
        """
        if len(self.position) > 0:
            if self.r == 0:
                return 0
            index = runindex(self.runs)
            # maxii = numpy.argmax(self.runs)
            maxii = numpy.argwhere(self.runs == numpy.max(self.runs)).flatten()[-1]
            runstart = index[maxii]
            runend = runstart + max(self.runs) - 1
            return self.position[runend] - self.position[runstart] + 2
        else:
            return 0

    #

    def get_gene_span(self):
        """Returns the number of nucleotides spanned by the gene.

        Returns:
            int: Number of nucleotides spanned by the gene's sites.
        """
        if len(self.position) > 0:
            return self.position[-1] - self.position[0] + 2
        return 0

    #

    def theta(self):
        """Return the insertion density ("theta") for the gene.

        Returns:
            float: Density of the gene (i.e. k/n )
        """
        if self.n:
            return float(self.k) / self.n
        else:
            return 0.0

    #

    def phi(self):
        """Return the non-insertion density ("phi") for the gene.

        Returns:
            float: Non-insertion density  (i.e. 1 - theta)
        """
        return 1.0 - self.theta()

    #

    def total_reads(self):
        """Return the total reads for the gene.

        Returns:
            float: Total sum of read-counts.
        """
        return numpy.sum(self.reads, 1)


#


class Genes:
    """Class defining a list of Gene objects with useful attributes for TnSeq
    analysis.

    This class helps define a list of Gene objects with attributes that
    facilitate TnSeq analysis. Includes methods that calculate useful statistics
    and even rudamentary analysis of essentiality.

    Attributes:
        wig_list: List of paths to datasets in .wig format.
        protTable: String with path to annotation in .prot_table format.
        norm: String with the normalization used/
        reps: String with information on how replicates were handled.
        minread: Integer with the minimum magnitude of read-count considered.
        ignore_codon: Boolean defining whether to ignore the start/stop codon.
        n_terminus: Float number of the fraction of the N-terminus to ignore.
        c_terminus: Float number of the fraction of the C-terminus to ignore.
        include_nc: Boolean determining whether to include non-coding areas.
        orf2index: Dictionary of orf id to index in the genes list.
        genes: List of the Gene objects.


    :Example:

        >>> import pytransit.tnseq_tools as tnseq_tools
        >>> G = tnseq_tools.Genes(["transit/data/glycerol_H37Rv_rep1.wig", "transit/data/glycerol_H37Rv_rep2.wig"], "transit/genomes/H37Rv.prot_table", norm="TTR")
        >>> print(G)
        Genes Object (N=3990)
        >>> print(G.global_theta())
        0.40853707222816626
        >>> print(G["Rv0001"]   # Lookup like dictionary)
        Rv0001  (dnaA)  k=0 n=31    r=31    theta=0.00000
        >>> print(G[2]          # Lookup like list)
        Rv0003  (recF)  k=5 n=35    r=14    theta=0.14286
        >>> print(G[2].reads)
        [[  62.            0.            0.            0.            0.            0.
         0.            0.            0.            0.            0.            0.
         0.            0.           63.            0.            0.           13.
        46.            0.            1.            0.            0.            0.
         0.            0.            0.            0.            0.            0.
         0.            0.            0.            0.            0.        ]
         [   3.14314432   67.26328843    0.            0.            0.            0.
         0.            0.            0.           35.20321637    0.            0.
         0.            0.           30.80281433    0.          101.20924707
         0.           23.25926796    0.           16.97297932    8.17217523
         0.            0.            2.51451546    3.77177318    0.62862886
         0.            0.           69.14917502    0.            0.            0.
         0.            0.        ]]


        .. seealso:: :class:`Gene`
        """

    #

    def __getitem__(self, i):
        """Defines __getitem__ method so that it works as dictionary and list.

        Arguments:
            i (int): Integer or string defining index or orf ID desired.
        Returns:
            Gene: A gene with the index or ID equal to i.
        """
        if isinstance(i, int):
            return self.genes[i]

        if isinstance(i, str):
            return self.genes[self.orf2index[i]]

    #

    def __contains__(self, item):
        """Defines __contains__ to check if gene exists in the list.

        Arguments:
            item (str): String with the id of the gene.

        Returns:
            bool: Boolean with True if item is in the list.
        """
        return item in self.orf2index

    #

    def __len__(self):
        """Defines __len__ returning number of genes.

        Returns:
            int: Number of genes in the list.
        """
        return len(self.genes)

    #

    def __str__(self):
        """Defines __str__ to print(a generic str with the size of the list.)

        Returns:
            str: Human readable string with number of genes in object.
        """
        return "Genes Object (N=%d)" % len(self.genes)
    
    def __repr__(self):
        return f"""{LazyDict(
            wig_list=self.wig_list,
            annotation=self.annotation,
            norm=self.norm,
            reps=self.reps,
            minread=self.minread,
            ignore_codon=self.ignore_codon,
            n_terminus=self.n_terminus,
            c_terminus=self.c_terminus,
            include_nc=self.include_nc,
        )}"""

    #

    def __init__(
        self,
        wig_list,
        annotation,
        norm="nonorm",
        reps="All",
        minread=1,
        ignore_codon=True,
        n_terminus=0.0,
        c_terminus=0.0,
        include_nc=False,
        data=[],
        position=[],
        genome="",
        transposon="himar1",
    ):
        """Initializes the gene list based on the list of wig files and a prot_table.

        This class helps define a list of Gene objects with attributes that
        facilitate TnSeq analysis. Includes methods that calculate useful statistics
        and even rudamentary analysis of essentiality.

        Arguments:
            wig_list (list): List of paths to datasets in .wig format.
            protTable (str): String with path to annotation in .prot_table format.
            norm (str): String with the normalization used/
            reps (str): String with information on how replicates were handled.
            minread (int): Integer with the minimum magnitude of read-count considered.
            ignore_codon (bool): Boolean defining whether to ignore the start/stop codon.
            n_terminus (float): Float number of the fraction of the N-terminus to ignore.
            c_terminus (float): Float number of the fraction of the C-terminus to ignore.
            include_nc (bool): Boolean determining whether to include non-coding areas.
            data (list): List of data. Used to define the object without files.
            position (list): List of position of sites. Used to define the object without files.


        """
        self.wig_list = wig_list
        self.annotation = annotation
        self.norm = norm
        self.reps = reps
        self.minread = minread
        self.ignore_codon = ignore_codon
        self.n_terminus = n_terminus
        self.c_terminus = c_terminus
        self.include_nc = include_nc

        isProt = True
        filename, file_extension = os.path.splitext(self.annotation)
        if file_extension.lower() in [".gff", ".gff3"]:
            isProt = False

        self.orf2index = {}
        self.genes = []

        orf2info = get_gene_info(self.annotation)
        if not numpy.any(data):
            if transposon.lower() == "himar1" and not genome:
                (data, position) = CombinedWig.gather_wig_data(self.wig_list)
            elif genome:
                (data, position) = get_data_w_genome(self.wig_list, genome)
            else:
                (data, position) = get_data_zero_fill(self.wig_list)

        ii_min = data < self.minread
        data[ii_min] = 0

        hash = get_pos_hash(self.annotation)

        if not noNorm:
            (data, factors) = norm_tools.normalize_data(
                data, norm, self.wig_list, self.annotation
            )
        else:
            factors = []

        if reps.lower() != "all":
            data = numpy.array([combine_replicates(data, method=reps)])

        K, N = data.shape

        self.data = data
        orf2posindex = {}
        visited_list = []
        for i in range(N):
            genes_with_coord = hash.get(position[i], [])
            for gene in genes_with_coord:
                if gene not in orf2posindex:
                    visited_list.append(gene)
                if gene not in orf2posindex:
                    orf2posindex[gene] = []

                name, desc, start, end, strand = orf2info.get(gene, ["", "", 0, 0, "+"])

                if strand == "+":
                    if self.ignore_codon and position[i] > end - 3:
                        continue
                else:
                    if self.ignore_codon and position[i] < start + 3:
                        continue

                if (position[i] - start) / float(end - start) < (self.n_terminus / 100.0):
                    continue

                if (position[i] - start) / float(end - start) > (
                    (100 - self.c_terminus) / 100.0
                ):
                    continue

                orf2posindex[gene].append(i)

        count = 0
        for line in open(self.annotation):
            if line.startswith("#"):
                continue
            tmp = line.split("\t")

            if isProt:
                gene = tmp[8].strip()
                name, desc, start, end, strand = orf2info.get(gene, ["", "", 0, 0, "+"])
            else:
                features = dict(
                    [
                        tuple(f.split("=", 1))
                        for f in filter(lambda x: "=" in x, tmp[8].split(";"))
                    ]
                )
                gene = features["ID"]
                name, desc, start, end, strand = orf2info.get(gene, ["", "", 0, 0, "+"])
            posindex = orf2posindex.get(gene, [])
            if posindex:
                pos_start = orf2posindex[gene][0]
                pos_end = orf2posindex[gene][-1]
                self.genes.append(
                    Gene(
                        gene,
                        name,
                        desc,
                        data[:, pos_start : pos_end + 1],
                        position[pos_start : pos_end + 1],
                        start,
                        end,
                        strand,
                    )
                )
            else:
                self.genes.append(
                    Gene(
                        gene,
                        name,
                        desc,
                        numpy.array([[]]),
                        numpy.array([]),
                        start,
                        end,
                        strand,
                    )
                )
            self.orf2index[gene] = count
            count += 1

    #

    def local_insertions(self):
        """Returns numpy array with the number of insertions, 'k', for each gene.

        Returns:
            narray: Numpy array with the number of insertions for all genes.
        """
        G = len(self.genes)
        K = numpy.zeros(G)
        for i in range(G):
            K[i] = self.genes[i].k
        return K

    #

    def local_sites(self):
        """Returns numpy array with total number of TA sites, 'n', for each gene.

        Returns:
            narray: Numpy array with the number of sites for all genes.
        """
        G = len(self.genes)
        N = numpy.zeros(G)
        for i in range(G):
            N[i] = self.genes[i].n
        return N

    #

    def local_runs(self):
        """Returns numpy array with maximum run of non-insertions, 'r', for each gene.

        Returns:
            narray: Numpy array with the max run of non-insertions for all genes.
        """
        G = len(self.genes)
        R = numpy.zeros(G)
        for i in range(G):
            R[i] = self.genes[i].r
        return R

    #

    def local_gap_span(self):
        """Returns numpy array with the span of nucleotides of the largest gap,
        's', for each gene.

        Returns:
            narray: Numpy array with the span of gap for all genes.
        """
        G = len(self.genes)
        S = numpy.zeros(G)
        for i in range(G):
            S[i] = self.genes[i].s
        return S

    #

    def local_gene_span(self):
        """Returns numpy array with the span of nucleotides of the gene,
        't', for each gene.

        Returns:
            narray: Numpy array with the span of gene for all genes.
        """
        G = len(self.genes)
        T = numpy.zeros(G)
        for i in range(G):
            T[i] = self.genes[i].t
        return T

    #

    def local_reads(self):
        """Returns numpy array of lists containing the read counts for each gene.

        Returns:
            narray: Numpy array with the list of reads for all genes.
        """
        all_reads = []
        G = len(self.genes)
        for i in range(G):
            all_reads.extend(self.genes[i].reads)
        return numpy.array(all_reads)

    #

    def local_thetas(self):
        """Returns numpy array of insertion frequencies, 'theta', for each gene.

        Returns:
            narray: Numpy array with the density for all genes.
        """
        G = len(self.genes)
        theta = numpy.zeros(G)
        for i in range(G):
            theta[i] = self.genes[i].theta()
        return theta

    #

    def local_phis(self):
        """Returns numpy array of non-insertion frequency, 'phi', for each gene.

        Returns:
            narray: Numpy array with the complement of density for all genes.
        """
        return 1.0 - self.theta()

    #

    def global_insertion(self):
        """Returns total number of insertions, i.e. sum of 'k' over all genes.

        Returns:
            float: Total sum of reads across all genes.
        """
        G = len(self.genes)
        total = 0
        for i in range(G):
            total += self.genes[i].k
        return total

    #

    def global_sites(self):
        """Returns total number of sites, i.e. sum of 'n' over all genes.

        Returns:
            int: Total number of sites across all genes.
        """
        G = len(self.genes)
        total = 0
        for i in range(G):
            total += self.genes[i].n
        return total

    #

    def global_run(self):
        """Returns the run assuming all genes were concatenated together.

        Returns:
            int: Max run across all genes.
        """
        return maxrun(self.tosses())

    #

    def global_reads(self):
        """Returns the reads among the library.

        Returns:
            list: List of all the data.
        """
        return self.data

    #

    def global_theta(self):
        """Returns global insertion frequency, of the library.

        Returns:
            float: Total sites with insertions divided by total sites.
        """
        return float(self.global_insertion()) / self.global_sites()

    #

    def global_phi(self):
        """Returns global non-insertion frequency, of the library.

        Returns:
            float: Complement of global theta i.e. 1.0-theta
        """
        return 1.0 - self.global_theta()

    #

    def total_reads(self):
        """Returns total reads among the library.

        Returns:
            float: Total sum of read-counts accross all genes.
        """
        reads_total = 0
        for g in self.genes:
            reads_total += g.total_reads()
        return reads_total

    #

    def tosses(self):
        """Returns list of bernoulli trials, 'tosses', representing insertions in the gene.

        Returns:
            list: Sites represented as bernoulli trials with insertions as true.
        """
        all_tosses = []
        for g in self.genes:
            all_tosses.extend(g.tosses)
        return all_tosses


#


def tossify(data):
    """Reduces the data into Bernoulli trials (or 'tosses') based on whether counts were observed or not.

    Arguments:
        data (list): List of numeric data.

    Returns:
        list: Data represented as bernoulli trials with >0 as true.
    """
    K, N = data.shape
    reduced = numpy.sum(data, 0)
    return numpy.zeros(N) + (numpy.sum(data, 0) > 0)


#


def runs(data):
    """Return list of all the runs of consecutive non-insertions.

    Arguments:
        data (list): List of numeric data.

    Returns:
        list: List of the length of the runs of non-insertions. Non-zero sites are treated as runs of zero.
    """
    runs = []
    current_r = 0
    for read in data:
        if read > 0:  # If ending a run of zeros
            if current_r > 0:  # If we were in a run, add to list
                runs.append(current_r)
            current_r = 0
            runs.append(current_r)
        else:
            current_r += 1
    # If we ended in a run, add it
    if current_r > 0:
        runs.append(current_r)

    if not runs:
        return [0]
    return runs


#


def runindex(runs):
    """Returns a list of the indexes of the start of the runs; complements runs().

    Arguments:
        runs (list): List of numeric data.

    Returns:
        list: List of the index of the runs of non-insertions. Non-zero sites are treated as runs of zero.
    """
    index = 0
    index_list = []
    runindex = 0
    for r in runs:
        for i in range(r):
            if i == 0:
                runindex = index
            index += 1
        if r == 0:
            runindex = index
            index += 1
        index_list.append(runindex)
    return index_list


#


def get_file_types(wig_list):
    """Returns the transposon type (himar1/tn5) of the list of wig files.

    Arguments:
        wig_list (list): List of paths to wig files.

    Returns:
        list: List of transposon type ("himar1" or "tn5").
    """
    if not wig_list:
        return []

    types = ["tn5" for i in range(len(wig_list))]
    for i, wig_filename in enumerate(wig_list):
        with open(wig_filename) as wig_file:
            prev_pos = 0
            for line in wig_file:
                if line[0] not in "0123456789":
                    continue
                tmp = line.split()
                pos = int(tmp[0])
                rd = float(tmp[1])
                if pos != prev_pos + 1:
                    types[i] = "himar1"
                    break
                prev_pos = pos
    return types


def check_wig_includes_zeros(wig_list):
    """Returns boolean list showing whether the given files include empty sites
    (zero) or not.

    Arguments:
        wig_list (list): List of paths to wig files.

    Returns:
        list: List of boolean values.
    """
    if not wig_list:
        return []
    includes = [False for i in range(len(wig_list))]
    for i, wig_filename in enumerate(wig_list):
        with open(wig_filename) as wig_file:
            for line in wig_file:
                if line[0] not in "0123456789":
                    continue
                tmp = line.split()
                pos = int(tmp[0])
                rd = float(tmp[1])
                if rd == 0:
                    includes[i] = True
                    break
    return includes


#


def get_unknown_file_types(wig_list, transposons):
    """ """
    file_types = set(get_file_types(wig_list))
    method_types = set(transposons)
    extra_types = list(file_types - method_types)
    return extra_types


#

def get_data_zero_fill(wig_list):
    """ Returns a tuple of (data, position) containing a matrix of raw read counts,
        and list of coordinates. Positions that are missing are filled in as zero.

    Arguments:
        wig_list (list): List of paths to wig files.

    Returns:
        tuple: Two lists containing data and positions of the wig files given.
    """

    K = len(wig_list)
    T = 0

    if not wig_list:
        return (numpy.zeros((1, 0)), numpy.zeros(0), [])

    # NOTE:  This might be slow as we need to find the last insertion site
    #       over all the replicates. This might be an area to attempt to optimize.
    last_line = ""
    for wig_name in wig_list:
        for line in open(wig_name):
            if line[0] not in "0123456789":
                continue
            last_line = line
        pos = int(last_line.split()[0])
        T = max(T, pos)

    if T == 0:
        return (numpy.zeros((1, 0)), numpy.zeros(0), [])

    data = numpy.zeros((K, T))
    position = numpy.array(range(T)) + 1  # numpy.zeros(T)
    for j, path in enumerate(wig_list):
        reads = []
        i = 0
        for line in open(path):
            if line[0] not in "0123456789":
                continue
            tmp = line.split()
            pos = int(tmp[0])
            rd = float(tmp[1])
            prev_pos = pos
            data[j, pos - 1] = rd
            i += 1
    return (data, position)


def get_data_w_genome(wig_list, genome):

    X = read_genome(genome)
    N = len(X)
    positions = []
    pos2index = {}
    count = 0
    for i in range(N - 1):
        if X[i : i + 2].upper() == "TA":
            pos = i + 1
            positions.append(pos)
            pos2index[pos] = count
            count += 1

    positions = numpy.array(positions)
    T = len(positions)
    K = len(wig_list)
    data = numpy.zeros((K, T))
    for j, path in enumerate(wig_list):
        for line in open(path):
            if line[0] not in "0123456789":
                continue
            tmp = line.split()
            pos = int(tmp[0])
            rd = float(tmp[1])
            if pos in pos2index:
                index = pos2index[pos]
                data[j, index] = rd
            else:
                print(
                    "Warning: Coordinate %d did not match a TA site in the genome. Ignoring counts."
                    % (pos)
                )
    return (data, positions)


#


def combine_replicates(data, method="Sum"):
    """Returns list of data merged together.

    Arguments:
        data (list): List of numeric (replicate) data to be merged.
        method (str): How to combine the replicate dataset.

    Returns:
        list: List of numeric dataset now merged together.
    """

    if method == "Sum":
        combined = numpy.round(numpy.sum(data, 0))
    elif method == "Mean":
        combined = numpy.round(numpy.mean(data, 0))
    elif method == "TTRMean":
        # factors = norm_tools.TTR_factors(data)
        # data = factors * data
        (data, factors) = norm_tools.normalize_data(data, "TTR")
        target_factors = norm_tools.norm_to_target(data, 100)
        data = target_factors * data
        combined = numpy.round(numpy.mean(data, 0))
    else:
        combined = data[0, :]

    return combined


#


def get_data_stats(reads):
    density = numpy.mean(reads > 0)
    meanrd = numpy.mean(reads)
    nzmeanrd = numpy.mean(reads[reads > 0])
    nzmedianrd = numpy.median(reads[reads > 0])
    maxrd = numpy.max(reads)
    totalrd = numpy.sum(reads)

    skew = scipy.stats.skew(reads[reads > 0])
    kurtosis = scipy.stats.kurtosis(reads[reads > 0])
    return (density, meanrd, nzmeanrd, nzmedianrd, maxrd, totalrd, skew, kurtosis)


def get_wig_stats(path):
    """Returns statistics for the given wig file with read-counts.

    Arguments:
        path (str): String with the path to the wig file of interest.

    Returns:
        tuple: Tuple with the following statistical measures:
            - density
            - mean read
            - non-zero mean
            - non-zero median
            - max read
            - total reads
            - skew
            - kurtosis
    """
    (data, position) = CombinedWig.gather_wig_data([path])
    reads = data[0]
    return get_data_stats(reads)


#


def get_extended_pos_hash_pt(path, N=None):
    

    hash = {}
    maxcoord = float("-inf")
    data = []
    for line in open(path):
        if line.startswith("#"):
            continue
        tmp = line.split("\t")
        orf = tmp[8]
        start = int(tmp[1])
        end = int(tmp[2])
        maxcoord = max(maxcoord, start, end)
        data.append((orf, start, end))

    genome_start = 1
    if N:
        genome_end = maxcoord
    else:
        genome_end = N

    for i, (orf, start, end) in enumerate(data):
        if genome_start > start:
            genome_start = start

        prev_orf = ""
        if i > 0:
            prev_orf = data[i - 1][0]

        next_orf = ""
        if i < len(data) - 1:
            next_orf = data[i + 1][0]

        for pos in range(genome_start, end + 1):
            if pos not in hash:
                hash[pos] = {"current": [], "prev": [], "next": []}

            hash[pos]["prev"].append(prev_orf)

            if pos >= start:
                hash[pos]["next"].append(next_orf)
                hash[pos]["current"].append(orf)
            else:
                hash[pos]["next"].append(orf)
        genome_start = end + 1

    if N:
        for pos in range(maxcoord, genome_end + 1):
            if pos not in hash:
                hash[pos] = {"current": [], "prev": [], "next": []}
            hash[pos]["prev"].append(prev_orf)
    return hash


def get_extended_pos_hash_gff(path, N=None):
    

    hash = {}
    maxcoord = float("-inf")
    data = []
    for line in open(path):
        if line.startswith("#"):
            continue
        tmp = line.strip().split("\t")
        features = dict([tuple(f.split("=")) for f in tmp[8].split(";")])
        if "ID" not in features:
            continue
        orf = features["ID"]
        chr = tmp[0]
        type = tmp[2]
        start = int(tmp[3])
        end = int(tmp[4])
        maxcoord = max(maxcoord, start, end)
        data.append((orf, start, end))

    genome_start = 1
    if N:
        genome_end = maxcoord
    else:
        genome_end = N

    for i, (orf, start, end) in enumerate(data):

        if genome_start > start:
            genome_start = start

        prev_orf = ""
        if i > 0:
            prev_orf = data[i - 1][0]

        next_orf = ""
        if i < len(data) - 1:
            next_orf = data[i + 1][0]

        for pos in range(genome_start, end + 1):
            if pos not in hash:
                hash[pos] = {"current": [], "prev": [], "next": []}

            hash[pos]["prev"].append(prev_orf)

            if pos >= start:
                hash[pos]["next"].append(next_orf)
                hash[pos]["current"].append(orf)
            else:
                hash[pos]["next"].append(orf)
        genome_start = end + 1

    if N:
        for pos in range(maxcoord, genome_end + 1):
            if pos not in hash:
                hash[pos] = {"current": [], "prev": [], "next": []}
            hash[pos]["prev"].append(prev_orf)
    return hash


def get_pos_hash_pt(path):
    """Returns a dictionary that maps coordinates to a list of genes that occur at that coordinate.

    Arguments:
        path (str): Path to annotation in .prot_table format.

    Returns:
        dict: Dictionary of position to list of genes that share that position.
    """
    hash = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        tmp = line.strip().split("\t")
        orf = tmp[8]
        start = int(tmp[1])
        end = int(tmp[2])
        for pos in range(start, end + 1):
            if pos not in hash:
                hash[pos] = []
            hash[pos].append(orf)
    return hash


#


def get_pos_hash_gff(path):
    """Returns a dictionary that maps coordinates to a list of genes that occur at that coordinate.

    Arguments:
        path (str): Path to annotation in GFF3 format.

    Returns:
        dict: Dictionary of position to list of genes that share that position.
    """
    hash = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        tmp = line.strip().split("\t")
        features = dict(
            [
                tuple(f.split("=", 1))
                for f in filter(lambda x: "=" in x, tmp[8].split(";"))
            ]
        )
        if "ID" not in features:
            continue
        orf = features["ID"]
        chr = tmp[0]
        type = tmp[2]
        start = int(tmp[3])
        end = int(tmp[4])
        for pos in range(start, end + 1):
            if pos not in hash:
                hash[pos] = []
            hash[pos].append(orf)
    return hash


#


def get_pos_hash(path):
    """Returns a dictionary that maps coordinates to a list of genes that occur at that coordinate.

    Arguments:
        path (str): Path to annotation in .prot_table or GFF3 format.

    Returns:
        dict: Dictionary of position to list of genes that share that position.
    """
    filename, file_extension = os.path.splitext(path)
    if file_extension.lower() in [".gff", ".gff3"]:
        return get_pos_hash_gff(path)
    else:
        return get_pos_hash_pt(path)


#


def get_gene_info_pt(path):
    """Returns a dictionary that maps gene id to gene information.

    Arguments:
        path (str): Path to annotation in .prot_table format.

    Returns:
        dict: Dictionary of gene id to tuple of information:
            - name
            - description
            - start coordinate
            - end coordinate
            - strand

    """
    orf2info = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        tmp = line.strip().split("\t")
        orf = tmp[8]
        name = tmp[7]
        desc = tmp[0]
        start = int(tmp[1])
        end = int(tmp[2])
        strand = tmp[3]
        orf2info[orf] = (name, desc, start, end, strand)
    return orf2info


#


def get_gene_info_gff(path):
    """Returns a dictionary that maps gene id to gene information.

    Arguments:
        path (str): Path to annotation in GFF3 format.

    Returns:
        dict: Dictionary of gene id to tuple of information:
            - name
            - description
            - start coordinate
            - end coordinate
            - strand

    """
    orf2info = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        tmp = line.strip().split("\t")
        chr = tmp[0]
        type = tmp[2]
        start = int(tmp[3])
        end = int(tmp[4])
        length = ((end - start + 1) / 3) - 1
        strand = tmp[6]
        features = dict(
            [
                tuple(f.split("=", 1))
                for f in filter(lambda x: "=" in x, tmp[8].split(";"))
            ]
        )
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
        if desc == "-":
            desc = features.get("product", "-")

        orf2info[orf] = (name, desc, start, end, strand)
    return orf2info


#


def get_gene_info(path):
    """Returns a dictionary that maps gene id to gene information.

    Arguments:
        path (str): Path to annotation in .prot_table or GFF3 format.

    Returns:
        dict: Dictionary of gene id to tuple of information:
            - name
            - description
            - start coordinate
            - end coordinate
            - strand

    """
    filename, file_extension = os.path.splitext(path)
    if file_extension.lower() in [".gff", ".gff3"]:
        return get_gene_info_gff(path)
    else:
        return get_gene_info_pt(path)


#


def get_coordinate_map(galign_path, reverse=False):
    """Attempts to get mapping of coordinates from galign file.

    Arguments:
        path (str): Path to .galign file.
        reverse (bool): Boolean specifying whether to do A to B or B to A.

    Returns:
        dict: Dictionary of coordinate in one file to another file.
    """
    c1Toc2 = {}
    for line in open(galign_path):
        if line.startswith("#"):
            continue
        tmp = line.split()
        star = line.strip().endswith("*")
        leftempty = tmp[0].startswith("-")
        rightempty = tmp[1].endswith("-")
        if leftempty:
            left = -1
        else:
            left = int(tmp[0])
        if rightempty:
            right = -1
        elif leftempty:
            right = int(tmp[1])
        else:
            right = int(tmp[2])

        if not reverse:
            if not leftempty:
                c1Toc2[left] = right
        else:
            if not rightempty:
                c1Toc2[right] = left
    return c1Toc2


#


def read_genome(path):
    """Reads in FASTA formatted genome file.

    Arguments:
        path (str): Path to .galign file.

    Returns:
        string: String with the genomic sequence.
    """
    seq = ""
    for line in open(path):
        if line.startswith(">"):
            continue
        seq += line.strip()
    return seq


#


def maxrun(lst, item=0):
    """Returns the length of the maximum run an item in a given list.

    Arguments:
        lst (list): List of numeric items.
        item (float): Number to look for consecutive runs of.

    Returns:
        int: Length of the maximum run of consecutive instances of item.
    """
    best = 0
    i, n = 0, len(lst)
    while i < n:
        if lst[i] == item:
            j = i + 1
            while j < n and lst[j] == item:
                j += 1
            r = j - i
            if r > best:
                best = r
            i = j
        else:
            i += 1
    return best


#


def getR1(n):
    """Small Correction term. Defaults to 0.000016 for now"""
    return 0.000016


#


def getR2(n):
    """Small Correction term. Defaults to 0.00006 for now"""
    return 0.00006


#


def getE1(n):
    """Small Correction term. Defaults to 0.01 for now"""
    return 0.01


#


def getE2(n):
    """Small Correction term. Defaults to 0.01 for now"""
    return 0.01


#


def getGamma():
    """Euler-Mascheroni constant ~ 0.577215664901 """
    return 0.5772156649015328606


#


def ExpectedRuns(n, pnon):
    """Expected value of the run of non=insertions (Schilling, 1990):

        ER_n =  log(1/p)(nq) + gamma/ln(1/p) -1/2 + r1(n) + E1(n)

    Arguments:
        n (int): Integer representing the number of sites.
        pins (float): Floating point number representing the probability of non-insertion.

    Returns:
        float: Size of the expected maximum run.

    """
    if n < 20:  # use exact calculation for genes with less than 20 TA sites
        # Eqn 17-20 in Boyd, https://www.math.ubc.ca/~boyd/bern.runs/bernoulli.html
        #  recurrence relations for F(n,k) = prob that max run has length k
        p, q = 1 - pnon, pnon
        F = numpy.ones((n + 1, n + 1))
        for k in range(n):
            F[k + 1, k] = 1 - numpy.power(q, k + 1)
        for k in range(n + 1):
            for n in range(n + 1):
                if n >= k + 2:
                    F[n, k] = F[n - 1, k] - p * numpy.power(q, k + 1) * F[n - k - 2, k]
        ERn = 0
        for k in range(1, n + 1):
            ERn += k * (F[n, k] - F[n, k - 1])
        return ERn

    pins = 1 - pnon
    gamma = getGamma()
    r1 = getR1(n)
    E1 = getE1(n)
    A = math.log(n * pins, 1.0 / pnon)
    B = gamma / math.log(1.0 / pnon)
    ER = A + B - 0.5 + r1 + E1
    return ER


#


def VarR(n, pnon):
    """Variance of the expected run of non-insertons (Schilling, 1990):

    .. math::

        VarR_n =  (pi^2)/(6*ln(1/p)^2) + 1/12 + r2(n) + E2(n)


    Arguments:
        n (int): Integer representing the number of sites.
        pnon (float): Floating point number representing the probability of non-insertion.

    Returns:
        float: Variance of the length of the maximum run.
    """
    r2 = getR2(n)
    E2 = getE2(n)
    A = math.pow(math.pi, 2.0) / (6 * math.pow(math.log(1.0 / pnon), 2.0))
    V = A + 1 / 12.0 + r2 + E2
    return V


#


def GumbelCDF(x, u, B):
    """CDF of the Gumbel distribution:

        e^(-e^( (u-x)/B))

    Arguments:
        x (int): Length of the max run.
        u (float): Location parameter of the Gumbel dist.
        B (float): Scale parameter of the Gumbel dist.

    Returns:
        float: Cumulative probability o the Gumbel distribution.
    """
    return math.exp(-1 * math.exp((u - x) / B))


#


def griffin_analysis(genes_obj, pins):
    """Implements the basic Gumbel analysis of runs of non-insertion, described in Griffin et al. 2011.

    This analysis method calculates a p-value of observing the maximun run of
    TA sites without insertions in a row (i.e. a "run", r). Unusually long
    runs are indicative of an essential gene or protein domain. Assumes that
    there is a constant, global probability of observing an insertion
    (tantamount to a Bernoulli probability of success).

    Arguments:
        genes_obj (Genes): An object of the Genes class defining the genes.
        pins (float): The probability of insertion.

    Returns:
        list: List of lists with results and information for the genes. The elements of the list are as follows:
            - ORF ID.
            - Gene Name.
            - Gene Description.
            - Number of TA sites with insertions.
            - Number of TA sites.
            - Length of largest run of non-insertion.
            - Expected run for a gene this size.
            - p-value of the observed run.
    """

    pnon = 1.0 - pins
    results = []
    for gene in genes_obj:
        if gene.n == 0:
            results.append(
                [gene.orf, gene.name, gene.desc, gene.k, gene.n, gene.r, 0.0, 1.000]
            )
        else:
            # do I need to estimate B better (using exact calc for variance) for small genes too?
            B = 1.0 / math.log(
                1.0 / pnon
            )  # beta param of Gumbel distn; like tau in our Bioinfo paper
            # u = math.log(gene.n*pins, 1.0/pnon) # instead, calculate this based on estimate of ExpectedRun() length below
            exprun = ExpectedRuns(gene.n, pnon)
            # u is mu of Gumbel (mean=mu+gamma*beta); matching of moments; like Eq 5 in Schilling, but subtract off unneeded terms
            u = exprun - getGamma() / math.log(1.0 / pnon)
            pval = 1.0 - GumbelCDF(gene.r, u, B)
            results.append(
                [gene.orf, gene.name, gene.desc, gene.k, gene.n, gene.r, exprun, pval]
            )
    return results


#


def runs_w_info(data):
    """Return list of all the runs of consecutive non-insertions with the start and end locations.

    Arguments:
        data (list): List of numeric data to check for runs.

    Returns:
        list: List of dictionary from run to length and position information of the tun.
    """
    runs = []
    start = 1
    current_r = 0
    for read in data:
        if read > 0:  # If ending a run of zeros
            if current_r > 0:  # If we were in a run, add to list
                end = start + current_r - 1
                runs.append(dict(length=current_r, start=start, end=end))
            start = start + (current_r + 1)
            current_r = 0
        else:
            current_r += 1

    # If we ended in a run, add it
    if current_r > 0:
        end = start + current_r - 1
        runs.append(dict(length=current_r, start=start, end=end))
    return runs


#


def get_genes_in_range(pos_hash, start, end):
    """Returns list of genes that occur in a given range of coordinates.

    Arguments:
        pos_hash (dict): Dictionary of position to list of genes.
        start (int): Start coordinate of the desired range.
        end (int): End coordinate of the desired range.

    Returns:
        list: List of genes that fall within range.

    """

    genes = set()
    for pos in range(start, end + 1):
        if pos in pos_hash:
            genes.update(pos_hash[pos])

    return list(sorted(genes))


if __name__ == "__main__":

    G = Genes(sys.argv[1].split(","), sys.argv[2], norm="TTR")
    theta = G.global_theta()
    print("#Insertion: %s" % G.global_insertion())
    print("#Sites: %s" % G.global_sites())
    print("#Run: %s" % G.global_run())
    print("#Theta: %1.4f" % theta)
    print("#Phi: %1.4f" % G.global_phi())
    print("#")

    griffin_results = griffin_analysis(G, theta)
    for i, gene in enumerate(G):
        pos = gene.position
        exprun, pval = griffin_results[i][-2:]
        print(
            "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%1.1f\t%1.5f"
            % (
                gene.orf,
                gene.name,
                gene.k,
                gene.n,
                gene.r,
                gene.s,
                gene.t,
                exprun,
                pval,
            )
        )
