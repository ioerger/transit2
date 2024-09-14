import sys
import os
import math
import warnings
from functools import total_ordering
from collections import namedtuple
from os.path import isabs, isfile, isdir, join, dirname, basename, exists, splitext, relpath

import numpy
import scipy.stats
from pytransit.generic_tools import csv, misc
from pytransit.generic_tools.lazy_dict import LazyDict, stringify
from pytransit.generic_tools.named_list import named_list, NamedListBase
from pytransit.generic_tools.misc import line_count_of, flatten_once, no_duplicates, indent, cache
from pytransit.globals import logging
from super_hash import super_hash

try:
    from pytransit.specific_tools import norm_tools
    no_norm = False
except ImportError:
    no_norm = True
    warnings.warn("Problem importing the norm_tools.py module. Read-counts will not be normalized. Some functions may not work.")

class Wig:
    '''
        self.id # id from the metadata file
        self.fingerprint # the "File" column from the metadata 
        self.condition_names # a list of strings
        self.ta_sites # list of ints
        self.insertion_counts # list of numbers
        self.rows # each element is always [position_number, insertion_count]
        self.column_index # int (column inside combined wig)
        self.extra_data.count
        self.extra_data.sum
        self.extra_data.non_zero_mean
        self.extra_data.non_zero_median
        self.extra_data.density
        self.extra_data.mean
        self.extra_data.max
        self.extra_data.skew
        self.extra_data.kurtosis
    '''
    def __init__(self, *, rows=None, id=None, fingerprint=None, condition_names=tuple(), column_index=None, extra_data=None, path=None):
        from random import random
        self.rows            = rows or []
        self.comments        = []
        
        self.fingerprint     = fingerprint or path or self.path
        self.column_index    = column_index
        self.condition_names = condition_names
        self.id              = id or f"{basename(self.fingerprint)}_" + f"{super_hash((self.fingerprint, rows))}"[0:6]

        self.ta_sites         = [ each[0] for each in self.rows ]
        self.insertion_counts = [ each[1] for each in self.rows ]

        self.extra_data = LazyDict(extra_data or {})
        density, mean_reads, non_zero_mean_reads, non_zero_median_reads, max_reads, sum_of_reads, skew, kurtosis = get_data_stats(numpy.array(self.insertion_counts))
        self.extra_data.update(dict(
            count=len(self.rows),
            sum=sum_of_reads,
            non_zero_mean=non_zero_mean_reads,
            non_zero_median=non_zero_median_reads,
            density=density,
            mean=mean_reads,
            max=max_reads,
            skew=skew,
            kurtosis=kurtosis,
        ))
        
    def __repr__(self):
        return f"""Wig(
            fingerprint={self.fingerprint},
            column_index={self.column_index},
            condition_names={self.condition_names},
            rows_shape=({len(self.rows)}, {len(self.rows[0])}),
            extra_data={indent(self.extra_data, by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")
    
    def __hash__(self):
        return hash((self.fingerprint, self.column_index))
    
    def __eq__(self, other):
        return hash(other) == self.__hash__()

    @staticmethod
    def selected_as_gathered_data(wig_objects):
        """
            Returns:
                insertion_counts, ta_sites
        """
        import numpy
        from pytransit.globals import logging, gui, cli, root_folder, debugging_enabled
        
        # fail fast
        if len(wig_objects) == 0:
            return numpy.zeros((0,0)), numpy.array([])
        
        # get positions
        positions = wig_objects[0].ta_sites
        read_counts_per_wig = [ each_wig.insertion_counts for each_wig in wig_objects ]
        return numpy.array(read_counts_per_wig), numpy.array(positions)
    
    @staticmethod
    def read(path):
        comments, headers, rows = csv.read(path, seperator="\t", first_row_is_column_names=False, comment_symbol="#")
        rows = tuple( each_row for each_row in rows if len(each_row) != 0 )
        assert len(rows) > 1, f"it looks like the wig file is empty (no rows): {path}"
        if rows[0][0].startswith("variableStep"):
            rows = rows[1:]
        
        # if they used spaces instead of tabs, change it, then reload the file
        if len(rows[0]) == 1 and " " in f"""{rows[0][0]}""":
            with open(path,'r') as f:
                output = f.read()
            
            lines = output.split("\n")
            output = ""
            for each_line in lines:
                if not each_line.startswith("#"):
                    each_line = each_line.replace(" ", "\t")
                output += f"{each_line}\n"
            
            with open(path, 'w') as the_file:
                the_file.write(output)
            return Wig.read(path)
        
        return Wig(rows=rows, id=None, fingerprint=None, condition_names=tuple(), column_index=None, extra_data=None, path=path)

class Condition:
    def __init__(self, name, extra_data=None):
        self.name = name
        self.extra_data = LazyDict(extra_data or {})
    
    def __repr__(self):
        return f"""Condition(
            name={self.name},
            extra_data={indent(self.extra_data, by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")
    
    def __hash__(self):
        return hash(self.name)
    
    def __eq__(self, other):
        return self.name == other.name

class CombinedWigMetadata:
    special_headers = [ "Condition", "Filename", "Id" ]
    # can contain more data
    def __init__(self, path=None, rows=None, headers=None, comments=None):
        self.path     = path
        self.headers  = headers or []
        self.rows     = rows or []
        self.comments = comments or []
        self._conditions_by_wig_fingerprint        = None
        self._covariates_by_wig_fingerprint_list   = None
        self._interactions_by_wig_fingerprint_list = None
        self._ordering_metadata                    = None
        if path:
            self.comments, self.headers, self.rows = csv.read(self.path, seperator="\t", first_row_is_column_names=True, comment_symbol="#")
            if any([ each_special_header not in self.headers for each_special_header in CombinedWigMetadata.special_headers]):
                logging.error(f'''
                    For metadataa file: {path}, I expected (at least) these headers: {CombinedWigMetadata.special_headers}, but got these headers: {self.headers}
                ''')
        
        self.covariate_names = [ each for each in self.headers if each not in CombinedWigMetadata.special_headers ]
        if not path:
            # 
            # generate ordering_metadata, conditions_by_wig_fingerprint
            # 
            self._covariates_by_wig_fingerprint_list   = []
            self._interactions_by_wig_fingerprint_list = []
            self._conditions_by_wig_fingerprint        = {}
            self._ordering_metadata                    = { "condition": [] }
            for row in self.rows:
                wig_fingerprint = row["Filename"]
                condition_name  = row["Condition"]
                self._conditions_by_wig_fingerprint[wig_fingerprint] = condition_name
                self._ordering_metadata["condition"].append(condition_name)
        
        self.conditions = no_duplicates(
            Condition(
                name=each_row["Condition"],
            )
                for each_row in self.rows
        )
    
    def copy(self):
        from copy import deepcopy
        new_rows = []
        for each_row in self.rows:
            new_rows.append(deepcopy(each_row))
        return CombinedWigMetadata(rows=new_rows, headers=self.headers, comments=self.comments)
        
    # TODO: maybe rename to "select_by"
    def with_only(self, condition_names=None, wig_fingerprints=None):
        from copy import deepcopy
        new_rows = []
        for each_row in self.rows:
            if condition_names and each_row["Condition"] not in condition_names:
                continue
            if wig_fingerprints and each_row["Filename"] not in wig_fingerprints:
                continue
            new_rows.append(deepcopy(each_row))
        
        return CombinedWigMetadata(rows=new_rows, headers=self.headers, comments=self.comments)
    
    def column(self, column_name):
        return [ each[column_name] for each in rows ]
    
    @property
    def condition_names(self, *, wig_fingerprint=None, id=None):
        return no_duplicates([ each_row["Condition"] for each_row in self.rows ])
    
    def condition_names_for(self, *, wig_fingerprint=None, id=None):
        conditions = []
        if wig_fingerprint:
            for each_row in self.rows:
                if each_row["Filename"] == wig_fingerprint:
                    conditions.append(each_row["Condition"])
        elif id:
            for each_row in self.rows:
                if each_row["Id"] == id:
                    conditions.append(each_row["Condition"])
        return no_duplicates(conditions)
    
    def row_for(self, wig_fingerprint=None):
        for each_row in self.rows:
            if each_row["Filename"] == wig_fingerprint:
                return each_row
    
    def covariates_dict_for(self, wig_fingerprint=None):
        row = self.row_for(wig_fingerprint=wig_fingerprint)
        if row:
            return {
                covariate_name : row[covariate_name]
                    for covariate_name in self.covariate_names
            }
    
    def id_for(self, wig_fingerprint=None):
        return self.row_for(wig_fingerprint=wig_fingerprint)["Id"]
    
    def fingerprint_for(self, wig_id=None):
        for each_row in self.rows:
            if each_row["Id"] == wig_id:
                return each_row["Filename"]
    
    def fingerprints_for(self, condition):
        from pytransit.generic_tools import misc
        return misc.no_duplicates([
            each_row["Filename"]
                for each_row in self.rows
                    if each_row["Condition"] == condition
        ])
    
    @property
    def wig_fingerprints(self):
        from pytransit.generic_tools import misc
        return misc.no_duplicates([ each_row["Filename"] for each_row in self.rows ])
    
    @property
    def wig_ids(self):
        from pytransit.generic_tools import misc
        return misc.no_duplicates([ each_row["Id"] for each_row in self.rows ])
    
    def __repr__(self):
        return f"""CWigMetadata(
            path={self.path},
            rows_shape=({len(self.rows)}, {len(self.headers)}),
            conditions={indent(self.conditions, by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")
    
    @staticmethod
    def read_condition_data(path, covars_to_read=[], interactions_to_read=[], column_name_for_condition="Condition"):
        conditions_by_wig_fingerprint        = {}
        covariates_by_wig_fingerprint_list   = [{} for _ in covars_to_read]
        interactions_by_wig_fingerprint_list = [{} for _ in interactions_to_read]
        ordering_metadata         = {"condition": [], "interaction": []}
        covars_to_read = [ each.lower() for each in covars_to_read ]
        interactions_to_read = [ each.lower() for each in interactions_to_read ]
        with open(path) as file:
            lines = file.readlines()
            column_names              = lines[0].split()
            column_names              = [ each.lower() for each in column_names ]
            try:
                index_for_condition       = column_names.index(column_name_for_condition.lower())
                index_for_wig_fingerprint = column_names.index("Filename".lower())
            except Exception as error:
                print(f'''column_names = {column_names}''')
                print(f'''path = {path}''')
                raise error
            # validate
            for each in covars_to_read:
                if each not in column_names:
                    raise Exception(f'''Tried to select {each} as a covariate, but the available covariates are: {column_names}''')
            for each in interactions_to_read:
                if each not in column_names:
                    raise Exception(f'''Tried to select {each} as an interaction, but the available interactions are: {column_names}''')
            covar_indexes             = [ column_names.index(each.lower()) for each in covars_to_read       ]
            interaction_indexes       = [ column_names.index(each.lower()) for each in interactions_to_read ]

            for each_line in lines[1:]:
                # skip comments
                if each_line[0] == "#":
                    continue
                row_as_list = each_line.rstrip().split('\t') # allow spaces in Filenames
                condition       = row_as_list[index_for_condition]
                wig_fingerprint = row_as_list[index_for_wig_fingerprint]
                
                conditions_by_wig_fingerprint[wig_fingerprint] = condition
                ordering_metadata["condition"].append(condition)
                for index, each_covar_column_index in enumerate(covar_indexes):
                    covariates_by_wig_fingerprint_list[index][wig_fingerprint] = row_as_list[each_covar_column_index]
                for index, each_interaction_column_index in enumerate(interaction_indexes):
                    interactions_by_wig_fingerprint_list[index][wig_fingerprint] = row_as_list[each_interaction_column_index]
                    
                    # TODO
                    # This makes sense only if there is only 1 interaction variable
                    # For multiple interaction vars, may have to rethink ordering.
                    ordering_metadata["interaction"].append(row_as_list[interaction_indexes[index]])

        return (
            conditions_by_wig_fingerprint,
            covariates_by_wig_fingerprint_list,
            interactions_by_wig_fingerprint_list,
            ordering_metadata,
        )
    
    
    # 
    # conditions_by_wig_fingerprint
    # 
    @property
    def conditions_by_wig_fingerprint(self): 
        if type(self._conditions_by_wig_fingerprint) == type(None):
            self._conditions_by_wig_fingerprint, self._covariates_by_wig_fingerprint_list, self._interactions_by_wig_fingerprint_list, self._ordering_metadata = self.read_condition_data(path=self.path)
        return self._conditions_by_wig_fingerprint
            
    # 
    # covariates_by_wig_fingerprint_list
    # 
    @property
    def covariates_by_wig_fingerprint_list(self): 
        if type(self._covariates_by_wig_fingerprint_list) == type(None):
            self._conditions_by_wig_fingerprint, self._covariates_by_wig_fingerprint_list, self._interactions_by_wig_fingerprint_list, self._ordering_metadata = self.read_condition_data(path=self.path)
        return self._covariates_by_wig_fingerprint_list
            
    # 
    # interactions_by_wig_fingerprint_list
    # 
    @property
    def interactions_by_wig_fingerprint_list(self): 
        if type(self._interactions_by_wig_fingerprint_list) == type(None):
            self._conditions_by_wig_fingerprint, self._covariates_by_wig_fingerprint_list, self._interactions_by_wig_fingerprint_list, self._ordering_metadata = self.read_condition_data(path=self.path)
        return self._interactions_by_wig_fingerprint_list
            
    # 
    # ordering_metadata
    # 
    @property
    def ordering_metadata(self): 
        if type(self._ordering_metadata) == type(None):
            self._conditions_by_wig_fingerprint, self._covariates_by_wig_fingerprint_list, self._interactions_by_wig_fingerprint_list, self._ordering_metadata = self.read_condition_data(path=self.path)
        return self._ordering_metadata
            
class CombinedWigData(NamedListBase):
    _names_to_index = {'sites':0,'counts_by_wig':1, 'wig_fingerprints':2 }
    @staticmethod
    @cache(watch_filepaths=lambda file_path: [file_path])
    def load(file_path):
        """
            Read the combined wig-file generated by Transit
            :: Filename -> Tuple([Site], [WigData], [Filename])
            Site :: Integer
            WigData :: [Number]
            Filename :: String
        """
        import ez_yaml
        from pytransit.specific_tools import transit_tools
        
        sites, counts_by_wig, wig_fingerprints, extra_data = [], [], [], {}
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
                        yaml_string += f"\n{line[1:]}"
                        continue
                    else:
                        yaml_mode_is_on = False
                        # add to the extra_data dict when its done
                        if len(yaml_string) > 0:
                            an_object = ez_yaml.to_object(string=yaml_string)
                            extra_data.update(an_object["extra_data"] or {})
                            if len(wig_fingerprints) == 0:
                                wig_fingerprints += extra_data.get('wig_fingerprints',[])
                    # 
                    # handle older file method
                    # 
                    old_header = "#File: "
                    if line.startswith(old_header) or line.startswith(old_header.lower()):
                        wig_fingerprints.append(line.rstrip()[len(old_header):])  # allows for spaces in filenames
                        continue
            
            # 
            # handle body
            # 
            if len(set(wig_fingerprints)) != len(wig_fingerprints):
                logging.warn(f"{file_path} contains duplicate file entries")
            counts_by_wig = [ [] for _ in wig_fingerprints ]
            for index, line in enumerate(lines):
                if index % 150 == 0: # 150 is arbitrary, bigger = slower visual update but faster read
                    percent_done = round(index/number_of_lines * 100, ndigits=2)
                    logging.log(f"\rreading lines: {percent_done}%          ", end="")
                
                # lines to skip
                if line.startswith("#") or len(line) == 0:
                    continue
                
                #
                # actual parsing
                #
                cols = line.split("\t")[0 : 1+len(wig_fingerprints)]
                cols = cols[: 1+len(wig_fingerprints)]  # additional columns at end could contain gene info
                # Read in position as int, and readcounts as float
                try:
                    cols = [
                        int(each_t_iv) if index == 0 else float(each_t_iv)
                            for index, each_t_iv in enumerate(cols)
                    ]
                except Exception as error:
                    print(f'''cols = {cols}''')
                    print(f'''len(cols) = {len(cols)}''')
                    print(f'''wig_fingerprints = {wig_fingerprints}''')
                    print(f'''len(wig_fingerprints) = {len(wig_fingerprints)}''')
                    import sys
                    sys.exit()
                position, wig_counts = cols[0], cols[1:]
                sites.append(position)
                for index, count in enumerate(wig_counts):
                    counts_by_wig[index].append(count)
            logging.log(f"\rreading lines: 100%          ")
        
        output_object = CombinedWigData((numpy.array(sites), numpy.array(counts_by_wig), wig_fingerprints))
        logging.log(f"created combined wig data object")
        
        return output_object

from pytransit.generic_tools.named_list import named_list
class CombinedWig:
    '''
        self.as_tuple         # (numpy.array(sites), numpy.array(counts_by_wig), wig_fingerprints)
        self.rows             # equivalent to the CSV rows of .comwig file; a list of lists, can contain numbers and strings
        self.wig_ids          # same order as columns/wig_fingerprints
        self.wig_fingerprints # same order as #File: columns
        self.read_counts_array[row_index, wig_index]
        self.main_path
        self.metadata_path # to get all these it would be [ each.metadata_path for each in gui.combined_wigs ]
        self.samples # list of Wig objects
        self.metadata # CombinedWigMetadata object
        self.metadata.path
        self.metadata.headers
        self.metadata.rows
        self.metadata.conditions
        self.metadata.condition_names
        self.metadata.wig_ids
        self.metadata.wig_fingerprints
        self.metadata.with_only(condition_names=[], wig_fingerprints=[])
        self.metadata.condition_names_for(wig_fingerprint="")
        self.metadata.condition_names_for(wig_id="")
        self.metadata.id_for(wig_fingerprint)
        self.metadata.fingerprints_for(condition_name)
    '''
    PositionsAndReads = named_list(["read_counts", "positions"])
    
    @staticmethod
    @cache(watch_filepaths=lambda *, main_path, metadata_path=None, annotation_path=None, comments=None, extra_data=None: [ main_path, metadata_path, annotation_path ])
    def load(*args, **kwargs):
        return CombinedWig(*args, **kwargs)
    
    def __init__(self, *, main_path, metadata_path=None, annotation_path=None, comments=None, extra_data=None):
        self.main_path       = main_path
        self.metadata_path   = metadata_path
        self.annotation_path = annotation_path
        self.metadata        = CombinedWigMetadata(self.metadata_path) if self.metadata_path else None
        self.as_tuple        = CombinedWigData.load(self.main_path) # for backwards compatibility (otherwise just used self.rows and helper methods)
        self.rows            = []
        self.comments        = comments or []
        self.extra_data      = LazyDict(extra_data or {})
        self.samples         = []
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
    # positions
    # 
    @property
    def ta_sites(self):
        return self.as_tuple.sites
    
    # 
    # files (same order as columns)
    # 
    @property
    def wig_fingerprints(self): return self.extra_data["wig_fingerprints"]
    @wig_fingerprints.setter
    def wig_fingerprints(self, value):
        self.extra_data["wig_fingerprints"] = value
    # 
    # wig_ids (same order as columns/wig_fingerprints)
    # 
    @property
    def wig_ids(self):
        return [ self.metadata.id_for(wig_fingerprint=each) for each in self.wig_fingerprints]
    
    # 
    # conditions
    # 
    @property
    def conditions(self):
        return self.metadata.conditions
    @property
    def condition_names(self):
        return self.metadata.condition_names
    
    # 
    # read_counts_array
    # 
    @property
    def read_counts_array(self):
        return self.as_tuple.counts_by_wig.transpose()
    
    # 
    # read_counts_wig
    # 
    @property
    @cache(watch_attributes=lambda self:[ self.wig_fingerprints ])
    def read_counts_by_wig_fingerprint(self):
        counts_for_wig = { each_path: [] for each_path in self.wig_fingerprints }
        for each_row in self.rows:
            for each_wig_fingerprint in self.wig_fingerprints:
                counts_for_wig[each_wig_fingerprint].append(
                    each_row[each_wig_fingerprint]
                )
        return counts_for_wig
    
    def copy(self):
        from copy import deepcopy
        # create a new shell
        class CombinedWigHelper(CombinedWig):
            def __init__(*args): pass # only difference it disabling the init 
        
        new_combined_wig = CombinedWigHelper()
        
        # copy the easy stuff
        new_combined_wig.main_path       = None
        new_combined_wig.metadata_path   = None
        new_combined_wig.annotation_path = self.annotation_path
        new_combined_wig.comments        = list(self.comments)
        new_combined_wig.extra_data      = deepcopy(dict(self.extra_data))
        new_combined_wig.rows            = []
        new_combined_wig.samples         = []
        new_combined_wig.metadata        = self.metadata and self.metadata.copy()
        
        new_combined_wig.as_tuple = CombinedWigData(deepcopy(each) for each in self.as_tuple)
        new_combined_wig.CWigRow = named_list([ "position", *new_combined_wig.wig_fingerprints ])
        for each_row in self.rows:
            new_combined_wig.rows.append(new_combined_wig.CWigRow(each_row))
        
        read_counts_by_wig_fingerprint = new_combined_wig.read_counts_by_wig_fingerprint
        for column_index, wig_fingerprint in enumerate(new_combined_wig.wig_fingerprints):
            new_combined_wig.samples.append(
                Wig(
                    rows=list(zip(new_combined_wig.ta_sites, read_counts_by_wig_fingerprint[wig_fingerprint])),
                    id=new_combined_wig.metadata and new_combined_wig.metadata.id_for(wig_fingerprint=wig_fingerprint),
                    fingerprint=wig_fingerprint,
                    column_index=column_index,
                    condition_names=new_combined_wig.metadata and new_combined_wig.metadata.condition_names_for(wig_fingerprint=wig_fingerprint),
                    extra_data=LazyDict(
                        is_part_of_cwig=True,
                    ),
                )
            )
        
        return new_combined_wig
        
        
    def with_only(self, condition_names=None, wig_fingerprints=None, wig_ids=None):
        import numpy
        from copy import deepcopy
        
        assert wig_ids == None or self.metadata != None, "Tried to filter by wig_ids, but the combined wig didnt have metadata attached (e.g. IDK what id's exist)"
        
        included_wig_fingerprints = self.wig_fingerprints
        if wig_ids != None:
            fingerprints_from_ids = [ self.metadata.fingerprint_for(each_id) for each_id in wig_ids ]
            included_wig_fingerprints = [ each for each in included_wig_fingerprints if each in fingerprints_from_ids ]
        if wig_fingerprints != None:
            included_wig_fingerprints = [ each for each in included_wig_fingerprints if each in wig_fingerprints ]
        if condition_names != None:
            fingerprints_from_conditions = flatten_once([ self.metadata.fingerprints_for(each_name) for each_name in condition_names ])
            included_wig_fingerprints = [ each for each in included_wig_fingerprints if each in fingerprints_from_conditions ]
            
        
        new_wig_fingerprints = included_wig_fingerprints
        
        # create a new shell
        class CombinedWigHelper(CombinedWig):
            def __init__(*args): pass # only difference it disabling the init 
        
        new_combined_wig = CombinedWigHelper()
        
        # copy the easy stuff
        new_combined_wig.main_path       = None
        new_combined_wig.metadata_path   = None
        new_combined_wig.annotation_path = self.annotation_path
        new_combined_wig.comments        = self.comments
        new_combined_wig.extra_data      = deepcopy(dict(self.extra_data))
        new_combined_wig.rows            = []
        new_combined_wig.samples         = []
        new_combined_wig.metadata        = self.metadata and self.metadata.with_only(condition_names=condition_names, wig_fingerprints=new_wig_fingerprints)
        
        # extract only data relating to new_wig_fingerprints
        sites, counts_by_wig_array, old_wig_fingerprints = self.as_tuple
        new_counts_by_wig = numpy.zeros((len(new_wig_fingerprints), sites.shape[0], ))
        for index, each_fingerprint in enumerate(new_wig_fingerprints):
            old_index = old_wig_fingerprints.index(each_fingerprint)
            new_counts_by_wig[index,:] = counts_by_wig_array[old_index,:]
        
        removed_wig_fingerprints = set(old_wig_fingerprints) - set(new_wig_fingerprints)
        # create new data object
        new_combined_wig.as_tuple = CombinedWigData((sites, new_counts_by_wig, new_wig_fingerprints))
        
        # generate named lists for rows
        new_combined_wig.extra_data["wig_fingerprints"] = new_wig_fingerprints
        new_combined_wig.CWigRow = named_list([ "position", *new_wig_fingerprints ])
        for position, row in zip(sites, new_counts_by_wig.transpose()):
            #
            # extract
            #
            read_counts = row[ 0: len(new_wig_fingerprints) ]
            other      = row[ len(new_wig_fingerprints) :  ] # often empty
            
            # force types
            position   = int(position)
            read_counts = [ float(each) for each in read_counts ]
            
            # save
            new_combined_wig.rows.append(new_combined_wig.CWigRow([position]+read_counts+[each for each in other]))
            
        for column_index, wig_fingerprint in enumerate(new_combined_wig.wig_fingerprints):
            new_combined_wig.samples.append(
                Wig(
                    rows=list(zip(new_combined_wig.ta_sites, new_combined_wig.read_counts_by_wig_fingerprint[wig_fingerprint])),
                    id=new_combined_wig.metadata and new_combined_wig.metadata.id_for(wig_fingerprint=wig_fingerprint),
                    fingerprint=wig_fingerprint,
                    column_index=column_index,
                    condition_names=new_combined_wig.metadata and new_combined_wig.metadata.condition_names_for(wig_fingerprint=wig_fingerprint),
                    extra_data=LazyDict(
                        is_part_of_cwig=True,
                    ),
                )
            )
        
        # Helpful check to prevent far down-the-line errors
        if len(new_combined_wig.as_tuple.counts_by_wig) == 0:
            if condition_names  == None: condition_names = "*All*"
            if wig_fingerprints == None: wig_fingerprints = "*All*"
            if wig_ids          == None: wig_ids = "*All*"
            logging.error(f"""
                when calling .with_only() these were the restrictions:
                    condition_names={condition_names}
                    wig_fingerprints={wig_fingerprints}
                    wig_ids={wig_ids}
                However these are the available values
                    condition_names={self.metadata.condition_names}
                    wig_fingerprints={self.wig_fingerprints}
                    wig_ids={self.metadata.wig_ids}
                So after performing a .with_only()
                there are no samples in the resulting combined_wig
            """)
        return new_combined_wig
        
    def summed(self, *, by_conditions=False):
        if by_conditions:
            assert self.metadata != None, "Tried to sum by conditions, but the combined wig didnt have metadata attached (e.g. IDK what conditions exist)"
            new_wigs = []
            new_read_counts = []
            condition_names = tuple(self.metadata.condition_names)
            for each_condition_name in condition_names:
                combined_reads_for_numpy = []
                
                # 
                # find all the related wigs (overlaps are possible)
                # 
                for column_index, wig_fingerprint in enumerate(self.wig_fingerprints):
                    if each_condition_name in self.metadata.condition_names_for(wig_fingerprint=wig_fingerprint):
                        combined_reads_for_numpy.append(self.read_counts_by_wig_fingerprint[wig_fingerprint])
                # 
                # Create new Wig objects and read counts
                # 
                import numpy
                if len(combined_reads_for_numpy) > 0:
                    read_counts = numpy.array(combined_reads_for_numpy).sum(axis=0).tolist()
                    new_read_counts.append(read_counts)
                    new_wigs.append(
                        Wig(
                            rows=list(zip(self.ta_sites, read_counts)),
                            id=each_condition_name,
                            fingerprint=each_condition_name,
                            column_index=len(new_wigs),
                            condition_names=[each_condition_name],
                            extra_data=LazyDict(
                                is_part_of_cwig=True,
                            ),
                        )
                    )
            
            copy_of_self = self.with_only()
            copy_of_self.as_tuple = CombinedWigData((self.as_tuple.sites, numpy.array(new_read_counts), condition_names))
            copy_of_self.samples = new_wigs
            copy_of_self.rows = numpy.hstack(( self.as_tuple.sites, copy_of_self.as_tuple.counts_by_wig.transpose() ))
            copy_of_self.metadata = CombinedWigMetadata(
                headers=headers,
                rows=[
                    { "Id": each_name, "Condition": each_name, "Filename": each_name }
                        for each_name in condition_names
                ],
                comments=self.metadata.comments,
            )
            
            return copy_of_self
        
        return self
    
    def averaged(self, by_genes=False, by_conditions=False, n_terminus=0, c_terminus=0):
        # 
        # averaged by genes only
        # 
        if by_genes and not by_conditions:
            genes = tnseq_tools.Genes(
                wig_list=[],
                data=self.as_tuple.counts_by_wig,
                annotation=self.annotation_path,
                position=self.as_tuple.sites,
                n_terminus=n_terminus,
                c_terminus=c_terminus,
            )
            new_read_counts = numpy.vstack(each.reads for each in genes)
            
            # wrap inside a CombinedWig to regain helper methods
            copy_of_self = self.with_only()
            copy_of_self.as_tuple = CombinedWigData((self.as_tuple.sites, numpy.array(new_read_counts), condition_names))
            copy_of_self.is_by_genes = True
            # TODO: copy_of_self.samples = new_wigs
            return copy_of_self
        
        # 
        # averaged by both
        # 
        # both is not as straightforward as applying both operations because of the Simpsons Paradox
        if by_conditions and by_genes: 
            condition_names = tuple(self.metadata.condition_names)
            count_per_condition = [ len(self.metadata.fingerprints_for(each_name)) for each_name in condition_names ]
            
            summed_by_condition = self.summed(by_conditions=True)
            genes = tnseq_tools.Genes(
                wig_list=[],
                data=summed_by_condition.as_tuple.counts_by_wig,
                annotation=summed_by_condition.annotation_path,
                position=summed_by_condition.as_tuple.sites,
                n_terminus=n_terminus,
                c_terminus=c_terminus,
            )
            reads_of_first_gene = numpy.array(next(iter(genes)).reads)
            assert reads_of_first_gene.shape[1] == len(condition_names), f"columns=condition names but, {len(condition_names)} != 2nd dim of: {reads_of_first_gene.shape}"
            
            new_read_counts = numpy.vstack(numpy.array(each.reads).mean(axis=1) for each in genes)
            # divide the sums so that it becomes an average
            for condition_index, number_of_samples_in_condition in enumerate(count_per_condition):
                new_read_counts[:,condition_index] = new_read_counts[:,condition_index] / number_of_samples_in_condition
            
            copy_of_self = self.with_only()
            copy_of_self.is_by_genes = True
            copy_of_self.as_tuple = CombinedWigData((self.as_tuple.sites, numpy.array(new_read_counts), condition_names))
            # TODO: copy_of_self.samples = new_wigs
            copy_of_self.metadata = CombinedWigMetadata(
                headers=headers,
                rows=[
                    { "Id": each_name, "Condition": each_name, "Filename": each_name }
                        for each_name in condition_names
                ],
                comments=self.metadata.comments,
            )
        
        # 
        # by conditions only
        # 
        if by_conditions and not by_genes: # redundant if statement check for clarity     
            new_wigs = []
            new_read_counts = []
            condition_names = tuple(self.metadata.condition_names)
            for each_condition_name in condition_names:
                combined_reads_for_numpy = []
                
                # 
                # find all the related wigs (overlaps are possible)
                # 
                for column_index, wig_fingerprint in enumerate(self.wig_fingerprints):
                    if each_condition_name in self.metadata.condition_names_for(wig_fingerprint=wig_fingerprint):
                        combined_reads_for_numpy.append(self.read_counts_by_wig_fingerprint[wig_fingerprint])
                # 
                # Create new Wig objects and read counts
                # 
                import numpy
                if len(combined_reads_for_numpy) > 0:
                    read_counts = numpy.array(combined_reads_for_numpy).sum(axis=0).tolist()
                    new_read_counts.append(read_counts)
                    new_wigs.append(
                        Wig(
                            rows=list(zip(self.ta_sites, read_counts)),
                            id=each_condition_name,
                            fingerprint=each_condition_name,
                            column_index=len(new_wigs),
                            condition_names=[each_condition_name],
                            extra_data=LazyDict(
                                is_part_of_cwig=True,
                            ),
                        )
                    )
            
            copy_of_self = self.with_only()
            copy_of_self.as_tuple = CombinedWigData((self.as_tuple.sites, numpy.array(new_read_counts), condition_names))
            copy_of_self.samples = new_wigs
            copy_of_self.rows = numpy.hstack(( self.as_tuple.sites, copy_of_self.as_tuple.counts_by_wig.transpose() ))
            copy_of_self.metadata = CombinedWigMetadata(
                headers=["Id", "Condition", "Filename"],
                comments=self.metadata.comments,
                rows=[
                    { "Id": each_name, "Condition": each_name, "Filename": each_name }
                        for each_name in condition_names
                ],
            )
            
            return copy_of_self
        
        return self
    
    def normalized_with(self, kind="TTR"):
        from pytransit.specific_tools import norm_tools
        new_combined_wig = self.copy()
        (read_counts_per_wig, factors) = norm_tools.normalize_data(
            new_combined_wig.as_tuple.counts_by_wig,
            kind,
            self.wig_fingerprints,
            self.annotation_path,
        )
        new_combined_wig.as_tuple = CombinedWigData((self.as_tuple.sites, numpy.array(read_counts_per_wig), self.as_tuple.wig_fingerprints))
        return new_combined_wig
    
    def with_loess_correction(self):
        from pytransit.specific_tools import stat_tools
        new_combined_wig = self.copy()
        for wig_index, _ in range(len(new_combined_wig.as_tuple.counts_by_wig)):
            new_combined_wig.as_tuple.counts_by_wig[wig_index] = stat_tools.loess_correction(
                new_combined_wig.ta_sites,
                new_combined_wig.as_tuple.counts_by_wig[wig_index]
            )
        
        return new_combined_wig
    
    def get_genes(self, ignore_codon=True, n_terminus=0.0, c_terminus=0.0, include_nc=False, norm="nonorm", reps="All", minread=1, genome="", transposon="himar1"):
        return Genes(
            self.wig_fingerprints,
            self.annotation_path,
            ignore_codon=ignore_codon,
            n_terminus=n_terminus,
            c_terminus=c_terminus,
            data=self.as_tuple.counts_by_wig,
            position=self.as_tuple.sites,
            include_nc=include_nc,
            norm=norm,
            reps=reps,
            minread=minread,
            genome=genome,
            transposon=transposon,
        )
    
    def wig_with_id(self, id):
        for each in self.samples:
            if each.id == id:
                return each
    
    def _load_main_path(self):
        from pytransit.generic_tools.informative_iterator import ProgressBar
        import ez_yaml
        comments, headers, rows = csv.read(self.main_path, seperator="\t", first_row_is_column_names=False, comment_symbol="#")
        comment_string = "\n".join(comments)

        sites, wig_fingerprints, extra_data = [], [], {}
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
                wig_fingerprints.append(line.rstrip()[6:])  # allows for spaces in filenames
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
            wig_fingerprints += self.extra_data.get('wig_fingerprints',[])
            no_duplicates = []
            # remove any duplicate entries while preserving order
            for each_fingerpint in wig_fingerprints:
                if each_fingerpint not in no_duplicates:
                    no_duplicates.append(each_fingerpint)
            wig_fingerprints = no_duplicates
        
        # 
        # process data
        # 
        extra_data["wig_fingerprints"] = wig_fingerprints
        self.extra_data.update(extra_data)
        self.CWigRow = named_list([ "position", *wig_fingerprints ])
        for progress, row in ProgressBar(rows, disable_logging=True):
            #
            # extract
            #
            position   = row[0]
            read_counts = row[ 1: 1+len(wig_fingerprints) ]
            other      = row[ 1+len(wig_fingerprints) :  ] # often empty
            
            # force types
            position   = int(position)
            read_counts = [ float(each) for each in read_counts ]
            
            # save
            self.rows.append(self.CWigRow([position]+read_counts+other))
            
            sites.append(position)

            if progress.updated:
                logging.log(f"formatting data {progress.percent:.0f}%", end="\r")

        logging.log(f"formatting data 100%", end="\r")
        
        for column_index, wig_fingerprint in enumerate(self.wig_fingerprints):
            self.samples.append(
                Wig(
                    rows=list(zip(self.ta_sites, self.read_counts_by_wig_fingerprint[wig_fingerprint])),
                    id=self.metadata.id_for(wig_fingerprint=wig_fingerprint) if self.metadata else None,
                    fingerprint=wig_fingerprint,
                    column_index=column_index,
                    condition_names=self.metadata.condition_names_for(wig_fingerprint=wig_fingerprint) if self.metadata else None,
                    extra_data=LazyDict(
                        is_part_of_cwig=True,
                    ),
                )
            )
        
        return self
    
    @staticmethod
    def gather_wig_data(list_of_paths):
        """ 
    
        Arguments:
            wig_list (list): List of paths to wig files.

        Returns:
            (read_counts, ta_site_positions)

        :Example:

            >>> from pytransit.specific_tools.tnseq_tools import CombinedWig
            >>> (data, position) = CombinedWig.gather_wig_data(["data/cholesterol_glycerol.transit/glycerol_rep1.wig", "data/cholesterol_glycerol.transit/glycerol_rep2.wig"])
            >>> print(data)
            array([[ 0.,  0.,  0., ...,  0.,  0.,  0.],
                [ 0.,  0.,  0., ...,  0.,  0.,  0.]])

        .. seealso:: :class:`get_file_types` :class:`combine_replicates` :class:`get_data_zero_fill` :class:`pytransit.specific_tools.norm_tools.normalize_data`
        """
        # if its already been processed, just return it
        if isinstance(list_of_paths, CombinedWig.PositionsAndReads):
            return list_of_paths
            
        # If empty just quickly return empty lists
        if not list_of_paths:
            return CombinedWig.PositionsAndReads((numpy.zeros((1, 0)), numpy.zeros(0)))

        # Check size of all wig file matches
        size_list = []
        for path in list_of_paths:
            line_count = 0
            with open(path) as file:
                for line in file:
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
        position_per_line = numpy.zeros(line_count, dtype=int)
        for path_index, path in enumerate(list_of_paths):
            line_index = 0
            prev_pos   = 0
            with open(path) as file:
                for line in file:
                    if line[0] not in "0123456789": # not sure why this is here -- Jeff
                        continue
                    
                    tmp           = line.split()
                    each_position = int(tmp[0])
                    read_count    = float(tmp[1])
                    prev_pos      = each_position

                    try:
                        data_per_path[path_index, line_index] = read_count
                    except Exception as error:
                        raise Exception(f'''
                            
                            Make sure that all wig files have the same number of TA sites (i.e. same strain)
                            
                            Original Error:\n{error}
                            
                        ''')
                    position_per_line[line_index] = each_position
                    line_index += 1
        
        return CombinedWig.PositionsAndReads((data_per_path, position_per_line))
    
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

        >>> from pytransit.specific_tools import tnseq_tools
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

    def __getitem__(self, i):
        """Return read-counts at position i.

        Arguments:
            i (int): integer of the index of the desired site.

        Returns:
            list: Reads at position i.
        """
        return self.reads[:, i]

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

    def __eq__(self, other):
        """Compares against other gene object.

        Returns:
            bool: True if the gene objects have same orf id.
        """
        return self.orf == other.orf

    def __lt__(self, other):
        """Compares against other gene object.

        Returns:
            bool: True if the gene object id is less than the other.
        """
        return self.orf < other.orf

    def get_gap_span(self):
        """Returns the span of the max_run of the gene (i.e. number of nucleotides).

        Returns:
            int: Number of nucleotides spanned by the max run.
        """
        if len(self.position) > 0:
            if self.r == 0:
                return 0
            index = run_index(self.runs)
            # maxii = numpy.argmax(self.runs)
            maxii = numpy.argwhere(self.runs == numpy.max(self.runs)).flatten()[-1]
            runstart = index[maxii]
            runend = runstart + max(self.runs) - 1
            return self.position[runend] - self.position[runstart] + 2
        else:
            return 0

    def get_gene_span(self):
        """Returns the number of nucleotides spanned by the gene.

        Returns:
            int: Number of nucleotides spanned by the gene's sites.
        """
        if len(self.position) > 0:
            return self.position[-1] - self.position[0] + 2
        return 0

    def theta(self):
        """Return the insertion density ("theta") for the gene.

        Returns:
            float: Density of the gene (i.e. k/n )
        """
        if self.n:
            return float(self.k) / self.n
        else:
            return 0.0

    def phi(self):
        """Return the non-insertion density ("phi") for the gene.

        Returns:
            float: Non-insertion density  (i.e. 1 - theta)
        """
        return 1.0 - self.theta()

    def total_reads(self):
        """Return the total reads for the gene.

        Returns:
            float: Total sum of read-counts.
        """
        return numpy.sum(self.reads, 1)

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

        >>> from pytransit.specific_tools import tnseq_tools
        >>> G = tnseq_tools.Genes(["transit/data/cholesterol_glycerol.transit/glycerol_rep1.wig", "transit/data/cholesterol_glycerol.transit/glycerol_rep2.wig"], "transit/data/genomes/H37Rv.prot_table", norm="TTR")
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

    def __contains__(self, item):
        """Defines __contains__ to check if gene exists in the list.

        Arguments:
            item (str): String with the id of the gene.

        Returns:
            bool: Boolean with True if item is in the list.
        """
        return item in self.orf2index

    def __len__(self):
        """Defines __len__ returning number of genes.

        Returns:
            int: Number of genes in the list.
        """
        return len(self.genes)

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

        is_prot = True
        filename, file_extension = os.path.splitext(self.annotation)
        if file_extension.lower() in [".gff", ".gff3"]:
            is_prot = False

        self.orf2index = {}
        self.genes = []

        orf2info = AnnotationFile(path=self.annotation).orf_to_info
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

        if not no_norm:
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
        with open(self.annotation) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                tmp = line.split("\t")

                if is_prot:
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

    def local_phis(self):
        """Returns numpy array of non-insertion frequency, 'phi', for each gene.

        Returns:
            narray: Numpy array with the complement of density for all genes.
        """
        return 1.0 - self.theta()

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

    def global_run(self):
        """Returns the run assuming all genes were concatenated together.

        Returns:
            int: Max run across all genes.
        """
        return max_run(self.tosses())

    def global_reads(self):
        """Returns the reads among the library.

        Returns:
            list: List of all the data.
        """
        return self.data

    def global_theta(self):
        """Returns global insertion frequency, of the library.

        Returns:
            float: Total sites with insertions divided by total sites.
        """
        return float(self.global_insertion()) / self.global_sites()

    def global_phi(self):
        """Returns global non-insertion frequency, of the library.

        Returns:
            float: Complement of global theta i.e. 1.0-theta
        """
        return 1.0 - self.global_theta()

    def total_reads(self):
        """Returns total reads among the library.

        Returns:
            float: Total sum of read-counts accross all genes.
        """
        reads_total = 0
        for g in self.genes:
            reads_total += g.total_reads()
        return reads_total

    def tosses(self):
        """Returns list of bernoulli trials, 'tosses', representing insertions in the gene.

        Returns:
            list: Sites represented as bernoulli trials with insertions as true.
        """
        all_tosses = []
        for g in self.genes:
            all_tosses.extend(g.tosses)
        return all_tosses

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

def run_index(runs):
    """Returns a list of the indexes of the start of the runs; complements runs().

    Arguments:
        runs (list): List of numeric data.

    Returns:
        list: List of the index of the runs of non-insertions. Non-zero sites are treated as runs of zero.
    """
    index = 0
    index_list = []
    the_run_index = 0
    for r in runs:
        for i in range(r):
            if i == 0:
                the_run_index = index
            index += 1
        if r == 0:
            the_run_index = index
            index += 1
        index_list.append(the_run_index)
    return index_list

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

def get_unknown_file_types(wig_list, transposons):
    """ """
    file_types = set(get_file_types(wig_list))
    method_types = set(transposons)
    extra_types = list(file_types - method_types)
    return extra_types

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
        with open(wig_name) as file:
            for line in file:
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
        with open(path) as file:
            for line in file:
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
        with open(path) as file:
            for line in file:
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

def get_pos_hash_pt(path):
    """Returns a dictionary that maps coordinates to a list of genes that occur at that coordinate.

    Arguments:
        path (str): Path to annotation in .prot_table format.

    Returns:
        dict: Dictionary of position to list of genes that share that position.
    """
    hash = {}
    with open(path) as file:
        for line in file:
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

def get_pos_hash_gff(path):
    """Returns a dictionary that maps coordinates to a list of genes that occur at that coordinate.

    Arguments:
        path (str): Path to annotation in GFF3 format.

    Returns:
        dict: Dictionary of position to list of genes that share that position.
    """
    hash = {}
    with open(path) as file:
        for line in file:
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

def get_coordinate_map(galign_path, reverse=False):
    """Attempts to get mapping of coordinates from galign file.

    Arguments:
        path (str): Path to .galign file.
        reverse (bool): Boolean specifying whether to do A to B or B to A.

    Returns:
        dict: Dictionary of coordinate in one file to another file.
    """
    c1Toc2 = {}
    
    with open(galign_path) as file:
        for line in file:
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

def read_genome(path):
    """Reads in FASTA formatted genome file.

    Arguments:
        path (str): Path to .galign file.

    Returns:
        string: String with the genomic sequence.
    """
    seq = ""
    with open(path) as file:
        for line in file:
            if line.startswith(">"):
                continue
            seq += line.strip()
    return seq

def max_run(lst, item=0):
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

def get_r1(n):
    """Small Correction term. Defaults to 0.000016 for now"""
    return 0.000016

def get_r2(n):
    """Small Correction term. Defaults to 0.00006 for now"""
    return 0.00006

def get_e1(n):
    """Small Correction term. Defaults to 0.01 for now"""
    return 0.01

def get_e2(n):
    """Small Correction term. Defaults to 0.01 for now"""
    return 0.01

def get_gamma():
    """Euler-Mascheroni constant ~ 0.577215664901 """
    return 0.5772156649015328606

def expected_runs(n, pnon):
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
    gamma = get_gamma()
    r1 = get_r1(n)
    E1 = get_e1(n)
    A = math.log(n * pins, 1.0 / pnon)
    B = gamma / math.log(1.0 / pnon)
    ER = A + B - 0.5 + r1 + E1
    return ER

def variance_run(n, pnon):
    """Variance of the expected run of non-insertons (Schilling, 1990):

    .. math::

        variance_run_n =  (pi^2)/(6*ln(1/p)^2) + 1/12 + r2(n) + E2(n)


    Arguments:
        n (int): Integer representing the number of sites.
        pnon (float): Floating point number representing the probability of non-insertion.

    Returns:
        float: Variance of the length of the maximum run.
    """
    r2 = get_r2(n)
    E2 = get_e2(n)
    A = math.pow(math.pi, 2.0) / (6 * math.pow(math.log(1.0 / pnon), 2.0))
    V = A + 1 / 12.0 + r2 + E2
    return V

def gumbel_cdf(x, u, B):
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
            exprun = expected_runs(gene.n, pnon)
            # u is mu of Gumbel (mean=mu+gamma*beta); matching of moments; like Eq 5 in Schilling, but subtract off unneeded terms
            u = exprun - get_gamma() / math.log(1.0 / pnon)
            pval = 1.0 - gumbel_cdf(gene.r, u, B)
            results.append(
                [gene.orf, gene.name, gene.desc, gene.k, gene.n, gene.r, exprun, pval]
            )
    return results

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

def rv_site_indexes_map(genes, ta_site_index_map, n_terminus=0.0, c_terminus=0.0):
    """
    ([Gene], {TAsite: Siteindex}) -> {Rv: Siteindex}
    """
    rv_site_indexes_map = {}
    for g, gene in enumerate(genes):
        site_indexes = []
        start = gene["start"] if gene["strand"] == "+" else gene["start"] + 3
        end = gene["end"] - 3 if gene["strand"] == "+" else gene["end"]
        for i in range(start, end + 1):
            co = i
            if (co - start) / float(end - start) < (n_terminus / 100.0):
                continue
            if (co - start) / float(end - start) > ((100 - c_terminus) / 100.0):
                continue
            if co in ta_site_index_map:
                site_indexes.append(ta_site_index_map[co])
        rv_site_indexes_map[gene["rv"]] = site_indexes
    return rv_site_indexes_map

def extract_yaml_data(comment_lines):
    import ez_yaml
    
    extra_data = {}
    contained_yaml_data = False
    yaml_mode_is_on = False
    yaml_string = "extra_data:\n"
    # 
    # handle header/comments
    # 
    for line in comment_lines:
        if line.startswith("#"):
            if line.startswith("#yaml:"):
                yaml_mode_is_on = True
                contained_yaml_data = True
                continue
            if yaml_mode_is_on and line.startswith("# "):
                yaml_string += f"\n{line[1:]}"
                continue
            else:
                yaml_mode_is_on = False
                # add to the extra_data dict when its done
                if len(yaml_string) > 0:
                    an_object = ez_yaml.to_object(string=yaml_string)
                    extra_data.update(an_object["extra_data"] or {})
    return extra_data

def read_results_file(path):
    from pytransit.generic_tools import csv
    # 
    # get column names
    # 
    comments, headers, rows = csv.read(path, seperator="\t", skip_empty_lines=True, comment_symbol="#")
    comments_string = "\n".join([ "#"+each for each in comments])
    if len(comments) == 0:
        raise Exception(f'''No comments in file, and I expected the last comment to be the column names, while trying to load file "{path}"''')
    column_names = comments[-1].split("\t")
    standardized_column_names = [ col for col in column_names]
    
    extra_data = extract_yaml_data([ "#"+each for each in comments])
    
    # 
    # get rows
    #
    rows_of_dicts = []
    for each_row in rows:
        rows_of_dicts.append({
            each_column_name: each_cell
                for each_column_name, each_cell in zip(standardized_column_names, each_row)
        })
    
    return standardized_column_names, rows_of_dicts, extra_data, comments_string

def wig_id_suggestions_from_filepaths(filepaths):
    import os
    from pytransit.generic_tools import misc
    
    # 
    # Pascal Case of Basenames
    # 
    basenames     = [ os.path.basename(each) for each in filepaths ]
    basenames     = [ each if not each.endswith(".wig") else each[:-4] for each in basenames ]
    first_attempt = [ misc.pascal_case_with_spaces(each) for each in basenames ]
    if misc.all_different(first_attempt):
        return first_attempt
    
    # 
    # Pascal Case of Dirnames
    # 
    dirnames = [ misc.pascal_case_with_spaces(os.path.dirname(each)) for each in filepaths ]
    minimal_dirnames = misc.remove_common_prefix(dirnames)
    fingerprints_base = [
        f"{minimal_dirname}_{minimal_basename}".replace("/","_").replace("\\","_")
            for minimal_dirname, minimal_basename in zip(minimal_dirnames, basenames) 
    ]
    second_attempt = [ misc.pascal_case_with_spaces(each) for each in fingerprints_base ]
    if misc.all_different(second_attempt):
        return second_attempt
    
    # 
    # Minimal simplified dirnames
    # 
    dirnames = [ os.path.dirname(each) for each in filepaths ]
    minimal_dirnames = misc.remove_common_prefix(dirnames)
    fingerprints_base = [
        f"{minimal_dirname}_{basename}".replace("/","_").replace("\\","_")
            for minimal_dirname, basename in zip(minimal_dirnames, basenames) 
    ]
    third_attempt = [ each for each in fingerprints_base ]
    if misc.all_different(third_attempt):
        return third_attempt
    
    # 
    # Minimal dirnames
    # 
    return [
        f"{minimal_dirname}/{minimal_basename}"
            for minimal_dirname, minimal_basename in zip(minimal_dirnames, basenames) 
    ]

class AnnotationFile:
    Info = named_list((
        "name",
        "description",
        "start_coordinate",
        "end_coordinate",
        "strand",
    ))
    
    def __init__(self, *, path):
        """
            Examples:
                print(AnnotationFile(path="./somewhere.gff3").orf_to_info)
                # dictionary with ORF id as keys, and named lists as values
            Summary:
                can load either a .gff or prot_table file
                uses file extension to figure out the difference
        """
        filename, file_extension = os.path.splitext(path)
        file_extension = file_extension.lower()
        self.path = path
        if file_extension in GffFile.file_extensions:
            self.orf_to_info = {
                key : AnnotationFile.Info(value)
                    for key, value in GffFile.extract_gene_info(path).items()
            }
        elif file_extension in ProtTable.file_extensions:
            self.orf_to_info = {
                key : AnnotationFile.Info(value)
                    for key, value in ProtTable.extract_gene_info(path).items()
            }
        else:
            raise Exception(f'''File extension {file_extension} not recognized. Please use .gff or .prot_table file extension''')
    
    def gene_description(self, *, orf_id, fallback_value=None):
        gene_info = self.orf_to_info.get(orf_id, None)
        if gene_info == None:
            return fallback_value
        else:
            return gene_info.description
    
    @property
    def as_list_of_dicts(self):
        gene_list = []
        for each_orf, each_value in self.orf_to_info.items():
            gene_list.append(dict(
                start=each_value.start_coordinate,
                end=each_value.end_coordinate,
                rv=each_orf,
                gene=each_value.name,
                strand=each_value.strand,
            ))
        return gene_list

@misc.singleton
class ProtTable:
    index_of_description = 0
    index_of_gene_start  = 1
    index_of_gene_end    = 2
    index_of_gene_strand = 3
    index_of_gene_name   = 7
    index_of_orf         = 8
    file_extensions = [".prot_table"]
    
    @staticmethod
    def extract_gene_info(path):
        """
            Returns a dictionary that maps gene id to gene information.

            Arguments:
                path (str): Path to annotation in .prot_table format.

            Returns:
                dict: 
                    keys are the ORF id's
                    values are tuples of:
                    - name
                    - description
                    - start coordinate
                    - end coordinate
                    - strand
        """
        orf2info = {}
        with open(path) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                row           = line.strip().split("\t")
                orf           = row[ProtTable.index_of_orf]
                name          = row[ProtTable.index_of_gene_name]
                desc          = row[ProtTable.index_of_description]
                start         = int(row[ProtTable.index_of_gene_start])
                end           = int(row[ProtTable.index_of_gene_end])
                strand        = row[ProtTable.index_of_gene_strand]
                orf2info[orf] = (name, desc, start, end, strand)
        return orf2info

nucleotides_group_size = 3
class GffFile:
    GffRow = named_list((
        # http://gmod.org/wiki/GFF3
        # all of these can be "." or the value mentioned in the comments below
        'seqid',      # name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seq ID must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
        'source',     # name of the program that generated this feature, or the data source (database or project name)
        'type',       # type of feature. Must be a term or accession from the SOFA sequence ontology
        'start',      # Start position of the feature, with sequence numbering starting at 1.
        'end',        # End position of the feature, with sequence numbering starting at 1.
        'score',      # A floating point value.
        'strand',     # defined as + (forward) or - (reverse).
        'phase',      # One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
        'attributes', # A semicolon-separated list of tag-value pairs, providing additional information about each feature. Some of these tags are predefined, e.g. ID, Name, Alias, Parent - see the GFF documentation for more details.
    ))
    file_extensions = [".gff", ".gff3"]
    max_number_of_columns = 9
    
    seqid_index      = 0
    source_index     = 1
    type_index       = 2
    start_index      = 3
    end_index        = 4
    score_index      = 5
    strand_index     = 6
    phase_index      = 7
    attributes_index = 8
    
    def __init__(self, path):
        self.comments, self.columns, self.rows = csv.read(
            path=path,
            seperator="\t",
            first_row_is_column_names=False,
            column_names=GffFile.GffRow.names,
            comment_symbol="#"
        )
        
        # convert "ID=operon001;Name=superOperon" to { "ID": "operon001", "Name": "superOperon" }
        attributes_per_row = (
            dict(
                tuple(each_assignment.split("="))
                    for each_assignment in row.attributes.split(";")
            )
                for row in self.rows
        )
        # overwrite the attribute string with an attribute dictionary
        for each_row, attributes in zip(self.rows, attributes_per_row):
            each_row.attributes = attributes

    @staticmethod
    def is_definitely_not_gff3(path):
        line = ""
        with open(path) as in_file:
            for line in in_file:
                if line.startswith("##gff-version 3"):
                    return False
                if not line.startswith("#"):
                    break
        
        row = line.split("\t")
        if len(row) != GffFile.max_number_of_columns:
            return True
        
        # convert into a named list
        row = GffFile.GffRow(row)
        
        # 
        # check phase is valid
        # 
        if row.phase not in [ '.', '1', '2', '3' ]:
            return True
        # 
        # check strand is valid
        # 
        if row.strand not in [ '.', '+', '-' ]:
            return True
        # 
        # check score is valid
        # 
        if row.score != '.':
            try:
                row.score = float(row.score)
            except ValueError as error:
                return True
        # 
        # check start is valid
        # 
        if row.start != '.':
            try:
                row.start = int(row.start)
            except ValueError as error:
                return True
        # 
        # check end is valid
        # 
        if row.end != '.':
            try:
                row.end = int(row.end)
            except ValueError as error:
                return True
        
        # 
        # consider having two attributes to be a "giveaway" of gff3
        # 
        attributes_list = row.attributes.split(";")
        attribute_assignments = [ each.split("=") for each in attributes_list if len(each) > 0 ]
        # each assignment should be len([varname, value]) == 2
        if len(attribute_assignments) > 0:
            if all(len(assignment) == 2 for assignment in attribute_assignments):
                return False
            else:
                return True
        
        return None
    
    @staticmethod
    def extract_gene_info(path):
        """
            Returns a dictionary that maps gene id to gene information.

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
        gff_data = GffFile(path=path)
        orf2info = {}
        for description, start, end, strand, size, _, _, gene_name, orf_id, _ in gff_data.as_prot_table_rows():
            orf2info[orf] = (gene_name, description, start, end, strand)
        return orf2info
    
    def as_prot_table_rows(self, allowed_types=["CDS"], max_row_length=max_number_of_columns):
        """
            NOTE:
                if you also want tRNAs and rRNAs, modify the allowed_types to include relevent records
        """
        new_rows = []
        for row in self.rows:
            if len(row) > max_row_length:
                continue
            if row.type not in allowed_types:
                continue
            if "locus_tag" not in row.attributes:
                continue
            
            orf_id      = row.attributes["locus_tag"].strip()
            gene_name   = row.attributes.get("gene", "").strip() or "-"
            description = row.attributes.get("product", "")
            size        = int(abs(row.end - row.start + 1) / nucleotides_group_size)
            strand      = row.strand.strip()
            
            new_rows.append(
                [ description, row.start, row.end, strand, size, "-", "-", gene_name, orf_id, "-" ]
            )
        return new_rows


def read_result(path):
    from pytransit.generic_tools.named_list import named_list
    comments, headers, rows = csv.read(path, seperator="\t", skip_empty_lines=True, comment_symbol="#")
    # real headers are in first comment
    headers   = comments[-1].split("\t")
    named_row = named_list([ each for each in headers if each ])
    rows = [ named_row(each) for each in rows ]
    
    return comments, headers, rows
