from random import random

import file_system_py as FS
import ez_yaml

import pytransit.basics.csv as csv
from pytransit.basics.lazy_dict import LazyDict, stringify, indent
from pytransit.basics.named_list import named_list

universal = LazyDict()

class Wig:
    def __init__(self, path, rows=None, extra_data=None):
        self.path            = path
        self.rows            = rows or []
        self.comments        = []
        
        self.extra_data = LazyDict(extra_data)
        if self.extra_data.get("condition", None) is None:
            self.extra_data["condition"] = FS.basename(path)
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

class SessionData:
    def __init__(self):
        self.annotation = None
        self.wigs = []
        self.conditions = []
    
    @property
    def condition_names(self):
        return set([ each.name for each in self.conditions ])
    
    def add_wig(self, path, condition=None):
        self.wigs.append(
            Wig(
                path=path,
                extra_data=LazyDict(
                    is_part_of_cwig=False,
                    condition=condition,
                )
            )
        )
        condition_name = self.wigs[-1].extra_data.condition
        if condition_name not in self.condition_names:
            self.conditions.append(
                Condition(
                    name=condition_name,
                )
            )
    
    def add_cwig(self, cwig_path, metadata_path):
        wig_group = WigGroup.load_from(cwig_path=cwig_path, metadata_path=metadata_path)
        self.wigs += wig_group.wigs
        # add conditions
        for each_wig in wig_group.wigs:
            condition_name = each_wig.extra_data.condition
            if condition_name not in self.condition_names:
                self.conditions.append(
                    Condition(
                        name=condition_name,
                    )
                )
    
    @property
    def files(self):
        return [ self.path for each in self.wigs ]
    
    def import_session(self, path):
        pass # TODO
    
    def export_session(self, path):
        pass # TODO
    
    def __repr__(self):
        return f"""SessionData(
            wigs={indent(self.wigs, by="            ", ignore_first=True)},
            conditions={indent(self.conditions, by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")

        
class CombinedWig:
    def __init__(self, path, comments=None, extra_data=None):
        self.path          = path
        self.rows          = []
        self.comments      = comments or []
        self.extra_data    = extra_data or {}
        self.wigs          = []
        
        self.load()
    
    @property
    def sites(self):
        return [ each.position for each in self.rows ]
    
    @property
    def files(self):
        return self.extra_data["files"]
    
    @files.setter
    def files(self, value):
        self.extra_data["files"] = value
    
    @property
    def read_counts_array(self):
        return [
            each_row[1:len(self.files)] 
                for each_row in self.rows 
        ]
    
    @property
    def read_counts_by_wig(self):
        counts_for_wig = { each_path: [] for each_path in self.files }
        for each_row in self.rows:
            for each_wig_path in self.files:
                counts_for_wig[each_wig_path].append(
                    each_row[each_wig_path]
                )
        return counts_for_wig
    
    def load(self):
        comments, headers, rows = csv.read(self.path, seperator="\t", first_row_is_headers=False)
        comment_string = "\n".join(comments)
        
        sites, counts_by_wig, files, extra_data = [], [], [], {}
        yaml_switch_is_on = False
        yaml_string = "extra_data:\n"
        for line in comments:
            # 
            # handle yaml
            # 
            if line.startswith("#yaml:"):
                yaml_switch_is_on = True
                continue
            if yaml_switch_is_on and line.startswith("# "):
                yaml_string += f"\n"+line[1:]
                continue
            else:
                yaml_switch_is_on = False
                
            # 
            # handle older file method
            # 
            if line.startswith("#File: "):
                files.append(line.rstrip()[7:])  # allows for spaces in filenames
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
            self.wigs.append(
                Wig(
                    path=each_path,
                    rows=list(zip(self.sites, read_counts_by_wig[each_path])),
                    extra_data=LazyDict(
                        is_part_of_cwig=True,
                    ),
                )
            )
        
        return self
    
    def __repr__(self):
        return f"""CombinedWig(
            path={self.path},
            rows_shape=({len(self.rows)}, {len(self.rows[0])}),
            extra_data={indent(self.extra_data, by="            ", ignore_first=True)},
            wigs={      indent(self.wigs      , by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")
        
class CWigMetadata:
    # can contain more data
    def __init__(self, path):
        self.path = path
        self.comments = []
        self.headers = []
        self.rows = []
        self.comments, self.headers, self.rows = csv.read(self.path, seperator="\t", first_row_is_headers=True)
    
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
    
    @property
    def conditions(self):
        return list(set(
            each_row["Condition"]
                for each_row in self.rows
        ))
    
    def __repr__(self):
        return f"""CWigMetadata(
            path={self.path},
            rows_shape=({len(self.rows)}, {len(self.headers)}),
            conditions={indent(self.conditions, by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")
    
    @property
    def conditions(self):
        return list(set(each_row["Condition"] for each_row in self.rows))
    
class WigGroup(LazyDict):
    @staticmethod
    def load_from(cwig_path, metadata_path):
        return WigGroup(
            cwig=CombinedWig(path=cwig_path),
            metadata=CWigMetadata(path=metadata_path),
        )
    
    def __init__(self, cwig, metadata):
        self.cwig = cwig
        self.metadata = metadata
        
        # attach a condition to each wig
        for each_wig in self.cwig.wigs:
            each_wig.extra_data.condition = self.metadata.condition_for(each_wig.path)
            each_wig.extra_data.id = self.metadata.id_for(each_wig.path)
        
    
    @property
    def conditions(self):
        return self.metadata.conditions
    
    @property
    def wigs(self):
        return self.cwig.wigs
    
    def __repr__(self):
        return f"""WigGroup(
            self.cwig={indent(self.self.cwig, by="            ", ignore_first=True)},
            metadata={indent(self.metadata, by="            ", ignore_first=True)},
            conditions={indent(self.conditions, by="            ", ignore_first=True)},
        )""".replace("\n        ", "\n")