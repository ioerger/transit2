def read_cwig(path):
    sites, counts_by_wig, files, metadata = [], [], [], {}
    with open(path) as f:
        yaml_mode_is_on = False
        yaml_string = "metadata:\n"
        for line in f.readlines():
            # 
            # handle yaml
            # 
            if line.startswith("#yaml:"):
                yaml_mode_is_on = True
                continue
            if yaml_mode_is_on and line.startswith("# "):
                yaml_string += f"\n{line[-1:]}"
                continue
            else:
                yaml_mode_is_on = False
                # add to the metadata dict when its done
                if len(yaml_string) > 0:
                    metadata.update(ez_yaml.to_object(string=yaml_string)["metadata"])
                    files += metadata.get('files',[])
                    files = list(set(files)) # remove any duplicate entries
            
            # 
            # handle older file method
            # 
            if line.startswith("#File: "):
                files.append(line.rstrip()[7:])  # allows for spaces in filenames
                continue
            
            # 
            # handle comments and empty lines
            # 
            if len(line) == 0 or line[0] == "#":
                continue
            
            #
            # actual parsing
            #
            cols = line.split("\t")[0 : 1+len(files)]
            cols = cols[: 1+len(files)]  # additional columns at end could contain gene info
            # Read in position as int, and readcounts as float
            cols = list(
                map(
                    lambda t_iv: int(t_iv[1]) if t_iv[0] == 0 else float(t_iv[1]),
                    enumerate(cols),
                )
            )
            position, wig_counts = cols[0], cols[1:]
            sites.append(position)
            for index, count in enumerate(wig_counts):
                if len(counts_by_wig) < index+1:
                    counts_by_wig.append([])
                counts_by_wig[index].append(count)
    
    combined_wig = CombinedWig(numpy.array(sites), numpy.array(counts_by_wig), files)
    combined_wig.metadata = metadata
    return combined_wig

class Wig:
    def __init__(self, path, part_of_cwig, condition=None):
        self.condition = None

class ExperimentData:
    def __init__(self):
        self.wig_paths = []
        self.wig_conditions = {}
    
    def add_wig(self, wig_filepath, condition=None):
        pass # TODO: auto-assign condition if none given
    
    def add_cwig(self, cwig_filepath, metadata_filepath):
        read_cwig(cwig_filepath)
        pass # TODO:
    
    def import_session(self, path):
        pass # TODO
    
    def export(self, path):
        pass # TODO
    
    @property
    def conditions(self):
        pass # TODO
    
    @property
    def files(self):
        return self[2]
    
    
    
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
            return (numpy.zeros((1, 0)), numpy.zeros(0), [])

        # Check size of all wig file matches
        size_list = []
        for _, path in enumerate(list_of_paths):
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

        data = numpy.zeros((len(list_of_paths), line_count))
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
                    data[path_index, line_index] = rd
                except Exception as error:
                    raise Exception(f'''
                        
                        Make sure that all wig files have the same number of TA sites (i.e. same strain)
                        
                        Original Error:\n{error}
                        
                    ''')
                position[line_index] = pos
                line_index += 1
        
        return data, position
    
    @classmethod
    def load(cls, cwig_path, condition_data_path):
        combined_wig_object  = cls.read_cwig(cwig_path)
        conditions_by_file, covariates_by_file_list, interactions_by_file_list, ordering_metadata = cls.read_condition_data(condition_data_path)
        # attach the data to the object
        combined_wig_object.conditions = ConditionData(
            conditions_by_file=conditions_by_file,
            covariates_by_file_list=covariates_by_file_list,
            interactions_by_file_list=interactions_by_file_list,
            ordering_metadata=ordering_metadata,
        )
        return combined_wig_object
        
    @classmethod
    def read_condition_data(cls, path, covars_to_read=[], interactions_to_read=[], condition_name="Condition"):
        """
        Filename -> ConditionMap
        ConditionMap :: {Filename: Condition}, [{Filename: Covar}], [{Filename: Interaction}]
        Condition :: String
        Covar :: String
        Interaction :: String
        """
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

        return (
            conditions_by_file,
            covariates_by_file_list,
            interactions_by_file_list,
            ordering_metadata,
        )
    
    @classmethod
    def read_cwig(cls, fname):
        """
            Read the combined wig-file generated by Transit
            :: Filename -> Tuple([Site], [WigData], [Filename])
            Site :: Integer
            WigData :: [Number]
            Filename :: String
        """
        # TODO: 
        
        
        
cwig = CombinedWig()        
cwig.metadata["RefGenome"] # 