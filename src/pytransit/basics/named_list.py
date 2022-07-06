def named_list(names):
    """
    Example:
        Position = named_list(['x','y','z'])
        a = Position([1,2,3])
        print(a.x)   # 1
        a.x = 4
        print(a[0])  # 4
        a[0] = 9
        print(a.x)   # 9
    """
    
    class NamedList(list):
        def __getitem__(self, key):
            if isinstance(key, int):
                return super(NamedList, self).__getitem__(key)
            # assume its a name
            else:
                index = names.index(key)
                if index >= len(self):
                    return None
                return self[index]
        
        def __getattr__(self, key):
            if key in names:
                return self[key]
            else:
                super(NamedList, self).__getattr__(key)
        
        def __setattr__(self, key, value):
            if key in names:
                index = names.index(key)
                while index >= len(self):
                    super(NamedList, self).append(None)
                super(NamedList, self).__setitem__(index, value)
            else:
                super(NamedList, self).__setattr__(key, value)
        
        def __setitem__(self, key, value):
            if isinstance(key, int):
                super(NamedList, self).__setitem__(key, value)
            # assume its a name
            else:
                index = names.index(key)
                while index >= len(self):
                    super(NamedList, self).append(None)
                super(NamedList, self).__setitem__(index, value)
        
        def __repr__(self):
            import itertools
            out_string = '['
            named_values = 0
            for each_value, each_name in zip(self, names):
                named_values += 1
                out_string += f' {each_name}={each_value},'
            # has unnamed values
            length = len(self)
            if named_values < length:
                for index in itertools.count(named_values-1): # starting where we left off
                    if index >= length:
                        break
                    value = self[index]
                    out_string += f' {value},'
            
            out_string += ' ]'
            return out_string
            
    return NamedList

