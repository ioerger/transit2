from pytransit.basics.named_list import named_list

# reads .csv, .tsv, etc 
def read(path, *, seperator=",", use_headers=None, first_row_is_headers=False, skip_empty_lines=True):
    import json
    
    comments = []
    rows     = []
    headers  = []
    is_first_data_row = True
    
    with open(path,'r') as file:
        for each_line in file.readlines():
            # remove all weird whitespace as a precaution
            each_line = each_line.replace("\r", "").replace("\n", "")
            
            # 
            # comments
            # 
            if each_line.startswith("#"):
                comments.append(each_line[1:])
                continue
            
            # 
            # empty lines
            # 
            if skip_empty_lines and len(each_line.strip()) == 0:
                continue
            
            # 
            # cell data
            #
            cells = each_line.split(seperator)
            cells_with_types = []
            skip_to = 0
            for index, each_cell in enumerate(cells):
                if index < skip_to:
                    continue
                
                stripped = each_cell.strip()
                if len(stripped) == 0:
                    cells_with_types.append(None)
                else:
                    first_char = stripped[0]
                    if not (first_char == '"' or first_char == '[' or first_char == '{'):
                        # this converts scientific notation to floats, ints with whitespace to ints, null to None, etc
                        try: cells_with_types.append(json.loads(each_cell))
                        # if its not valid JSON, just treat it as a string
                        except Exception as error:
                            cells_with_types.append(each_cell)
                    else: # if first_char == '"' or first_char == '[' or first_char == '{'
                        remaining_end_indicies = reversed(list(range(index, len(cells))))
                        skip_to = 0
                        for each_remaining_end_index in remaining_end_indicies:
                            try:
                                cells_with_types.append(
                                    json.loads(seperator.join(cells[index:each_remaining_end_index]))
                                )
                                skip_to = each_remaining_index
                                break
                            except Exception as error:
                                pass
                        # continue the outer loop
                        if skip_to != 0:
                            continue
                        else:
                            # if all fail, go with the default of the shortest cell as a string
                            cells_with_types.append(each_cell)
            
            # 
            # headers
            # 
            if is_first_data_row:
                is_first_data_row = False
                if first_row_is_headers:
                    headers.append([ str(each) for each in cells_with_types ])
                    continue
            
            rows.append(cells_with_types)
    
    # if headers
    if first_row_is_headers or use_headers:
        RowItem = named_list(headers[0] or use_headers)
        # tranform each into a named list (backwards compatible with regular list)
        rows = [ RowItem(each_row) for each_row in rows ]
    
    return comments, headers, rows

def write(path, comments, headers, rows, seperator=","):
    import json
    
    def element_to_string(element):
        # strings are checked for seperators, if no seperators or whitespace, then unquoted
        if isinstance(element, str):
            if not (
                seperator in element or
                '\n' in element or
                '\r' in element or
                '\t' in element or
                ' ' in element
            ):
                # no need for quoting
                return element
        # all other values are stored in json format
        else:
            return json.dumps(element)
            
    
    with open(path, 'w') as the_file:
        # 
        # comments
        # 
        the_file.write(
            "\n".join([ f"#{each}" for each in comments ])
        )
        
        # 
        # headers
        # 
        if len(headers) > 0:
            the_file.write(
                seperator.join(tuple(
                    element_to_string(str(each)) for each in headers 
                ))
            )
        
        # 
        # rows
        # 
        for each_row in rows:
            the_file.write(
                seperator.join(tuple(
                    element_to_string(each_cell)
                    for each_cell in each_row 
                ))
            )