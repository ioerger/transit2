from pytransit.basics.named_list import named_list

# reads .csv, .tsv, etc 
def read(path, *, seperator=",", first_row_is_column_names=False, column_names=None, skip_empty_lines=True, comment_symbol=None):
    """
        Examples:
            comments, column_names, rows = csv.read("something/file.csv", first_row_is_column_names=True, comment_symbol="#")
            comments, _empty_list, rows = csv.read("something/file.csv", first_row_is_column_names=False)
            comments, column_names_from_file, rows = csv.read(
                "something/file.csv",
                column_names=["column1_new_name"],
                first_row_is_column_names=True,
            )
        Summary:
            Reads in CSV's
            - Converts numbers, null, booleans, etc into those types in accordance with JSON
              (e.g. null=>None, true=>True, 2.3e31=>float, "hi\n"=>str('hi\n'))
            - Anything that is not json-parsable is kept as a string
            - Comments can be enabled by with the comment_symbol arg
            - Comments must start as the first character of a line, no trailing comments
            - Blank spaces (e.g. ,,, ) are converted to None (e.g. ,null,null,)
            - Read() will sill parse even if some lines are missing columns
        Returns:
            value: tuple(comments, column_names, rows)
            rows:
                - Always returns a list
                - Each element is a named list
                - Named lists inherit from lists (full backwards compatibility)
                - Named lists may also be accessed using column_names
                  for example: rows[0]["column1"] and rows[0].column1 are both valid
            column_names:
                - Will always return an empty list when first_row_is_column_names=False
                - Will always be the column names according to the file (even if overridden)
                - Every element in the list will be a string
            comments:
                - A list of strings
                - One string per line
                - The comment_symbol itself is not removed from the string
            
        Arguments:
            path:
                - Any string path or path-object
                - Will throw error if file does not exist
            first_row_is_column_names:
                - Boolean, default is False
                - If true all elements in the first row will be parsed as strings (even if they look like numbers/null)
                - Not all columns need a name
                - However using the same name twice or more will cause problems
            column_names:
                - Optional, a list of strings
                - Will override the column_names within the file if provided
                - Doesn't need to cover all columns (trailing columns can be unnamed)
            skip_empty_lines:
                - Boolean, default is True
                - A line with spaces or tabs will still qualify as empty
            
    """
    import json
    
    comments     = []
    rows         = []
    file_column_names = []
    is_first_data_row = True
    
    with open(path,'r') as file:
        for each_line in file.readlines():
            # remove all weird whitespace as a precaution
            each_line = each_line.replace("\r", "").replace("\n", "")
            
            # 
            # comments
            # 
            if comment_symbol:
                if each_line.startswith(comment_symbol):
                    comments.append(each_line[len(comment_symbol):])
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
            # file_column_names
            # 
            if is_first_data_row:
                is_first_data_row = False
                if first_row_is_column_names:
                    file_column_names = [ str(each) for each in cells_with_types ]
                    continue
            
            rows.append(cells_with_types)
    
    # if file_column_names
    if first_row_is_column_names or column_names:
        RowItem = named_list(column_names or file_column_names)
        # tranform each into a named list (backwards compatible with regular list)
        rows = [ RowItem(each_row) for each_row in rows ]
    
    return comments, file_column_names, rows

def write(path, *, rows, column_names=[], seperator=",", eol="\n", comment_symbol=None, comments=[]):
    import json
    assert comment_symbol or len(comments) == 0, "Comments were provided,"
    def contains_comment_symbol(string):
        if not comment_symbol:
            return False
        else:
            return comment_symbol in string
    
    def element_to_string(element):
        # strings are checked for seperators, if no seperators or whitespace, then unquoted
        if isinstance(element, str):
            if not (
                contains_comment_symbol(element) or
                seperator in element or
                eol in element or
                '\n' in element or
                '\r' in element
            ):
                # no need for quoting
                return element
        # all other values are stored in json format
        return json.dumps(element)
    
    def break_up_comments(comments):
        for each in comments:
            yield from f"{each}".replace("\r", "").split("\n")
    
    with open(path, 'w') as the_file:
        # 
        # comments
        # 
        the_file.write(
            eol.join([ f"{comment_symbol}{each}" for each in break_up_comments(comments) ])
        )
        if len(comments) > 0:
            the_file.write(eol)
        
        # 
        # column_names
        # 
        if len(column_names) > 0:
            the_file.write(
                seperator.join(tuple(
                    element_to_string(str(each)) for each in column_names 
                ))+eol
            )
        
        # 
        # rows
        # 
        for each_row in rows:
            row_string_escaped = tuple(
                element_to_string(each_cell)
                    for each_cell in each_row 
            )
            line = seperator.join(row_string_escaped)+eol
            the_file.write(
                seperator.join(row_string_escaped)+eol
            )