def write_dat(path, heading, table, eol="\n"):
    if len(heading) != 0:
        heading = "#" + heading
    heading = heading.replace("\n", "\n#")
    body = eol.join([ "\t".join(each_row) for each_row in table ])
    string = heading + eol + body
    with open(path, 'w') as outfile:
        outfile.write(string, outfile)