try:
    import wx
    WX_VERSION = int(wx.version()[0])
    hasWx = True
except Exception as error:
    hasWx = False
    WX_VERSION = 0

def write_dat(path, heading, table, eol="\n"):
    if len(heading) != 0:
        heading = "#" + heading
    heading = heading.replace("\n", "\n#")
    body = eol.join([ "\t".join(each_row) for each_row in table ])
    string = heading + eol + body
    with open(path, 'w') as outfile:
        outfile.write(string, outfile)