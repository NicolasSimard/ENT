Dcolprinted = False

""" The use is as follows:

    comps-to-textable.py file maxell nbrD out
    
    where:
    - file: the file containing the data in the following format:
    D1
    f_D1(1)
    f_D1(2)
    ...
    f_D1(maxell)
    D2
    f_D2(1)
    ...
    f_DnbrD(maxell)
    - maxell: the maximal value of ell
    - nbrD : the number of discriminants
    - out: the file to which we print the markdown table
    
    The output is a Latex table with a specific format.
"""
def row(rowdata):
    global Dcolprinted
    if not Dcolprinted:
        col1 = "\\multicolumn{{1}}{{|c|}}{{\\multirow{{{0}}}{{*}}{{$D$}}}}".format(sys.argv[3])
        Dcolprinted = True
    else:
        col1 = "\multicolumn{1}{ |c| }{}"
    col2 = "& \\multicolumn{{1}}{{|c|}}{{{0}}}\n".format(rowdata.pop(0))
    rest = "& " + " & ".join(rowdata) + "\\\\\n"
    line = "\\cline{{2-{}}}\n".format(maxell+2)
    return col1+col2+rest+line
    
def header(maxell):
    line1 = "\\begin{{tabular}}{{cc|*{{{0}}}{{l|}}}}\n".format(maxell)
    line2 = "\\cline{{3-{}}}\n".format(maxell+2)
    line3 = "& & \\multicolumn{{{0}}}{{ c| }}{{$\\ell$}} \\\\ \\cline{{3-{1}}}\n".format(maxell,maxell+2)
    line4 = "& & "+" & ".join([str(i) for i in range(1,maxell+1)])+"\\\\ \\cline{{1-{0}}}\n".format(maxell+2)
    return line1+line2+line3+line4

def footer():
    return "\\end{tabular}"
    
if __name__ == "__main__":
    import sys
    file = open(sys.argv[1],'r')
    maxell = int(sys.argv[2])
    out = open(sys.argv[4],'w')
    out.write(header(maxell))
    rowdata = []
    for line in file:
        rowdata.append(line.strip())
        if len(rowdata) == maxell + 1:
            out.write(row(rowdata))
            rowdata = []
    out.write(footer())
    file.close()
    out.close()
    