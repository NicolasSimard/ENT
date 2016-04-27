""" The use is as follows:

    mat_to_md.py file p min max out
    
    where:
    - file: the file containing the matrix
    - p: the prime whose powers label the rows
    - min: the begining of the numbering of the top rows
    - max: the end of the numbering of the top rows (inclusive)
    - out: the file to which we print the markdown table
    
    ex:
    Suppose the file ex.data contains the matrix
    [a,b;c,d]
    
    Then the command
    mat_to_md.py ex.data 3 3 5 table.md
    
    will produce a file called table.md containing the following text:
    ***      |   3     |     4
    ---------|---------|---------
    3^1      |a        |b
    3^2      |b        |d
    
"""
def print_row(row,name=""):
    tmp = "{:9}|".format(name if len(name) != 0 else "***")
    data = list(map(lambda x:list(x.split(",")),list(row.split("["))))
    data.pop(0)
    data = [[int(tup[n].strip(" []")) for n in range(2)] for tup in data]
    return tmp + "|".join(["({0[0]}:{0[1]:5})".format(data[n]) for n in range(len(data))]) + "\n"
    
def print_header(min,max):
    return "***      |" + "|".join(["{:9}".format(n) for n in range(min,max+1)]) + "\n"\
    + "|".join(["---------".format(n) for n in range(max-min+2)]) + "\n"
    
if __name__ == "__main__":
    import sys
    file = open(sys.argv[1],'r')
    p = int(sys.argv[2])
    min = int(sys.argv[3])
    max = int(sys.argv[4])
    out = open(sys.argv[5],'w')
    for line in file:
        line=line[1:-2] # Drop the opening [ and closing ] of the matrix
        rows = line.split(";") # Split in rows
        out.write(print_header(min,max))
        n=1
        for row in rows:
            out.write(print_row(row,str(p)+"^"+str(n)))
            n+=1
    file.close()
    out.close()
    