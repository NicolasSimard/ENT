def print_row(row,name=""):
    tmp = "{:9}|".format(name if len(name) != 0 else "***")
    data = list(map(lambda x:list(x.split(",")),list(row.split("["))))
    data.pop(0)
    data = [[int(tup[n].strip(" []")) for n in range(2)] for tup in data]
    #print(len(data),*data,sep="::",end="\n\n")
    # print(*list(data[0].split(",")),sep="::")
    return tmp + "|".join(["({0[0]}:{0[1]:5})".format(data[n]) for n in range(len(data))]) + "\n"
    
def print_header(N):
    return "***      |" + "|".join(["{:9}".format(n) for n in range(N)]) + "\n"\
    + "|".join(["---------".format(n) for n in range(N+1)]) + "\n"
    
if __name__ == "__main__":
    import sys
    p = int(sys.argv[1])
    file = open(sys.argv[2],'r')
    output = open(sys.argv[2] + ".md",'w')
    for line in file:
        line=line[1:-2] # Drop the opening [ and closing ] of the matrix
        rows = line.split(";") # Split in rows
        output.write(print_header(len(rows[0].split("["))-1))
        # print(*rows,sep=";\n")
        n=1
        for row in rows:
            output.write(print_row(row,str(p)+"^"+str(n)))
            n+=1
    file.close()
    output.close()