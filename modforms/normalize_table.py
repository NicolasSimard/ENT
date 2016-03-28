def convert_line(entries):
    p,e = list(map(int,entries.pop(0).split("^")))
    phi = p**e-p**(e-1)
    temp = "{:2}^{}    ".format(p,e)
    for k_n in entries:
        try:
            temp += "|({}:{:5})".format((int(k_n)-w)//phi,int(k_n))
        except ValueError:
            temp += "|{:9}".format(k_n)
    temp += "\n"
    return temp

if __name__ == "__main__":
    w=2
    import sys
    file = open(sys.argv[-1],'r')
    output = open(sys.argv[-1] + ".norm",'w')
    for line in file:
        if line.find("|") != -1:
            entries = list(map(lambda x: x.strip(" *"),line.split("|")))
            try:
                int(entries[0].split("^")[0])
                output.write(convert_line(entries))
            except ValueError:
                output.write(line)
        else:
            output.write(line)
    file.close()
    output.close()
    