from sys import argv


with open(argv[1]) as infile:
    infile.readline()
    for line in infile:
        line = line.strip().split('.t1')
        gid = line[0]+'.pdb'
        line[1].lstrip(" ").split(" ")
        others = "\t".join(line[1])
        print(f"{gid}\t{others}")
