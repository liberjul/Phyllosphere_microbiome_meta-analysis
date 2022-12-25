import sys

set1 = set()
set2 = set()
with open(sys.argv[1], "r") as ifile1, open(sys.argv[2], "r") as ifile2:
    line = ifile1.readline()
    while line != "":
        set1.add(line.strip("\n").split("\t")[-1])
        line = ifile1.readline()
    line = ifile2.readline()
    while line != "":
        set2.add(line.strip("\n").split("\t")[-1])
        line = ifile2.readline()

set_diff = list(set1 - set2)
# print(F'{len(set1)}, {len(set2)}, {len(set_diff)}')
for i in set_diff:
    print(i)
