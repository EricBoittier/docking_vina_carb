# import os
# import re
#
# file =  open("sites.txt", "w")
#
# regex = re.compile("SITE\s\s\s\s\s\d")
#
# for protein in os.listdir("."):
#     if protein.__contains__(".pdb"):
#         file.write("\n")
#         t = False
#         file.write(protein+"\n")
#         print(protein)
#         with open(protein, "r") as f:
#             for line in f.readlines():
#                 for match in regex.findall(line):
#                     t = True
#                 if t:
#                     file.write(line)
#                 t = False
#         if t:
#             print("yes")
#         else:
#             print("\n\n\nno\n\n\n")

with open("sites.txt", "r") as file:
    for line in file.readlines():
        s = line.split()
        count = 0
        for split in s:
            print('{} {}'.format(count, split))
            count+=1
        print("")


# 4 7 10 13