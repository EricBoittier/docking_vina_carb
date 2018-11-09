file = open("/Users/ericboittier/Desktop/GAGTOOLS/1AXM.hb2", "r")
lines = file.readlines()

for x in lines[0:8]:
    print(x)

hbonds = {}

for line in lines:
    if line.__contains__("IDS") or line.__contains__("SGN"):
        split = line.split()
        print(line)

        if split[2] not in hbonds.keys():
            hbonds[split[2]] = []

        hbonds[split[2]].append([split[0], split[4]])

        print("Donor: {} Acceptor: {} Distance: {}".format(split[0], split[2], split[4]))


freq = {}

for key in hbonds.keys():
    res = key.split("-")[1]

    if res not in freq.keys():
        freq[res] = {}

    for hbond in hbonds[key]:
        if hbond[0].split("-")[1] not in freq[res].keys():
            freq[res][hbond[0].split("-")[1]] = 0
        freq[res][hbond[0].split("-")[1]] += 1
print("")
print(freq)