import dsi
det = dsi.DetInfo()

def main():

    parseLTDoc()

def parseLTDoc():

    with open("/Users/wisecg/dev/sandbox/vince/ds05-exp.txt") as f:
        lines = f.readlines()

    detList = det.allDets
    totExpo = {cpd:0 for cpd in detList}

    grandTot = 0

    for line in lines[2:]:

        tmp = line.rstrip().split()
        if len(tmp)==0: continue
        if tmp[0]=="DS":
            ds = int(tmp[1])
            print("new DS:",ds)
            continue
        if not tmp[0].isdigit(): continue

        if ds!=0:break

        chan = int(tmp[0])
        if chan%2==1:
            print("lg chan,",chan)
            chan -= 1 # hack to combine LG exposure with the HG exposure

        cpd = det.getChanCPD(ds, chan)
        expo = float(tmp[5])
        grandTot += expo

        print(chan, cpd, expo)

        totExpo[cpd] += expo

    # print("summary")
    # for cpd in sorted(totExpo):
    #     print(cpd, totExpo[cpd])

    print("grand total",grandTot)

if __name__=="__main__":
    main()