#!/usr/bin/env python3
import tinydb as db
from pprint import pprint

def main():

    latDB = db.TinyDB("../latDB.json")
    pars = db.Query()

    db_keys = []

    for item in latDB:
        rec = dict(item)
        key = rec["key"]
        vals = rec["vals"]
        # tmp = key.split("_")
        # print(key)
        db_keys.append(key)

    for key in sorted(db_keys):
        print(key)


    # example of iterating over the DB
    # calDB = db.TinyDB('../calDB.json')
    # for item in calDB:
    #     d = dict(item)
    #     key = d["key"]
    #     vals = d["vals"]
    #     tmp = key.split("_")
    #     tmp = [str(t) for t in tmp]
    #
    #     if tmp[0]=="fitSlo" and tmp[1]=="ds%d" % dsNum:
    #         print(tmp)
    #         # print(vals)
    #         # return
    #
    #     # nRec += 1
    #     # if nRec > 10: break
    #
    # use a regexp to print a buncha records
    # calDB = db.TinyDB('../calDB.json')
    # pars = db.Query()
    # recList = calDB.search(pars.key.matches("riseNoise"))
    # for idx in range(len(recList)):
    #     key = recList[idx]['key']
    #     vals = recList[idx]['vals']
    #     print(key)



if __name__=="__main__":
    main()