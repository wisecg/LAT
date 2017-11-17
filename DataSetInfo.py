# DataSetInfo.py
# C. Wiseman, B. Zhu
# v1. 20 June 2017
# v2. 12 Sept 2017
import numpy as np

# ==================================================================================
#                                RUN INFORMATION
# Taken from DataSetInfo.hh from the GAT version on Sep 12th 2017.
# ==================================================================================

# number of BG run ranges
dsMap = {0:75,1:51,2:7,3:24,4:18,5:112,6:26}

# runs must cover bg list and calibration runs
dsRanges = {
    0:[2571,7614],
    1:[9407,14502],
    2:[14699,15892],
    3:[16797,18589],
    4:[60000791,60002394],
    5:[18623,23958]
    6:[25672,100000]
    }

# calibration master list, and some lookup functions
class CalInfo:
    def __init__(self):
        self.master = {}
        self.master["ds0_m1"] = {
            0 : [[2571,2579],2571,2629],
            1 : [[2630,2643],2630,2649],
            2 : [[2650,2657],2650,2673],
            3 : [[2674,2681],2674,2930],
            4 : [[2931,2974],2931,3056],
            5 : [[3057,3077],3057,3129],
            6 : [[3130,3136],3130,3272],
            7 : [[3273,3275,3277,3278,3280,3280],3273,3280],
            8 : [[3281,3292],3281,3432],
            9 : [[3433,3460],3433,3645],
            10: [[3646,3662],3646,3687],
            11: [[3688,3708],3688,3982],
            12: [[3983,4003],3983,4134],
            13: [[4135,4170],4135,4520],
            14: [[4521,4546],4520,4831],
            15: [[4833,4840,4842,4853],4832,4907],
            16: [[4908,4937],4908,4981],
            17: [[4982,5006],4982,5061],
            18: [[5062,5089],5062,5252],
            19: [[5253,5276],5253,5330],
            20: [[5331,5368],5331,5414],
            21: [[5415,5448],5415,5501],
            22: [[5502,5524],5502,5850],
            23: [[5851,5871],5851,5910],
            24: [[5911,5921],5911,5993],
            25: [[5994,6014],5994,6168],
            26: [[6169,6189],6169,6317],
            27: [[6318,6330],6318,6352],
            28: [[6353,6365],6353,6544],
            29: [[6545,6552],6545,6577],
            30: [[6578,6598],6578,6753],
            31: [[6754,6774],6754,6853],
            32: [[6854,6886],6854,6903],
            33: [[6904,6925],6904,7274],
            34: [[7275,7279,7281,7295],7275,7614]
            }
        self.master["ds1_m1"] = {
            0: [[9407,9412,9414,9420],9407,9452],
            1: [[9453,9454,9456,9467],9453,9496],
            2: [[9497,9504],9497,9520],
            3: [[9521,9521,9523,9535],9521,9698],
            4: [[9699,9707],9699,9737],
            5: [[9738,9758],9738,9787],
            6: [[9788,9811],9788,9937],
            7: [[9938,9942,9944,9950],9938,9997],
            8: [[9998,10014],9998,10038],
            9: [[10039,10059],10039,10092],
            10: [[10093,10112],10093,10239],
            11: [[10240,10260],10240,10357],
            12: [[10358,10368],10358,10401],
            13: [[10402,10420],10402,10496],
            14: [[10497,10505],10497,10528],
            15: [[10529,10545],10529,10720],
            16: [[10721,10721,10723,10727,10729,10735],10721,10766],
            17: [[10767,10783],10767,10809],
            18: [[10810,10824],10810,10846],
            19: [[10847,10861],10847,10942],
            20: [[10943,10956],10943,10982],
            21: [[10983,10996],10983,11028],
            22: [[11029,11042],11029,11067],
            23: [[11068,11081],11068,11214],
            24: [[11215,11223],11215,11306],
            25: [[11307,11321],11307,11339],
            26: [[11340,11347],11340,11488],
            27: [[11489,11496],11489,11506],
            28: [[11507,11527],11507,11798],
            29: [[11799,11808],11799,12501],
            30: [[12502,12509],12502,12626],
            31: [[12627,12634],12627,12655],
            32: [[12656,12662],12656,12725],
            33: [[12726,12733],12726,12799],
            34: [[12800,12808],12800,12863],
            35: [[12864,12873],12864,12885],
            36: [[12886,12910],12886,12968],
            37: [[12969,12998],12969,13057],
            38: [[13058,13064],13058,13138],
            39: [[13139,13146],13139,13352],
            40: [[13353,13360],13353,13384],
            41: [[13385,13392],13385,13419],
            42: [[13420,13427],13420,13549],
            43: [[13550,13556],13550,13689],
            44: [[13690,13697],13690,13705],
            45: [[13706,13713],13706,13739],
            46: [[13740,13747],13740,13772],
            47: [[13773,13780],13773,13904],
            48: [[13905,13912],13905,13953],
            49: [[13954,13962,13964,13965],13954,13986],
            50: [[13987,13994],13987,14092],
            51: [[14093,14100],14093,14148],
            52: [[14149,14156],14149,14179],
            53: [[14180,14188],14180,14273],
            54: [[14274,14282],14274,14316],
            55: [[14317,14325],14317,14375],
            56: [[14376,14384],14376,14502]
            }
        self.master["ds2_m1"] = {
            0: [[14699,14706],14699,14854],
            1: [[14855,14861],14855,14926],
            2: [[14927,14933],14927,15053],
            3: [[15054,15060],15054,15247],
            4: [[15248,15254],15248,15329],
            5: [[15330,15336],15330,15489],
            6: [[15490,15506],15490,15626],
            7: [[15627,15632],15627,15706],
            8: [[15707,15713],15707,15788],
            9: [[15789,15795],15789,15844],
            10: [[15845,15858],15845,15892]
            }
        self.master["ds3_m1"] = {
            0: [[16836,16854],16797,16910],
            1: [[16911,16929],16911,17016],
            2: [[17017,17033],17017,17182],
            3: [[17183,17203],17183,17422],
            4: [[17423,17441],17423,17520],
            5: [[17521,17529],17521,17687],
            6: [[17688,17696],17688,17956],
            7: [[17957,17965],17957,18327],
            8: [[18328,18345],18328,18589]
            }
        self.master["ds4_m2"] = {
            0: [[60000791,60000800],60000791,60001013],
            1: [[60001014,60001028],60001014,60001147],
            2: [[60001148,60001161],60001148,60001206],
            3: [[60001207,60001237],60001207,60001351],
            4: [[60001352,60001365],60001352,60001441],
            5: [[60001446,60001461],60001446,60001507],
            6: [[60001508,60001521],60001508,60001577],
            7: [[60001578,60001592],60001578,60001719],
            8: [[60001720,60001732],60001720,60001854],
            9: [[60001855,60001868],60001855,60001913],
            10: [[60001914,60001926],60001914,60002394]
            }
        self.master["ds5_m1"] = {
            0: [[18713,18732],18623,19093],
            1: [[19094,19109],19094,19242],
            2: [[19243,19258],19243,19519],
            3: [[19520,19542],19520,19750],
            4: [[19751,19761],19751,20041],
            5: [[20042,20056],20042,20200],
            6: [[20201,20216],20201,20446],
            7: [[20447,20464],20447,21179],
            8: [[21180,21193],21180,21436],
            9: [[21437,21451],21437,21801],
            10: [[21802,21817],21802,21969],
            11: [[21970,21998],21970,22281],
            12: [[22282,22292],22282,22512],
            13: [[22513,22566],22513,22853],
            14: [[22854,22865],22854,23481],
            15: [[23482,23496],23482,23886],
            16: [[23887,23887,23889,23889,23891,23907,23909,23928],23883,23958]
            }
        self.master["ds5_m2"] = {
            0: [[18746,18750,18753,18756,18759,18759],18623,19054],
            1: [[19055,19069],19055,19585],
            2: [[19586,19597],19586,19761],
            3: [[19762,19769],19762,20064],
            4: [[20065,20072],20065,20465],
            5: [[20466,20481],20466,20555],
            6: [[20556,20560,20562,20563,20565,20569,20571,20572,20574,20578,20580,20581,20583,20583],20556,21202],
            7: [[21203,21214],21203,21453],
            8: [[21454,21467],21454,21818],
            9: [[21819,21831],21819,22293],
            10: [[22294,22302],22294,22567],
            11: [[22568,22635],22568,22841],
            12: [[22842,22851],22842,23498],
            13: [[23498,23509],23482,23958]
            }

        # Track all the 'hi' run coverage numbers for fast run range lookups
        self.covIdx = {}
        for key in self.master:
            tmp = []
            for idx in self.master[key]:
                tmp.append(self.master[key][idx][2])
            self.covIdx[key] = np.asarray(tmp)

    def GetMasterList(self): return self.master
    def GetCovArr(self,key): return self.covIdx[key]
    def GetIdxs(self,key): return len(self.covIdx[key])

    def GetKeys(self,dsNum=None):
        keyList = sorted(self.master.keys())
        if dsNum==None:
            return keyList
        else:
            thisDSList = []
            for key in keyList:
                if "ds%d" % dsNum in key: thisDSList.append(key)
            return thisDSList

    def GetCalIdx(self,key,run):
        """ Look up the calibration index corresponding to a particular run. """
        if key not in self.covIdx:
            print "Key %s not found in master list!" % key
            return None
        else:
            idx = np.searchsorted(self.covIdx[key], run)
            if idx not in self.master[key]:
                print "Run %d out of range of key %s.  calIdx was %d" % (run, key, idx)
                return None
            lst = self.master[key][idx]
            lo, hi = lst[1], lst[2]
            if lo <= run <= hi:
                return idx
            else:
                print "Run %d not found with key %s, lo=%d hi=%d" % (run,key,lo,hi)
                return None

    def GetCalList(self,key,idx,runLimit=None):
        """ Generate a list of runs for a given calibration index. """
        if key not in self.master:
            print "Key %s not found in master list!" % key
            return None
        else:
            runList = []
            if idx not in self.master[key]:
                return None
            lst = self.master[key][idx][0]
            for i in xrange(0,len(lst),2):
                lo, hi = lst[i], lst[i+1]
                runList += range(lo, hi+1)
            if runLimit is not None:
                del runList[10:]
            return runList


# ==================================================================================
#                          RAW CHANNEL & DETECTOR INFORMATION
# These are 'raw' lists of ALL ENABLED CHANNELS throughout the DS's.
# Generated by skim-checks.cc on 20 June 2017.
# (Could move that routine to $GATDIR/Apps/check-data.cc)
# ==================================================================================

DetID = [0,1,2,3,4,5]
CPD = [0,1,2,3,4,5]
PMon = [0,1,2,3,4,5]

# DS-0
DetID[0] = {576:1426650, 577:1426650, 592:28469, 593:28469, 594:28465, 595:28465, 598:28470, 599:28470, 600:28463, 601:28463, 608:28455, 609:28455, 610:1425730, 611:1425730, 614:1425381, 615:1425381, 624:1425740, 625:1425740, 626:1426611, 627:1426611, 628:1425742, 629:1425742, 640:1425380, 641:1425380, 642:1426621, 643:1426621, 644:1425741, 645:1425741, 646:28482, 647:28482, 656:1426610, 657:1426610, 662:1425751, 663:1425751, 664:28477, 665:28477, 672:1000020, 674:1426640, 675:1426640, 677:28480, 678:1000012, 679:1000013, 680:1426622, 681:1426622, 688:1426612, 689:1426612, 690:1425750, 691:1425750, 692:1426981, 693:1426981, 696:1425731, 697:1425731}

CPD[0] = {576:123, 577:123, 592:145, 593:145, 594:144, 595:144, 598:142, 599:142, 600:143, 601:143, 608:141, 609:141, 610:134, 611:134, 614:133, 615:133, 624:163, 625:163, 626:162, 627:162, 628:161, 629:161, 640:114, 641:114, 642:173, 643:173, 644:172, 645:172, 646:171, 647:171, 656:153, 657:153, 662:152, 663:152, 664:151, 665:151, 672:2, 674:122, 675:122, 677:131, 678:1, 679:1, 680:124, 681:124, 688:113, 689:113, 690:112, 691:112, 692:111, 693:111, 696:154, 697:154}

PMon[0] = [678, 679, 680, 681, 677, 672]


# DS-1
DetID[1] = {578:1425380, 579:1425380, 580:1426612, 581:1426612, 582:1425750, 583:1425750, 592:1425370, 593:1425370, 594:1426621, 595:1426621, 596:0, 598:1425741, 599:1425741, 600:28482, 601:28482, 608:1425381, 609:1425381, 610:1426980, 611:1426980, 612:0, 614:28469, 615:28469, 616:28480, 617:28480, 624:28455, 625:28455, 626:1425740, 627:1425740, 628:28470, 629:28470, 632:1425742, 633:1425742, 640:1426650, 641:1426650, 644:0, 648:1426640, 649:1426640, 664:1425730, 665:1425730, 672:1426610, 673:1426610, 674:0, 675:0, 676:0, 677:0, 678:1425751, 679:1425751, 690:1426620, 691:1426620, 692:28474, 693:28474, 694:28465, 695:28465}

CPD[1] = {578:114, 579:114, 580:113, 581:113, 582:112, 583:112, 592:174, 593:174, 594:173, 595:173, 596:1, 598:172, 599:172, 600:171, 601:171, 608:133, 609:133, 610:132, 611:132, 612:1, 614:145, 615:145, 616:131, 617:131, 624:141, 625:141, 626:163, 627:163, 628:142, 629:142, 632:161, 633:161, 640:123, 641:123, 644:1, 648:122, 649:122, 664:134, 665:134, 672:153, 673:153, 674:0, 675:0, 676:1, 677:0, 678:152, 679:152, 690:164, 691:164, 692:121, 693:121, 694:144, 695:144}

PMon[1] = [644, 612, 596, 676, 674, 675, 677] # 674,675,677 are not in the MJTChannelMap's due to a bug.


# DS-2
DetID[2] = {578:1425380, 579:1425380, 580:1426612, 581:1426612, 582:1425750, 592:1425370, 593:1425370, 594:1426621, 596:0, 598:1425741, 599:1425741, 600:28482, 601:28482, 608:1425381, 609:1425381, 610:1426980, 611:1426980, 612:0, 616:28480, 617:28480, 626:1425740, 627:1425740, 632:1425742, 633:1425742, 640:1426650, 641:1426650, 644:0, 648:1426640, 649:1426640, 664:1425730, 665:1425730, 672:1426610, 673:1426610, 676:0, 690:1426620, 691:1426620, 692:28474, 693:28474}

CPD[2] = {578:114, 579:114, 580:113, 581:113, 582:112, 592:174, 593:174, 594:173, 596:1, 598:172, 599:172, 600:171, 601:171, 608:133, 609:133, 610:132, 611:132, 612:1, 616:131, 617:131, 626:163, 627:163, 632:161, 633:161, 640:123, 641:123, 644:1, 648:122, 649:122, 664:134, 665:134, 672:153, 673:153, 676:1, 690:164, 691:164, 692:121, 693:121}

PMon[2] = [644, 612, 596, 676]


# DS-3
DetID[3] = {578:1425380, 579:1425380, 580:1426612, 581:1426612, 582:1425750, 583:1425750, 592:1425370, 593:1425370, 594:1426621, 596:0, 598:1425741, 599:1425741, 600:28482, 601:28482, 608:1425381, 609:1425381, 610:1426980, 611:1426980, 612:0, 614:28469, 615:28469, 616:28480, 617:28480, 624:28455, 625:28455, 626:1425740, 627:1425740, 628:28470, 629:28470, 632:1425742, 633:1425742, 640:1426650, 641:1426650, 644:0, 648:1426640, 649:1426640, 664:1425730, 665:1425730, 672:1426610, 673:1426610, 676:0, 678:1425751, 679:1425751, 688:28463, 689:28463, 690:1426620, 691:1426620, 692:28474, 693:28474, 694:28465, 695:28465}

CPD[3] = {578:114, 579:114, 580:113, 581:113, 582:112, 583:112, 592:174, 593:174, 594:173, 596:1, 598:172, 599:172, 600:171, 601:171, 608:133, 609:133, 610:132, 611:132, 612:1, 614:145, 615:145, 616:131, 617:131, 624:141, 625:141, 626:163, 627:163, 628:142, 629:142, 632:161, 633:161, 640:123, 641:123, 644:1, 648:122, 649:122, 664:134, 665:134, 672:153, 673:153, 676:1, 678:152, 679:152, 688:143, 689:143, 690:164, 691:164, 692:121, 693:121, 694:144, 695:144}

PMon[3] = [644, 612, 596, 676]


# DS-4
DetID[4] = {1104:0, 1106:28594, 1107:28594, 1110:1427481, 1111:1427481, 1112:0, 1136:28466, 1137:28466, 1140:28459, 1141:28459, 1142:1426641, 1143:1426641, 1144:28576, 1145:28576, 1170:28607, 1171:28607, 1172:1427491, 1173:1427491, 1174:28481, 1175:28481, 1176:1427490, 1177:1427490, 1200:0, 1204:1427480, 1205:1427480, 1208:28456, 1209:28456, 1232:28717, 1233:28717, 1236:1429090, 1237:1429090, 1238:1427121, 1239:1427121, 1240:0, 1270:0, 1296:1235170, 1297:1235170, 1298:1429091, 1299:1429091, 1300:0, 1302:1427120, 1303:1427120, 1330:28487, 1331:28487, 1332:1428531, 1333:1428531, 1336:0}

CPD[4] = {1104:2, 1106:223, 1107:223, 1110:213, 1111:213, 1112:2, 1136:244, 1137:244, 1140:211, 1141:211, 1142:212, 1143:212, 1144:222, 1145:222, 1170:241, 1171:241, 1172:232, 1173:232, 1174:221, 1175:221, 1176:231, 1177:231, 1200:2, 1204:214, 1205:214, 1208:242, 1209:242, 1232:274, 1233:274, 1236:273, 1237:273, 1238:272, 1239:272, 1240:2, 1270:2, 1296:261, 1297:261, 1298:262, 1299:262, 1300:2, 1302:254, 1303:254, 1330:251, 1331:251, 1332:253, 1333:253, 1336:2}

PMon[4] = [1112, 1104, 1200, 1240, 1270, 1300, 1336]


# DS-5
DetID[5] = {584:1425730, 585:1425730, 592:1426610, 593:1426610, 596:0, 598:1425751, 599:1425751, 608:1425381, 609:1425381, 610:1426980, 611:1426980, 612:0, 614:28469, 615:28469, 616:28480, 617:28480, 624:28455, 625:28455, 626:1425740, 627:1425740, 628:28470, 629:28470, 632:1425742, 633:1425742, 640:1426650, 641:1426650, 644:0, 648:1426640, 649:1426640, 658:1425380, 659:1425380, 660:1426612, 661:1426612, 662:1425750, 663:1425750, 672:1425370, 673:1425370, 674:1426621, 675:1426621, 676:0, 678:1425741, 679:1425741, 680:28482, 681:28482, 688:28463, 689:28463, 690:1426620, 691:1426620, 692:28474, 693:28474, 694:28465, 695:28465, 1104:0, 1106:28594, 1107:28594, 1110:1427481, 1111:1427481, 1112:0, 1120:28466, 1121:28466, 1124:28459, 1125:28459, 1126:1426641, 1127:1426641, 1128:28576, 1129:28576, 1170:28607, 1171:28607, 1172:1427491, 1173:1427491, 1174:28481, 1175:28481, 1176:1427490, 1177:1427490, 1200:0, 1204:1427480, 1205:1427480, 1208:28456, 1209:28456, 1232:28717, 1233:28717, 1236:1429090, 1237:1429090, 1240:0, 1296:1235170, 1297:1235170, 1298:1429091, 1299:1429091, 1300:0, 1302:1427120, 1303:1427120, 1330:28487, 1331:28487, 1332:1428531, 1333:1428531, 1334:0, 1335:0, 1336:0}

CPD[5] = {584:134, 585:134, 592:153, 593:153, 596:1, 598:152, 599:152, 608:133, 609:133, 610:132, 611:132, 612:1, 614:145, 615:145, 616:131, 617:131, 624:141, 625:141, 626:163, 627:163, 628:142, 629:142, 632:161, 633:161, 640:123, 641:123, 644:0, 648:122, 649:122, 658:114, 659:114, 660:113, 661:113, 662:112, 663:112, 672:174, 673:174, 674:173, 675:173, 676:1, 678:172, 679:172, 680:171, 681:171, 688:143, 689:143, 690:164, 691:164, 692:121, 693:121, 694:144, 695:144, 1104:2, 1106:223, 1107:223, 1110:213, 1111:213, 1112:2, 1120:244, 1121:244, 1124:211, 1125:211, 1126:212, 1127:212, 1128:222, 1129:222, 1170:241, 1171:241, 1172:232, 1173:232, 1174:221, 1175:221, 1176:231, 1177:231, 1200:2, 1204:214, 1205:214, 1208:242, 1209:242, 1232:274, 1233:274, 1236:273, 1237:273, 1240:2, 1296:261, 1297:261, 1298:262, 1299:262, 1300:2, 1302:254, 1303:254, 1330:251, 1331:251, 1332:253, 1333:253, 1334:2, 1335:0, 1336:2}

PMon[5] = [612, 676, 596, 1112, 1104, 1200, 1240, 1334, 1300, 1336, 644]

# Copied from DataSetInfo.hh
BkgInfo = [0,1,2,3,4,5,6]
BkgInfo[0] = {0:[2580,2580,2582,2612], 1:[2614,2629,2644,2649,2658,2673],2:[2689,2715], 3:[2717,2750], 4: [2751,2757,2759,2784], 5:[2785,2820], 6:[2821,2855], 7:[2856,2890], 8:[2891,2907,2909,2920], 9: [3137,3166], 10:[3167,3196], 11:[3197,3199,3201,3226], 12:[3227,3256], 13:[3257,3271, 3293,3310], 14:[3311,3340],15:[3341,3370],16:[3371,3372, 3374,3400],17:[3401,3424, 3426,3428, 3431,3432],18:[3461,3462, 3464,3500],19:[3501,3530],20:[3531,3560],21:[3561,3580, 3596,3610],22:[3611,3644],23:[4034,4035, 4038,4039, 4045,4074],24:[4075,4104],25:[4105,4133],26:[4239,4245, 4248,4254, 4256,4268],27:[4270,4271, 4273,4283],28:[4285,4311],29:[4313,4318,4320,4320,4322,4326,4328,4336],30:[4338,4361],31:[4363,4382],32:[4384,4401],33:[4403,4409,4411,4427],34:[4436,4454],35:[4457,4457, 4459,4489],36:[4491,4493, 4573,4573, 4575,4590],37:[4591,4609,4611,4624],38:[4625,4635,4637,4654],39:[4655,4684],40:[4685,4714],41:[4715,4744],42:[4745,4777],43:[4789,4797,4800,4823, 4825,4831],44:[4854,4872],45:[4874,4883, 4885,4907],46:[4938,4945, 4947,4959],47:[4962,4962,4964,4964,4966,4968,4970,4980],48:[5007,5038],49:[5040,5053,5055,5056,5058,5061],50:[5090,5117],51:[5125,5154],52:[5155,5184],53:[5185,5224],54:[5225,5251],55:[5277,5284,5286,5300],56:[5301,5330],57:[5372,5376, 5378,5392, 5405,5414],58:[5449,5458, 5461,5479],59:[5480,5496, 5498,5501, 5525,5526, 5531,5534],60:[5555,5589],61:[5591,5608],62:[5610,5639],63:[5640,5669],64:[5670,5699],65:[5701,5729],66:[5730,5751, 5753,5764],67:[5766,5795],68:[5796,5822],69:[5826,5850],70:[5889,5890, 5894,5895],71:[6553,6575, 6577,6577, 6775,6775],72:[6776,6782, 6784,6809],73:[6811,6830],74:[6834,6853],75:[6887,6903, 6957,6963]}

BkgInfo[1] = {0:[9422,9440],1:[9471,9487, 9492, 9492],2:[9536,9565],3:[9638,9647, 9650, 9668],4:[9674,9676, 9678, 9678, 9711, 9727],5:[9763,9780],6:[9815,9821, 9823, 9832, 9848, 9849, 9851, 9854],7:[9856,9912],8:[9928,9928],9:[9952,9966, 10019, 10035],10:[10074,10090, 10114, 10125],11:[10129,10149],12:[10150,10171],13:[10173,10203],14:[10204,10231],15:[10262,10278, 10298,10299, 10301,10301, 10304,10308],16:[10312,10342],17:[10344,10350, 10378,10394, 10552,10558],18:[10608,10648],19:[10651,10677],20:[10679,10717],21:[10745,10761, 10788,10803],22:[10830,10845, 10963,10976],23:[11002,11008, 11010,11019, 11046,11066],24:[11083,11113],25:[11114,11144],26:[11145,11175],27:[11176,11200, 11403,11410],28:[11414,11417, 11419,11426, 11428,11432, 11434,11444, 11446,11451],29:[11453,11453, 11455,11458, 11466,11476, 11477,11483],30:[12521,12522, 12525,12526, 12528,12537, 12539,12539, 12541,12543, 12545,12547, 12549,12550],31:[12551,12551, 12553,12560, 12562,12575, 12577,12578, 12580,12580],32:[12607,12625, 12636,12647, 12652,12653],33:[12664,12675],34:[12677,12695, 12697,12724],35:[12736,12765],36:[12766,12798],37:[12816,12816, 12819,12819, 12824,12824, 12827,12827, 12829,12831, 12834,12838, 12842,12842, 12843,12861, 12875,12875],38:[13000,13003, 13005,13028],39:[13029,13053, 13055,13056],40:[13066,13070, 13076,13092, 13094,13096],41:[13100,13115, 13117,13119, 13123,13137],42:[13148,13150, 13154,13156, 13186,13189, 13191,13204, 13206,13211],43:[13212,13242],44:[13243,13275],45:[13276,13287, 13306,13311, 13313,13325],46:[13326,13350, 13362,13368],47:[13369,13383, 13396,13411],48:[13519,13548],49:[13699,13704, 13715,13719],50:[14010,14040, 14041,14041],51:[14342,14372, 14386,14387]}

BkgInfo[2] = {0:[14775,14786, 14788,14805],1:[14908,14925, 14936,14941, 14943,14948],2:[15043,15052, 15062,15083],3:[15188,15188, 15190,15193, 15195,15218],4:[15324,15326, 15338,15338, 15343,15364],5:[15471,15483, 15511,15519, 15613,15621, 15625,15625],6:[15635,15657],7:[15763,15767, 15769,15787, 15797,15803]}

BkgInfo[3] = {0:[16797,16826, 16827,16835],1:[16857,16886],2:[16887,16910, 16931,16935, 16947,16952],3:[16957,16959, 16970,16999],4:[17000,17009, 17035,17057],5:[17060,17090],6:[17091,17121],7:[17122,17127, 17129,17131, 17138,17156],8:[17159,17181, 17305,17318],9:[17322,17343],10:[17351,17381],11:[17382,17412, 17413,17422],12:[17448,17477],13:[17478,17493],14:[17500,17519],15:[17531,17553, 17555,17559],16:[17567,17597],17:[17598,17628],18:[17629,17659],19:[17660,17686],20:[17703,17717, 17720,17721],21:[17852,17882],22:[17883,17913],23:[17914,17944],24:[17945,17948, 17967,17980]}

BkgInfo[4] = {0:[60000802,60000821, 60000823,60000823, 60000827,60000828, 60000830,60000830],1:[60000970,60001000],2:[60001001,60001010],3:[60001033,60001054, 60001056,60001062],4:[60001063,60001086],5:[60001088,60001093],6:[60001094,60001124],7:[60001125,60001125, 60001163,60001181, 60001183,60001185],8:[60001187,60001205, 60001309,60001319],9:[60001331,60001350, 60001380,60001382],10:[60001384,60001414],11:[60001415,60001441],12:[60001463,60001489],13:[60001491,60001506],14:[60001523,60001542],15:[60001597,60001624],16:[60001625,60001655],17:[60001656,60001686],18:[60001687,60001714]}

BkgInfo[5] = {0:[18623,18624, 18628,18629, 18645,18652],1:[18654,18685],2:[18686,18703, 18707,18707],3:[18761,18783],4:[18808,18834],5:[18835,18838, 18844,18844, 18883,18914],6:[18915,18918, 18920,18951, 18952,18957],7:[19240,19240, 19264,19280, 19305,19318],8:[19320,19351],9:[19352,19383],10:[19384,19385, 19387,19415],11:[19416,19425, 19428,19430, 19436,19445],12:[19481,19496, 19502,19515],13:[19613,19644],14:[19645,19676],15:[19677,19677, 19696,19697, 19707,19722],16:[19733,19747, 19771,19773],17:[19775,19801, 19806,19806],18:[19834,19860],19:[19862,19893],20:[19894,19899, 19901,19907],21:[19968,19998],22:[19999,19999, 20021,20040],23:[20074,20105],24:[20106,20130, 20132,20134],25:[20136,20167],26:[20168,20199],27:[20218,20237],28:[20239,20270],29:[20271,20286, 20311,20316, 20319,20332],30:[20335,20365],31:[20366,20375, 20377,20397],32:[20398,20415],33:[20417,20445],34:[20483,20487, 20489,20491, 20494,20509],35:[20522,20537],36:[20611,20629, 20686,20691],37:[20755,20756, 20758,20786],38:[20787,20795, 20797,20828],39:[20829,20860],40:[20861,20876, 20877,20882],41:[20884,20915],42:[20916,20927,20929,20957],43:[20964,20995],44:[20996,21012],45:[21014,21045],46:[21046,21058],47:[21060,21091],48:[21092,21104],49:[21106,21136],50:[21158,21167,21169,21178,21201,21201],51:[21217,21248],52:[21249,21278],53:[21280,21311],54:[21312,21343],55:[21344,21375],56:[21376,21389,21391,21407],57:[21408,21424,21426,21435,21452,21453],58:[21469,21499],59:[21501,21532],60:[21533,21564],61:[21565,21585, 21587,21587],62:[21595,21614, 21617,21618, 21622,21628],63:[21630,21661],64:[21662,21674, 21691,21692, 21694,21705],65:[21747,21776],66:[21778,21800, 21833,21837],67:[21839,21853, 21856,21857, 21862,21879],68:[21891,21893, 21895,21908, 21922,21937],69:[21940,21940, 21953,21968],70:[22001,22032],71:[22033,22064],72:[22065,22095],73:[22097,22100, 22102,22122],74:[22127,22142],75:[22147,22171, 22173,22176],76:[22180,22213],77:[22214,22247],78:[22248,22250, 22266,22280,22304,22304,22316,22333],79:[22340,22356,22369,22392],80:[22400,22428],81:[22430,22463],82:[22464,22488],83:[22490,22512, 22636,22644, 22647,22650],84:[22652,22653, 22655,22670, 22673,22674],85:[22678,22711],86:[22712,22742],87:[22744,22750, 22753,22755, 22760,22763, 22765,22777, 22814,22815],88:[22817,22834,22838,22838,22840,22840,22853,22853,22867,22867],89:[22876,22909],90:[22910,22943],91:[22944,22946, 22952,22952, 22954,22954, 22959,22982, 22984,22986],92:[22993,22996, 23085, 23101],93:[23111,23144],94:[23145,23175,23211,23212],95:[23218,23232,23246,23260,23262,23262],96:[23282,23306],97:[23308,23334],98:[23338,23370],99:[23372,23405],100:[23406,23433],101:[23440,23458, 23461,23462, 23469,23480],102:[23511,23513,23520,23521,23525,23542,23548,23548],103:[23551,23584],104:[23585,23618],105:[23619,23642],106:[23645,23668,23675,23690],107:[23704,23715,23718,23719,23721,23721],108:[23725,23758],109:[23759,23792],110:[23793,23826],111:[23827,23849, 23851,23867],112:[23869,23881, 23939,23940, 23942,23958]}

BkgInfo[6] = {0:[25704,25725,  25728,25737],1:[25738,25756,  25759,25763,  25765,25771],2:[25772,25787,  25790,25800],3:[25801,25819,  25822,25830],4:[26023,26034,  26036,26038,  26052,26066],5:[26163,26169,  26171,26176,  26179,26190],6:[26191,26192,  26194,26194,  26313,26325,  26328,26344],7:[26465,26490,  26493,26495],8:[26591,26601,  26603,26616,  26617,26617,  26619,26622],9:[26742,26745,  26773,26773,  26775,26776,  26780,26789,  26791,26805],10:[26907,26918,  26920,26938],11:[27060,27070,  27074,27091],12:[27217,27224,  27227,27248],13:[27920,27922,  27924,27930,  27932,27936,  27958,27969],14:[27991,28013,  28015,28018],15:[28019,28035,  28037,28045,  28048,28049],16:[28050,28072,  28074,28076,  28079,28080],17:[28081,28108,  28111,28111],18:[28136,28160],19:[28161,28161,  28164,28169,  28171,28186],20:[28300,28320],21:[28391,28402],22:[28403,28406,  28409,28422,  28662,28673],23:[28674,28683,  28685,28688,  28690,28693,  28813,28816,  28818,28824],24:[28825,28837,  28840,28844,  28942,28955],25:[28956,28964,  28967,28967,  28992,28997],26:[29092,29105,  29109,29122,  29124,29124]}

def GetBkgIdx(dsNum, runNum):
    """ Not completely accurate way of finding the background idx """
    #TODO: Add way of finding missing runs in a subDS (does it even matter?)
    bkgidx = [key for key in BkgInfo[dsNum] if runNum <= BkgInfo[dsNum][key][-1] and runNum >= BkgInfo[dsNum][key][0]]
    try:
        return bkgidx[0]
    except:
        print "Run %d not found in Dataset, returning -1"%(runNum)
        return -1


# ==================================================================================
#                          VETO & BAD DETECTOR INFORMATION
# The first two functions, 'LoadBadDetectorMap' and 'LoadVetoDetectorMap', were
# taken from DataSetInfo.hh from the GAT version on June 21 2017.
# ==================================================================================

def LoadBadDetectorMap(dsNum):
    detIDIsBad = []
    if dsNum==0: detIDIsBad = [28474, 1426622, 28480, 1426980, 1426620, 1425370]
    if dsNum==1: detIDIsBad = [1426981, 1426622, 28455, 28470, 28463, 28465, 28469, 28477, 1425751, 1425731, 1426611]
    if dsNum==2: detIDIsBad = [1426981, 1426622, 28455, 28470, 28463, 28465, 28469, 28477, 1425731, 1426611]
    if dsNum==3: detIDIsBad = [1426981, 1426622, 28477, 1425731, 1426611]
    if dsNum==4: detIDIsBad = [28595, 28461, 1428530, 28621, 28473, 1426651, 1429092, 1426652, 28619]
    if dsNum==5: detIDIsBad = [1426981, 1426622, 28477, 1425731, 1426611, 28595, 28461, 1428530, 28621, 28473, 1426651, 1429092, 1426652, 28619, 1427121]
    if dsNum==6: detIDIsBad = [1426981, 28474, 1426622, 28477, 1425731, 1426611, 28595, 28461, 1428530, 28621, 28473, 1426651, 1429092, 1426652, 28619, 1427121]
    return detIDIsBad


def LoadVetoDetectorMap(dsNum):
    detIDIsVetoOnly = []
    if dsNum == 0: detIDIsVetoOnly = [1426621, 1425381, 1425742]
    if dsNum == 1: detIDIsVetoOnly = [1426621, 28480]
    if dsNum == 2: detIDIsVetoOnly = [1426621, 28480, 1425751]
    if dsNum == 3: detIDIsVetoOnly = [1426621, 28480, 28470, 28463]
    if dsNum == 4: detIDIsVetoOnly = [28459, 1426641, 1427481, 28456, 1427120, 1427121]
    if dsNum == 5: detIDIsVetoOnly = [1426621, 28480, 1426641, 1235170, 1429090]
    if dsNum == 6: detIDIsVetoOnly = [1426621, 28480, 1426641, 1235170, 1429090]
    return detIDIsVetoOnly


def GetGoodChanList(dsNum):
    badIDs = LoadBadDetectorMap(dsNum) + LoadVetoDetectorMap(dsNum)

    # make a list of the channels corresponding to the bad IDs.
    badChans = []
    for badID in badIDs:
        for ch, detID in DetID[dsNum].iteritems():
            if badID == detID: badChans.append(ch)
    # print sorted(badChans)

    # high-gain channels, without pulser monitors, without bad+veto channels.
    goodList = [key for key in DetID[dsNum] if key%2==0 and key not in PMon[dsNum] and key not in badChans]
    # print sorted(goodList)
    return sorted(goodList)



def GetThreshDicts(dsNum, threshCut=0.9):
    #TODO: Update with DB once DB parameters are implemented!
    import numpy as np
    import pandas as pd
    import os

    # homePath = os.path.expanduser('~') # glob doesn't know how to expand this
    homePath = "/Users/brianzhu"
    inDir = homePath + "/project/thresholds/"
    df = pd.read_hdf(inDir+'ThreshDS%d_Processed.h5' % dsNum, 'threshTree')

    goodRuns, badRuns, goodRunErfs = {}, {}, {}
    for column in df:
        col = int(column)
        goodRuns[col], badRuns[col], goodRunErfs[col] = [], [], []

        for idx, vals in enumerate(df.loc[:,column]):

            # skip NaN run ranges where data wasn't collected for the channel
            if np.isnan(vals).any(): continue

            thresh, sigma, hi, lo = vals[0], vals[1], int(vals[2]), int(vals[3])

            if vals[0] <= threshCut:
                goodRuns[col].append([hi,lo])
                goodRunErfs[col].append([hi,lo,thresh,sigma])
            else:
                badRuns[col].append([hi,lo])

    return goodRuns,badRuns,goodRunErfs
