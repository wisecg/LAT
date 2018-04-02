import imp
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')

cal = ds.CalInfo()
# print cal.GetSpecialRuns("extPulser",1)
# print cal.GetSpecialRuns("extPulser")
# print cal.GetSpecialNIdxs("extPulser")
# print cal.GetSpecialKeys()
# for key in cal.GetSpecialKeys():
    # print key, cal.GetSpecialNIdxs(key)
