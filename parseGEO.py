import urllib, os, gzip, csv, re

gseDir = "./tmp/"
outDir = "./output/"

def gzMatrixName(gsename):
  return gsename + "_series_matrix.txt.gz"

def txtMatrixName(gsename):
  return gsename + "_series_matrix.txt"

def rawName(gsename):
  return gsename + "_RAW.tar"

def gzMatrixPath(gsename):
  return gseDir + gsename + "/" + gzMatrixName(gsename)

def txtMatrixPath(gsename):
  return gseDir + gsename + "/" + txtMatrixName(gsename)

def matrixUrl(gsename):
  subUrl = gsename[:-3] + "nnn"
  return "https://ftp.ncbi.nlm.nih.gov/geo/series/" + subUrl + "/" + gsename + "/matrix/" + gzMatrixName(gsename)

def suppUrl(gsename):
  subUrl = gsename[:-3] + "nnn"
  return "https://ftp.ncbi.nlm.nih.gov/geo/series/{0}/{1}/suppl/".format(subUrl, gsename)

def suppFilelist(gsename):
  return suppUrl(gsename) + "filelist.txt"

orderedKeys = ["title","geo_accession","status","submission_date","last_update_date","pubmed_id","summary","overall_design","type","contributor","sample_id","contact_name","contact_email","contact_web_link","contact_phone","contact_fax","contact_laboratory","contact_department","contact_institute","contact_address","contact_city","contact_state","contact_zip/postal_code","contact_country","web_link","supplementary_file","platform_id","platform_taxid","sample_taxid","relation"]

def downloadMatrix(gsename):
  downloadPath = gseDir + gsename + "/"
  if not os.path.exists(downloadPath):
    os.makedirs(downloadPath)
  print "Downloading RAW file:", filepath
  urllib.urlretrieve(matrixUrl(gsename), downloadPath + gzMatrixName(gsename))

def hasIDAT(*gsenames):
  res = []
  for gsename in gsenames:
    downloadPath = gseDir + gsename + "/suppl/"
    if not os.path.exists(downloadPath):
      os.makedirs(downloadPath)
    filepath = downloadPath + "filelist.txt"
    if not os.path.isfile(filepath):
      url = suppFilelist(gsename)
      urllib.urlretrieve(url, filepath)
    # idat files exist if the '.idat' extension can be found in filelist.txt
    # Alternatively 'IDAT', 'Grn.idat', 'Red.idat' are also markers for idat
    res.append(gsename if '.idat' in open(filepath).read() else False)
  return res[0] if len(res) == 1 else res

def downloadRAW(gsename):
  if hasIDAT(gsename):
    downloadPath = gseDir + gsename + "/suppl/"
    if not os.path.exists(downloadPath):
      os.makedirs(downloadPath)
    url = suppUrl(gsename) + rawName(gsename)
    filepath = downloadPath + rawName(gsename)
    print "Downloading RAW file:", filepath
    urllib.urlretrieve(url, filepath)

def isValidGZ(fname):
  with gzip.open(fname) as f:
    try:
      f.readline()
      return True
    except:
      print "Cannot read file: ", fname
      return False

def getKey(l):
  return l[0][l[0].find("_")+1:]

def concatDicValues(key, val, dic):
  if key in dic:
    dic[key] += ";" + val
  else:
    dic[key] = val

def appendDicValues(key, val, dic):
  if key in dic:
    dic[key].append(val)
  else:
    dic[key] = val

def parseHeader(l, header):
  key = getKey(l)
  val = l[1] if len(l) > 1 else ""
  concatDicValues(key, val, header)

def parseSampleTable(l, sampleTable):
  key = getKey(l)
  if l[1].find(":") != -1:
    rk = l[1][:l[1].find(":")]
    rkIsCommon = True
    for e in l[1:]:
      if not e.startswith(rk):
        rkIsCommon = False
    if rkIsCommon:
      key = rk
      val = [x[l[1].find(":")+2:] for x in l[1:]]
    else:
      val = l[1:]
  else:
    val = l[1:]
  appendDicValues(key, val, sampleTable)

def parseGEO(f):
  header = {}
  sampleHeader = {}
  sampleTable = {}
  r = csv.reader(f, delimiter = '\t', quotechar = '"')
  for l in r:
    l = [x for x in l if x != '']
    if l == []:
      continue
    elif l[0].startswith("!Series_"):
      parseHeader(l, header)
    elif l[0].startswith("!Sample_"):
      if len(set(l[1:])) == 1:
        parseHeader(l, sampleHeader)
      else:
        parseSampleTable(l, sampleTable)
    else:
      return header, sampleHeader, sampleTable

def readCSVgz(fname):
  with gzip.open(fname) as f:
    return parseGEO(f)

def readCSVtxt(fname):
  with open(fname) as f:
    return parseGEO(f)

def allKeys(*dics):
  keys = []
  for d in dics:
    keys += d.keys()
  return set(keys)

def writeCSV(fname, dic):
  with open(fname, 'wb') as csv_file:
      writer = csv.writer(csv_file)
      # print allKeys(dic) - set(orderedKeys)
      for key in orderedKeys:
        if key in dic: 
          row = []
          row.append(key)
          row.extend(dic[key])
          writer.writerow(row)

def writeSampleTables(fname, *sampleTables):
  with open(fname, 'wb') as f:
    w = csv.writer(f)
    for i, table in enumerate(sampleTables):
      w.writerow([gsenames[i]])
      for k, v in table.iteritems():
        row = []
        row.append(k)
        row.extend(v)
        w.writerow(row)
      w.writerow("")
      w.writerow("")

def mergeDics(*dics):
  keys = allKeys(*dics)
  dic = {}
  for k in keys:
    dic[k] = []
  for d in dics:
    for k in keys:
      if k in d:
        dic[k].append(d[k])
      else:
        dic[k].append("")
  return dic


def parseMatrix(gsenames):
  headers = []
  sampleHeaders = []
  sampleTables = []
  for i in xrange(len(gsenames)):
    if not os.path.isfile(gzMatrixPath(gsenames[i])):
      downloadMatrix(gsenames[i])
    if isValidGZ(gzMatrixPath(gsenames[i])):
      a, b, c = readCSVgz(gzMatrixPath(gsenames[i]))
      headers.append(a)
      sampleHeaders.append(b)
      sampleTables.append(c)
  headers = mergeDics(*headers)
  sampleHeaders = mergeDics(*sampleHeaders)

  print "Writing Series Headers ..."
  writeCSV(outDir + "headers.csv", headers)
  print "Sample Headers ..."
  writeCSV(outDir + "sampleHeaders.csv", headers)
  print "Sample Tables ..."
  writeSampleTables(outDir + "sampleTables.csv", *sampleTables)



import time


# gsenames = ["GSE74432","GSE74486","GSE64380","GSE41273","GSE75248","GSE68747","GSE72120","GSE69270","GSE57361","GSE51921","GSE61431","GSE59685","GSE50586","GSE61107","GSE53191","GSE52588","GSE63347","GSE74193","GSE41169"]
# gsenames = ["GSE74432","GSE74486","GSE64380","GSE41273","GSE75248","GSE68747","GSE72120","GSE69270","GSE51921","GSE61431","GSE59685","GSE50586","GSE61107","GSE53191","GSE52588","GSE63347","GSE74193","GSE41169"]
gsenames = ['GSE74432', 'GSE75248', 'GSE72120', 'GSE61107', 'GSE74193'] # list with IDAT
# gsenames = ["GSE64380", "GSE68747", "GSE74486"]
# gsenames = ["GSE64380", "GSE68747"]
# gsenames = ["GSE57361"] # two GPLs, faulty matrix
# gsenames = ["GSE57361","GSE59685"] # two GPLs

# gsenames = ["GSE59685"]
# gsenames = ["GSE74432"]


start = time.time()

# parseMatrix(gsenames)

print hasIDAT(*gsenames)
# map(downloadRAW, gsenames)



end = time.time()
print(end - start)




