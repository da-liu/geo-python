import urllib, os, gzip, csv, re

DOWNLOAD_DIR = "./tmp/"
OUTPUT_DIR   = "./output/"
ORDERED_KEYS = ["title","geo_accession","status","submission_date","last_update_date","pubmed_id","summary","overall_design","type","contributor","sample_id","contact_name","contact_email","contact_web_link","contact_phone","contact_fax","contact_laboratory","contact_department","contact_institute","contact_address","contact_city","contact_state","contact_zip/postal_code","contact_country","web_link","supplementary_file","platform_id","platform_taxid","sample_taxid","relation"]

def gzMatrixName(gsename):
  return gsename + "_series_matrix.txt.gz"

def txtMatrixName(gsename):
  return gsename + "_series_matrix.txt"

def rawName(gsename):
  return gsename + "_RAW.tar"

def gzMatrixPath(gsename):
  return DOWNLOAD_DIR + gsename + "/" + gzMatrixName(gsename)

def txtMatrixPath(gsename):
  return DOWNLOAD_DIR + gsename + "/" + txtMatrixName(gsename)

def subUrl(gsename):
  '''
  Convert "GSE74432" => "GSE74nnn"
  Convert "GSE13" => "GSEnnn"
  '''
  return gsename[:3] + "nnn" if len(gsename[3:]) <= 3 else gsename[:-3] + "nnn"

def matrixUrl(gsename):
  return "https://ftp.ncbi.nlm.nih.gov/geo/series/" + subUrl(gsename) + "/" + gsename + "/matrix/" + gzMatrixName(gsename)

def suppUrl(gsename):
  return "https://ftp.ncbi.nlm.nih.gov/geo/series/{0}/{1}/suppl/".format(subUrl(gsename), gsename)

def suppFilelist(gsename):
  return suppUrl(gsename) + "filelist.txt"

def downloadMatrix(gsename):
  downloadPath = DOWNLOAD_DIR + gsename + "/"
  if not os.path.exists(downloadPath):
    os.makedirs(downloadPath)
  filepath = downloadPath + gzMatrixName(gsename)
  if not os.path.isfile(gzMatrixPath(gsename)) or not isValidGZ(gzMatrixPath(gsename)):
    print "Downloading matrix file:", filepath
    urllib.urlretrieve(matrixUrl(gsename), filepath)

def hasIDAT(*gsenames):
  res = []
  for gsename in gsenames:
    downloadPath = DOWNLOAD_DIR + gsename + "/suppl/"
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
    downloadPath = DOWNLOAD_DIR + gsename + "/suppl/"
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
  # consider cases of duplicate keys
  # k = "ftp", v = ['link1', 'link2', ...]
  # k = "ftp", v = ['link1', 'link2', ...]
  # k = "ftp", v = ['link1', 'link2', ...]
  # ...
  # We want to convert them to 
  # k = "ftp", v = ['link1', 'link2', ...]
  # k = "ftp1", v = ['link1', 'link2', ...]
  # k = "ftp2", v = ['link1', 'link2', ...]
  # ...
  if key in dic:
    index = 1
    newKey = key + str(index)
    while newKey in dic:
      index += 1
      newKey = key + str(index)
    dic[newKey] = val
  else:
    dic[key] = val

def parseHeader(l, header):
  key = getKey(l)
  val = l[1] if len(l) > 1 else ""
  concatDicValues(key, val, header)

def parseSampleTable(l, sampleTable):
  key = getKey(l)
  if len(l) > 1:
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
  else:
    val = ""
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
      # print allKeys(dic) - set(ORDERED_KEYS)
      for key in ORDERED_KEYS:
        if key in dic: 
          row = []
          row.append(key)
          row.extend(dic[key])
          writer.writerow(row)
      sortedKeys = sorted(allKeys(dic) - set(ORDERED_KEYS))
      for key in sortedKeys:
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

def cleanAttributes(fname, *sampleTables):
  listOfDics = []
  for i, table in enumerate(sampleTables):
    dic = {}
    for k, v in table.iteritems():
      if len(set(v)) == len(v):
        if v == "":
          v = '...'
        else:
          v = '"' + v[0] + '", ...'
      else:
        v = "; ".join(set(v))
      dic[k] = v
    listOfDics.append(dic)
  return listOfDics

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
  writeCSV(OUTPUT_DIR + "headers.csv", headers)
  
  print "Sample Headers ..."
  writeCSV(OUTPUT_DIR + "sampleHeaders.csv", headers)
  
  print "Sample Tables ..."
  writeSampleTables(OUTPUT_DIR + "sampleTables.csv", *sampleTables)

  print "Attributes ..."
  attrTables = cleanAttributes(OUTPUT_DIR + "attributes.csv", *sampleTables)
  attrTables = mergeDics(*attrTables)
  writeCSV(OUTPUT_DIR + "attrTables.csv", attrTables)


def openGSEList(fname):
  return [line.translate(None, '\n"') for line in open(fname)]




import time


# gsenames = ["GSE74432","GSE74486","GSE64380","GSE41273","GSE75248","GSE68747","GSE72120","GSE69270","GSE57361","GSE51921","GSE61431","GSE59685","GSE50586","GSE61107","GSE53191","GSE52588","GSE63347","GSE74193","GSE41169"]
# gsenames = ["GSE74432","GSE74486","GSE64380","GSE41273","GSE75248","GSE68747","GSE72120","GSE69270","GSE51921","GSE61431","GSE59685","GSE50586","GSE61107","GSE53191","GSE52588","GSE63347","GSE74193","GSE41169"]
# gsenames = ['GSE74432', 'GSE75248', 'GSE72120', 'GSE61107', 'GSE74193'] # list with IDAT
# gsenames = ["GSE64380", "GSE68747", "GSE74486"]
# gsenames = ["GSE64380", "GSE68747"]
# gsenames = ["GSE57361"] # two GPLs, faulty matrix
# gsenames = ["GSE57361","GSE59685"] # two GPLs

# gsenames = ["GSE59685"]
# gsenames = ["GSE74432"]

gsenames = openGSEList("random100.txt")


# gsenames = gsenames[13:14]
start = time.time()


parseMatrix(gsenames)

# print hasIDAT(*gsenames)
# map(downloadRAW, gsenames)


end = time.time()
print(end - start)




