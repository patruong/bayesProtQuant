import csv
import sys

from triqler import parsers

def main(argv):
  inputFile = argv[0]
  reader = csv.reader(open(inputFile, 'r'), delimiter = '\t')
  header = reader.next()
  
  for i, colName in enumerate(header):
    print i, colName
  
  '''
  0 R.Condition
  1 R.FileName
  2 PG.Organisms
  3 PG.ProteinAccessions
  4 PG.Cscore
  5 PG.NrOfStrippedSequencesIdentified
  6 PG.Qvalue
  7 EG.IsDecoy
  8 EG.PrecursorId
  9 EG.PEP
  10 EG.Qvalue
  11 EG.Cscore
  12 FG.NormalizedMS2PeakArea
  '''
  
  writer = csv.writer(open(argv[1], 'w'), delimiter = '\t')
  writer.writerow(parsers.TriqlerSimpleInputRowHeaders)
  
  conditionCol = header.index("R.Condition")
  if "PSSS3" in inputFile:
    runCol = header.index("R.FileName")
  else:
    runCol = header.index("Measurement")
  peptideCol = header.index("EG.PrecursorId")
  decoyCol = header.index("EG.IsDecoy")
  organismCol = header.index("PG.Organisms")
  proteinCol = header.index("PG.ProteinAccessions")
  scoreCol = header.index("EG.Cscore")
  intensityCol = header.index("FG.NormalizedMS2PeakArea")
    
  for spectrumIdx, row in enumerate(reader):
    peptide = row[peptideCol][1:-3].replace("C[+C2+H3+N+O]", "C").replace("M[+O]", "M[16]").replace("[+C2+H2+O]", "[42]")
    
    proteinSuffix = ""
    if "PSSS3" in inputFile:
      condition = row[conditionCol][-3:]
      run = row[runCol][21:24] + row[runCol][29:33]
      organism = row[organismCol]
      if organism == "Arabidopsis thaliana":
        proteinSuffix = "_ARATH"
      elif organism == "Caenorhabditis elegans OX=6239":
        proteinSuffix = "_CAEEL"
      elif organism == "Homo sapiens":
        proteinSuffix = "_HUMAN"
    elif "PSSS2" in inputFile:
      run = row[runCol][28:32]
      condition = row[conditionCol][-2:]
    elif "GSDS" in inputFile:
      condition = row[conditionCol][-2:]
      run = row[runCol][20:23] + row[runCol][28:31]
      organism = row[organismCol]
      if organism == "Bos taurus":
        proteinSuffix = "_BOVIN"
      elif organism == "Equus caballus":
        proteinSuffix = "_CAEEL"
      elif organism == "Gallus gallus":
        proteinSuffix = "_HORSE"
      elif organism == "Homo sapiens":
        proteinSuffix = "_HUMAN"
    else:
      run = row[runCol]
      condition = row[conditionCol]
    
    proteinPrefix = ""
    if row[decoyCol] == "True":
      proteinPrefix = "decoy_"
    
    proteins = [proteinPrefix + x + proteinSuffix for x in row[proteinCol].split(";")]
    if len(proteins) == 1:
      tir = parsers.TriqlerInputRow(run, condition, int(row[peptideCol][-1]), spectrumIdx, 0.0, -1, row[scoreCol], row[intensityCol], peptide, proteins)
      writer.writerow(tir.toSimpleList())

if __name__ == "__main__":
  main(sys.argv[1:])
