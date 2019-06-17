import pandas as pd 
import time
import numpy as np
start = time.time()
print("Starting process!")
df = pd.read_csv("500-PSSS3-equ decoy_Report.xls", sep = "\t", usecols = [1,2,3,6,7,9,12,14])

proteinId = df["PG.ProteinAccessions"].apply(lambda x: np.nan if ";" in x else x)
df["proteinId"] = proteinId
df = df[pd.notnull(df["proteinId"])]
proteinId = df["PG.ProteinAccessions"]

reformatted_df = pd.DataFrame()
reformatted_df["qvalue"] = df["PG.Qvalue"]
df.drop(["PG.Qvalue"], inplace = True, axis = 1)
reformatted_df["proteinQuantity"] = df["PG.Quantity"]
df.drop(["PG.Quantity"], inplace = True, axis = 1)
reformatted_df["decoy"] = df["EG.IsDecoy"]
df.drop(["EG.IsDecoy"], inplace = True, axis = 1)
specie = df["PG.Organisms"].apply(lambda x: "CAEEL" if x == "Caenorhabditis elegans OX=6239" else "ARATH" if x == "Arabidopsis thaliana" else "HUMAN" if x == "Homo sapiens" else x)
df.drop("PG.Organisms", inplace = True, axis = 1)
reformatted_df["peptideIntensity"] = df["FG.NormalizedMS2PeakArea"]
df.drop("FG.NormalizedMS2PeakArea", inplace = True, axis = 1)
reformatted_df["peptideQvalue"] = df["EG.Qvalue"]

samples = df["R.FileName"].str[-15:-12]
runs = df["R.FileName"].str[-6:-3]
df.drop("R.FileName", inplace = True, axis = 1)
del df
run_id = samples + ":" + samples + "_" + runs
run_id = run_id.apply(lambda x: "S03:S04_R05" if x == "MHR_MHR_X01" else x)
#reformatted_df["run_id"] = run_id
reformatted_df["sample"] = run_id.str[:3]
reformatted_df["run"] = run_id.str[-3]
reformatted_df["protein"] = proteinId
reformatted_df["specie"] = specie
reformatted_df.to_csv("rawPSSS3_relevantCols.csv", sep = "\t", index = False)
end = time.time()
print(end-start)
print("Done!")



