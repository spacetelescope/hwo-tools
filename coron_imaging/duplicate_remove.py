import sys
import numpy as np
import pandas as pd

filename = sys.argv[1]
outfilename = filename.split(".")[0] + "_dedup.txt"

df = pd.read_fwf(filename, colspecs="infer", header=None)
df.info()
print(df)
duplicate = df.duplicated(subset=0, keep="first")
print(duplicate[duplicate == True])
print(df[duplicate == True])

#duplicate = duplicate[duplicate == True]

for entry in df[duplicate == True].index:
    print(df.iloc[entry-1][2])
    print(df.iloc[entry][2])
    value = np.mean((df.iloc[entry-1][2], df.iloc[entry][2]))
    print(entry, value)
    df.iloc[entry][2] = value

dedup = df.drop_duplicates(subset=0)

dedup.to_csv(outfilename, header=False, index=False)
