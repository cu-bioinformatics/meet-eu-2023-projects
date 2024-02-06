import re
import pandas as pd
import os
import sys

path=sys.argv[1]
info={"header":{},"df":{}}
with open(path) as f:
    table=False
    n=0
    
    for line in f.readlines():
        line=line.strip()
        if line.startswith("mode"):
            header=[l.strip() for l in line.split("|")]
            info["header"][0]=header[0]
            info["header"][1]=header[1]
        elif line.startswith("|"):
            header=[l.strip() for l in line.split("|")]
            info["header"][2]=header[2]
            info["header"][3]=header[3]
        elif line.startswith("-"):
            table=True
        elif line.startswith("Writing output"):
            table=False

        if table==True and not line.startswith("-"):
            elements=line.split()
            info["df"][n]={info["header"][i]:elements[i] for i in range(len(elements))}
            n=n+1
            
df=pd.DataFrame(info["df"]).T
new_filename="results.csv"
new_path=os.path.dirname(path)+"/"+new_filename
df.to_csv(new_path, index = False, header=True)
