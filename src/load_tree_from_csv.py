import pandas as pd

column_dtypes = {}

column_dtypes['ID']= 'Int32'
column_dtypes['nodeType']= 'Int8'
column_dtypes['nrChildren']= 'Int8'
column_dtypes['IDfirstChild']= 'Int32'
column_dtypes['split']= 'Int8'
column_dtypes['T1_min']= 'Int32'
column_dtypes['T1_max']= 'Int32'
column_dtypes['V1_min']= 'Int32'
column_dtypes['V1_max']= 'Int32'
column_dtypes['T2_min']= 'Int32'
column_dtypes['T2_max']= 'Int32'
column_dtypes['V2_min']= 'Int32'
column_dtypes['V2_max']= 'Int32'
column_dtypes['A_min']= 'Int32'
column_dtypes['A_max']= 'Int32'
column_dtypes['P1_index']="Int8"
column_dtypes['P2_index']="Int8"
column_dtypes['P3_index']="Int8"
column_dtypes['Q1_index']="Int8"
column_dtypes['Q2_index']="Int8"
column_dtypes['Q3_index']="Int8"
column_dtypes['r']="Int32"
column_dtypes['sigma_Q']="Int8"
column_dtypes['wx_nominator']="Int64"
column_dtypes['wy_nominator']="Int64"
column_dtypes['w_denominator']="Int64"
column_dtypes['S_index']="Int8"

chunksize=10000
chunkdata=[]

for chunk in pd.read_csv("../data/solution_tree.csv", sep=',',  dtype=column_dtypes, chunksize=int(chunksize)):
    chunkdata.append(chunk)
    
TREE=pd.concat(chunkdata,ignore_index=True)
del(chunkdata)

# Save as parquet:
#TREE.to_parquet(
#    "../data/solution_tree.parquet",
#    engine="pyarrow",
#    compression="brotli"
#)










