"""One-million-row throughput benchmark using resampled public NHANES inputs.
This is a computational stress test, NOT one million independent participants.
"""
from pathlib import Path
import json,time,resource,platform
import numpy as np,pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import HistGradientBoostingRegressor
from run_analysis import FEATURES,DATA,OUT,get_file

def main():
 parts=[]
 for year,suffix in [(2013,'H'),(2015,'I')]:
  d=get_file((year,suffix,'DEMO'))[2][['SEQN','RIDAGEYR','RIAGENDR']]
  b=get_file((year,suffix,'BMX'))[2][['SEQN','BMXBMI','BMXWAIST']]
  p=get_file((year,suffix,'BPX'))[2]
  p['systolic']=p[['BPXSY1','BPXSY2','BPXSY3']].replace(0,np.nan).mean(axis=1)
  p['diastolic']=p[['BPXDI1','BPXDI2','BPXDI3']].mean(axis=1)
  t=get_file((year,suffix,'TRIGLY'))[2][['SEQN','LBXTR','WTSAF2YR']]
  d=d.merge(b,on='SEQN').merge(p[['SEQN','systolic','diastolic']],on='SEQN').merge(t,on='SEQN')
  parts.append(d[(d.RIDAGEYR>=20)&(d.LBXTR>0)&(d.WTSAF2YR>0)])
 d=pd.concat(parts);X=d[FEATURES];w=d.WTSAF2YR/d.WTSAF2YR.mean()
 pipe=Pipeline([('impute',SimpleImputer(strategy='median',add_indicator=True)),('scale',StandardScaler()),('model',HistGradientBoostingRegressor(max_iter=150,max_leaf_nodes=7,l2_regularization=10,learning_rate=.05,random_state=42))])
 pipe.fit(X,d.LBXTR,model__sample_weight=w)
 rng=np.random.default_rng(42);n=1000000;batch=10000;checksum=0.;start=time.perf_counter()
 for j in range(0,n,batch):
  z=X.iloc[rng.integers(0,len(X),size=min(batch,n-j))];pred=pipe.predict(z)
  assert np.isfinite(pred).all();checksum+=float(pred.sum())
 elapsed=time.perf_counter()-start
 r=dict(rows=n,batch_size=batch,seconds=elapsed,rows_per_second=n/elapsed,process_peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,checksum=checksum,source='Resampled training predictor rows; not new participants',timing_scope='Includes batch resampling and prediction; excludes model training and download',python=platform.python_version(),threads=2)
 (OUT/'scalability.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
if __name__=='__main__':main()
