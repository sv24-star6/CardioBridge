"""CardioBridge v0.1: temporally held-out NHANES triglyceride benchmark."""
from pathlib import Path
import concurrent.futures, hashlib, io, json, time, urllib.request, platform
import numpy as np
import pandas as pd
import sklearn
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor, HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'; OUT=ROOT/'outputs'
DATA.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
CYCLES=[(2013,'H'),(2015,'I'),(2017,'J')]
FEATURES=['RIDAGEYR','RIAGENDR','BMXBMI','BMXWAIST','systolic','diastolic']
def get_file(task):
 year,suffix,component=task
 name=f'{component}_{suffix}.xpt'
 url=f'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{name}'
 path=DATA/name
 if not path.exists():
  with urllib.request.urlopen(url,timeout=90) as r: content=r.read()
  if not content.startswith(b'HEADER RECORD'): raise ValueError(f'Not XPT: {url}')
  path.write_bytes(content)
 content=path.read_bytes()
 return (year,component,pd.read_sas(io.BytesIO(content),format='xport'),dict(file=name,url=url,sha256=hashlib.sha256(content).hexdigest(),bytes=len(content)))
def metrics(y,p,w):
 return dict(MAE_mg_dL=float(mean_absolute_error(y,p,sample_weight=w)),RMSE_mg_dL=float(np.sqrt(mean_squared_error(y,p,sample_weight=w))),R2=float(r2_score(y,p,sample_weight=w)))
def main():
 tasks=[(y,s,c) for y,s in CYCLES for c in ['DEMO','BMX','BPX','TRIGLY']]
 with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex: downloaded=list(ex.map(get_file,tasks))
 manifest=[v[3] for v in downloaded]
 frames={(y,c):df for y,c,df,_ in downloaded}
 cohorts=[]; flow=[]
 for year,_ in CYCLES:
  d=frames[year,'DEMO'][['SEQN','RIDAGEYR','RIAGENDR','SDMVPSU','SDMVSTRA']]
  b=frames[year,'BMX'][['SEQN','BMXBMI','BMXWAIST']]
  p=frames[year,'BPX'].copy()
  p['systolic']=p[['BPXSY1','BPXSY2','BPXSY3']].replace(0,np.nan).mean(axis=1)
  p['diastolic']=p[['BPXDI1','BPXDI2','BPXDI3']].mean(axis=1)
  t=frames[year,'TRIGLY'][['SEQN','LBXTR','WTSAF2YR']]
  d=d.merge(b,on='SEQN',how='left',validate='one_to_one').merge(p[['SEQN','systolic','diastolic']],on='SEQN',how='left',validate='one_to_one').merge(t,on='SEQN',how='left',validate='one_to_one')
  adult=d.RIDAGEYR>=20
  eligible=adult & d.LBXTR.notna() & (d.LBXTR>0) & (d.WTSAF2YR>0)
  flow.append(dict(cycle=year,total=len(d),adults=int(adult.sum()),eligible=int(eligible.sum())))
  d=d.loc[eligible].copy(); d['cycle']=year; cohorts.append(d)
 df=pd.concat(cohorts,ignore_index=True)
 train=df[df.cycle<2017].copy(); test=df[df.cycle==2017].copy()
 assert not set(train.SEQN)&set(test.SEQN)
 X=train[FEATURES]; Z=test[FEATURES]; y=train.LBXTR; yt=test.LBXTR
 w=train.WTSAF2YR.to_numpy(); w=w/w.mean(); wt=test.WTSAF2YR.to_numpy()
 models={
  'Elastic net':ElasticNet(alpha=1.0,l1_ratio=.1,max_iter=10000,random_state=42),
  'Random forest':RandomForestRegressor(n_estimators=250,min_samples_leaf=20,max_features=1.0,random_state=42,n_jobs=2),
  'Gradient boosting':HistGradientBoostingRegressor(max_iter=150,max_leaf_nodes=7,l2_regularization=10,learning_rate=.05,random_state=42)}
 results=[]; preds={}
 baseline=np.full(len(test),np.average(y,weights=w)); preds['Mean baseline']=baseline
 results.append(dict(model='Mean baseline',**metrics(yt,baseline,wt),train_seconds=0))
 for name,model in models.items():
  pipe=Pipeline([('impute',SimpleImputer(strategy='median',add_indicator=True)),('scale',StandardScaler()),('model',model)])
  start=time.perf_counter(); pipe.fit(X,y,model__sample_weight=w); elapsed=time.perf_counter()-start
  pred=pipe.predict(Z); preds[name]=pred
  results.append(dict(model=name,**metrics(yt,pred,wt),train_seconds=elapsed))
  print(name,results[-1],flush=True)
 table=pd.DataFrame(results); table.to_csv(OUT/'model_comparison.csv',index=False)
 pd.DataFrame(flow).to_csv(OUT/'cohort_flow.csv',index=False)
 missing=df.groupby('cycle')[FEATURES].agg(lambda x:x.isna().mean()); missing.to_csv(OUT/'missingness.csv')
 prediction=test[['SEQN','cycle','LBXTR','WTSAF2YR']].copy()
 for k,v in preds.items(): prediction[k]=v
 prediction.to_csv(OUT/'test_predictions.csv',index=False)
 subgroup=[]
 for label,mask in [('Female',test.RIAGENDR==2),('Male',test.RIAGENDR==1),('Age 20–49',test.RIDAGEYR<50),('Age 50+',test.RIDAGEYR>=50)]:
  for name,pred in preds.items(): subgroup.append(dict(group=label,n=int(mask.sum()),model=name,**metrics(yt[mask],pred[mask],wt[mask])))
 pd.DataFrame(subgroup).to_csv(OUT/'subgroup_metrics.csv',index=False)
 fig,axes=plt.subplots(1,2,figsize=(11,4.6))
 axes[0].barh(table.model,table.RMSE_mg_dL); axes[0].set_xlabel('Weighted RMSE (mg/dL); lower is better'); axes[0].set_title('2017–2018 held-out evaluation')
 axes[1].scatter(yt,preds['Gradient boosting'],s=9,alpha=.3); lim=max(yt.max(),preds['Gradient boosting'].max()); axes[1].plot([0,lim],[0,lim],'--'); axes[1].set_xlabel('Measured triglycerides (mg/dL)'); axes[1].set_ylabel('Predicted triglycerides (mg/dL)'); axes[1].set_title('Gradient boosting: all eligible test records')
 fig.tight_layout(); fig.savefig(OUT/'benchmark.png',dpi=180); plt.close(fig)
 (OUT/'provenance.json').write_text(json.dumps(dict(files=manifest,python=platform.python_version(),pandas=pd.__version__,sklearn=sklearn.__version__,features=FEATURES,train_n=len(train),test_n=len(test),seed=42),indent=2))
 print('COHORT',flow,flush=True); print(table.to_string(index=False),flush=True)
if __name__=='__main__': main()
