"""Clinical-utility evaluation for CardioBridge.

Requires a linked validation dataset with observed binary CAD outcome, clinical
predictors and a calculated PGS. This script evaluates incremental predictive
value; it does NOT by itself prove improved patient outcomes.

Expected TSV columns:
CAD, age, sex, clinical_score, PGS
clinical_score and PGS should be numeric. CAD must be 0/1.
"""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,brier_score_loss
from sklearn.model_selection import StratifiedKFold,cross_val_predict
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'outputs';OUT.mkdir(exist_ok=True)
def calibration(y,p):
 eps=1e-8;lp=np.log(np.clip(p,eps,1-eps)/(1-np.clip(p,eps,1-eps)))
 X=np.c_[np.ones(len(lp)),lp]
 beta=np.linalg.lstsq(X,np.log((y+.5)/(1-y+.5)),rcond=None)[0]
 return {'calibration_intercept_approx':float(beta[0]),'calibration_slope_approx':float(beta[1])}
def net_benefit(y,p,t):
 pred=p>=t;n=len(y);tp=np.sum(pred&(y==1));fp=np.sum(pred&(y==0))
 return (tp/n)-(fp/n)*(t/(1-t))
def nri_categories(y,p0,p1,cuts=(.05,.10,.20)):
 c0=np.digitize(p0,cuts);c1=np.digitize(p1,cuts);case=y==1;ctrl=~case
 case_nri=np.mean(c1[case]>c0[case])-np.mean(c1[case]<c0[case])
 ctrl_nri=np.mean(c1[ctrl]<c0[ctrl])-np.mean(c1[ctrl]>c0[ctrl])
 return float(case_nri+ctrl_nri)
def fit_cv(X,y):
 cv=StratifiedKFold(5,shuffle=True,random_state=42)
 m=LogisticRegression(max_iter=5000)
 return cross_val_predict(m,X,y,cv=cv,method='predict_proba')[:,1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--data',required=True);a=ap.parse_args()
 d=pd.read_csv(a.data,sep='\t',compression='infer').dropna(subset=['CAD','age','sex','clinical_score','PGS'])
 y=d.CAD.astype(int).to_numpy()
 if set(np.unique(y))-set([0,1]):raise ValueError('CAD must be binary 0/1')
 models={'Clinical':d[['age','sex','clinical_score']],'PGS':d[['PGS']],'Clinical + PGS':d[['age','sex','clinical_score','PGS']]}
 preds={k:fit_cv(x,y) for k,x in models.items()};rows=[]
 for k,p in preds.items():
  row={'model':k,'n':len(y),'events':int(y.sum()),'AUC':float(roc_auc_score(y,p)),'Brier':float(brier_score_loss(y,p))}
  row.update(calibration(y,p));rows.append(row)
 pd.DataFrame(rows).to_csv(OUT/'model_utility_metrics.csv',index=False)
 thresholds=np.arange(.01,.301,.01);nb=[]
 for t in thresholds:
  prevalence=y.mean();nb.append({'threshold':t,'model':'Treat all','net_benefit':float(prevalence-(1-prevalence)*t/(1-t))});nb.append({'threshold':t,'model':'Treat none','net_benefit':0.0})
  for k,p in preds.items():nb.append({'threshold':t,'model':k,'net_benefit':float(net_benefit(y,p,t))})
 pd.DataFrame(nb).to_csv(OUT/'decision_curve.csv',index=False)
 nri=nri_categories(y,preds['Clinical'],preds['Clinical + PGS'])
 summary={'n':len(y),'events':int(y.sum()),'clinical_auc':float(roc_auc_score(y,preds['Clinical'])),'clinical_plus_pgs_auc':float(roc_auc_score(y,preds['Clinical + PGS'])),'delta_auc':float(roc_auc_score(y,preds['Clinical + PGS'])-roc_auc_score(y,preds['Clinical'])),'category_NRI_5_10_20pct':nri,'interpretation':'Predictive/decision-curve evaluation only. Clinical utility requires context-specific evidence that use improves decisions or outcomes.'}
 (OUT/'utility_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
