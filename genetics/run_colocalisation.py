"""Wakefield ABF colocalisation under one causal variant per trait per region."""
from pathlib import Path
import concurrent.futures,hashlib,json
import numpy as np,pandas as pd
from scipy.special import logsumexp
from remote_tabix import fetch
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data';OUT=ROOT/'outputs'
FILES={'ANGPTL1':(1,'295397982100'),'APOE':(19,'295413323771'),'APOM':(6,'295410206292'),'APOL1':(22,'295410729967'),'IL6R':(1,'295731216403'),'CRP':(1,'295550192434')}
SHARE='u3flbp13zjydegrxjb2uepagp1vb6bj2'
def abf(beta,se,prior_sd):
 v=np.asarray(se)**2;r=prior_sd**2/(prior_sd**2+v);return .5*(np.log1p(-r)+r*(np.asarray(beta)/se)**2)
def posterior(l1,l2,p1=1e-4,p2=1e-4,p12=1e-5):
 n=len(l1);left=np.r_[-np.inf,np.logaddexp.accumulate(l2)[:-1]];right=np.r_[np.logaddexp.accumulate(l2[::-1])[::-1][1:],-np.inf];off=logsumexp(l1+np.logaddexp(left,right));logs=np.array([0,np.log(p1)+logsumexp(l1),np.log(p2)+logsumexp(l2),np.log(p1)+np.log(p2)+off,np.log(p12)+logsumexp(l1+l2)]);return np.exp(logs-logsumexp(logs))
def main():
 instruments=pd.read_csv(OUT/'candidate_instruments.csv');results=[];provenance=[];qc=[]
 def run(gene):
  ch,fid=FILES[gene];path=DATA/f'{gene}_chr{ch}.tsv.gz';url=f'https://app.box.com/index.php?rm=box_download_shared_file&shared_name={SHARE}&file_id=f_{fid}'
  if not path.exists():
   b=fetch(url)
   if b[:2]!=b'\x1f\x8b':raise ValueError('Not gzip data: '+gene)
   path.write_bytes(b)
  r=instruments[instruments.gene==gene].iloc[0];pos=int(r.pos37);chunks=[]
  for chunk in pd.read_csv(path,sep='\t',chunksize=100000):
   sub=chunk[chunk.position.between(pos-500000,pos+500000)]
   if len(sub):chunks.append(sub)
  p=pd.concat(chunks,ignore_index=True);lead=p[p.position==pos];lead=lead[lead.apply(lambda x:set([str(x.Allele1).upper(),str(x.Allele2).upper()])==set([r.EA,r.OA]),axis=1)]
  assert len(lead)==1;assert np.isclose(abs(lead.iloc[0].Effect),abs(r.bx),atol=.00015);assert np.isclose(lead.iloc[0].StdErr,r.sex,atol=.00015)
  p['a1']=p.Allele1.str.upper();p['a2']=p.Allele2.str.upper();p=p[p.a1.isin(list('ACGT'))&p.a2.isin(list('ACGT'))&(p.a1!=p.a2)].copy();p['key']=p.position.astype(str)+'_'+p.apply(lambda x:'_'.join(sorted([x.a1,x.a2])),axis=1)
  o=pd.read_csv(DATA/f'cad_region_{gene}.tsv',sep='\t');split=o.markername.str.split(r'[:_]',regex=True);o['pos37']=pd.to_numeric(split.str[1]);o['a1']=split.str[2].str.upper();o['a2']=split.str[3].str.upper();o=o[o.pos37.between(pos-500000,pos+500000)&o.a1.isin(list('ACGT'))&o.a2.isin(list('ACGT'))&(o.a1!=o.a2)].copy();o['key']=o.pos37.astype(int).astype(str)+'_'+o.apply(lambda x:'_'.join(sorted([x.a1,x.a2])),axis=1)
  pd_dups=p.key.duplicated(keep=False);od_dups=o.key.duplicated(keep=False);n_p=len(p);n_o=len(o);p=p[~pd_dups];o=o[~od_dups];z=p.merge(o[['key','beta','standard_error','p_value','rsid']],on='key',validate='one_to_one');z=z[np.isfinite(z[['Effect','StdErr','beta','standard_error']]).all(axis=1)&(z.StdErr>0)&(z.standard_error>0)].copy();assert len(z)>100 and z.position.eq(pos).any()
  l1=abf(z.Effect.to_numpy(),z.StdErr.to_numpy(),.15);l2=abf(z.beta.to_numpy(),z.standard_error.to_numpy(),.2);rows=[]
  for prior in [1e-6,1e-5,1e-4]:
   p1_eff=min(1e-4,1/(len(z)+1));p2_eff=min(1e-4,1/(len(z)+1));p12_eff=min(prior,1/(len(z)+1));pp=posterior(l1,l2,p1=p1_eff,p2=p2_eff,p12=p12_eff);rows.append(dict(gene=gene,n_shared=len(z),p1=p1_eff,p2=p2_eff,p12=prior,p12_effective=p12_eff,**dict(zip(['PP_H0','PP_H1','PP_H2','PP_H3','PP_H4'],pp))))
  z['protein_logABF']=l1;z['CAD_logABF']=l2;z['SNP_PP_given_H4']=np.exp(l1+l2-logsumexp(l1+l2));z.to_csv(OUT/f'coloc_variants_{gene}.csv.gz',index=False,compression='gzip')
  quality=dict(gene=gene,chr37=ch,start37=pos-500000,end37=pos+500000,protein_SNPs=n_p,CAD_SNPs=n_o,shared_SNPs=len(z),protein_coverage=len(z)/n_p,CAD_coverage=len(z)/n_o,protein_duplicate_rows=int(pd_dups.sum()),CAD_duplicate_rows=int(od_dups.sum()),lead_beta_agrees=True,lead_se_agrees=True);source=dict(gene=gene,url=url,filename=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest());return rows,quality,source
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
  for rows,q,s in ex.map(run,FILES):results+=rows;qc.append(q);provenance.append(s)
 pd.DataFrame(results).to_csv(OUT/'colocalisation_sensitivity.csv',index=False);pd.DataFrame(qc).to_csv(OUT/'colocalisation_qc.csv',index=False);(OUT/'colocalisation_sources.json').write_text(json.dumps(provenance,indent=2));print(pd.DataFrame(results).to_string(index=False),flush=True)
if __name__=='__main__':main()
