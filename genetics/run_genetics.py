"""Exploratory cis-pQTL Wald-ratio MR. No clinical/causal validation claim."""
from pathlib import Path
import concurrent.futures,json,hashlib,warnings
import numpy as np,pandas as pd
from scipy.stats import norm
from remote_tabix import RemoteTabix,fetch
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data';OUT=ROOT/'outputs';OUT.mkdir(exist_ok=True);DATA.mkdir(exist_ok=True)
URL='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/harmonised/GCST90132314.h.tsv.gz'
SUN='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/41586_2018_175_MOESM4_ESM.xlsx'
PANEL=['ANGPTL1','APOM','APOA5','APOF','APOE','APOL1','CRP','IL6R']
def prepare():
 if not (DATA/'sun2018.xlsx').exists():(DATA/'sun2018.xlsx').write_bytes(fetch(SUN))
 if not (DATA/'cad.tbi').exists():(DATA/'cad.tbi').write_bytes(fetch(URL+'.tbi'))
 with warnings.catch_warnings():
  warnings.simplefilter('ignore');d=pd.read_excel(DATA/'sun2018.xlsx',sheet_name='ST4 - pQTL summary',header=None).iloc[6:]
 d=d.rename(columns={1:'assay',2:'protein_label',3:'protein_name',4:'uniprot',5:'rsid',6:'chr37',7:'pos37',10:'EA',11:'OA',12:'EAF',14:'location',22:'bx',23:'sex',24:'px',28:'uncorrelated_PAV',29:'PAV_adjusted',30:'cis_eQTL'});d['gene']=d.assay.astype(str).str.split('.').str[0];d=d[d.gene.isin(PANEL)&(d.location=='cis')].copy()
 for c in ['bx','sex','pos37','chr37','EAF']:d[c]=pd.to_numeric(d[c])
 d['F']=(d.bx/d.sex)**2;d=d.sort_values('F',ascending=False).drop_duplicates('gene').sort_values('gene');cols=['gene','assay','protein_label','protein_name','uniprot','rsid','chr37','pos37','EA','OA','EAF','bx','sex','px','F','uncorrelated_PAV','PAV_adjusted','cis_eQTL'];d=d[cols].copy();d['status']='eligible'
 for idx,r in d.iterrows():
  if len(r.EA)!=1 or len(r.OA)!=1:d.loc[idx,'status']='excluded_indel_requires_normalisation'
  elif set([r.EA,r.OA]) in [set('AT'),set('CG')]:d.loc[idx,'status']='excluded_palindromic_conservative'
  elif r.F<10:d.loc[idx,'status']='excluded_weak_instrument'
 d.to_csv(OUT/'candidate_instruments.csv',index=False);return d
def main():
 d=prepare();t=RemoteTabix(URL,DATA/'cad.tbi',DATA/'ranges');header=t.header().split('\t')
 def lookup(r):
  ch=int(r.chr37);pos=int(r.pos37);cache=DATA/f'cad_region_{r.gene}.tsv'
  if not cache.exists():cache.write_text('\t'.join(header)+'\n'+'\n'.join(t.query(str(ch),max(0,pos-2000000),pos+2000000))+'\n')
  region=pd.read_csv(cache,sep='\t',dtype={'chromosome':str});original=region.markername.astype(str).str.split(r'[:_]',regex=True);return r,region[original.str[0].eq(str(ch))&original.str[1].eq(str(pos))].copy()
 rows=[];records=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
  for r,matches in ex.map(lookup,[r for _,r in d[d.status=='eligible'].iterrows()]):
   accepted=[]
   for _,o in matches.iterrows():
    ea,oa=str(o.effect_allele).upper(),str(o.other_allele).upper();comp=str.maketrans('ACGT','TGCA');pairs=[(r.EA,r.OA,1,'direct'),(r.OA,r.EA,-1,'swapped'),(r.EA.translate(comp),r.OA.translate(comp),1,'complement'),(r.OA.translate(comp),r.EA.translate(comp),-1,'complement_swapped')];hit=[(s,label) for a,b,s,label in pairs if (ea,oa)==(a,b)]
    if len(hit)!=1:continue
    sign,label=hit[0];eaf=float(o.effect_allele_frequency);aligned_eaf=eaf if sign==1 else 1-eaf
    if abs(aligned_eaf-r.EAF)<=.15:accepted.append((o,sign,label,aligned_eaf))
   if len(accepted)!=1:d.loc[d.gene==r.gene,'status']='excluded_no_unique_harmonised_match';continue
   o,sign,label,oeaf=accepted[0];assert str(o.rsid)==r.rsid;by=float(o.beta)*sign;sy=float(o.standard_error);theta=by/r.bx;se1=sy/abs(r.bx);se=np.sqrt(sy**2/r.bx**2+by**2*r.sex**2/r.bx**4)
   rows.append(dict(gene=r.gene,assay=r.assay,rsid=r.rsid,uniprot=r.uniprot,chr37=int(r.chr37),pos37=int(r.pos37),pos38=int(o.base_pair_location),EA=r.EA,OA=r.OA,protein_beta=r.bx,protein_se=r.sex,F=r.F,CAD_beta_aligned=by,CAD_se=sy,CAD_p=o.p_value,CAD_EAF_aligned=oeaf,protein_EAF=r.EAF,harmonisation=label,log_OR=theta,SE_delta=se,OR=np.exp(theta),CI_low=np.exp(theta-1.96*se),CI_high=np.exp(theta+1.96*se),p=2*norm.sf(abs(theta/se)),SE_first_order=se1,p_first_order=2*norm.sf(abs(theta/se1)),CAD_n=o.n,CAD_cases=o.cases));rec=o.to_dict();rec['gene']=r.gene;records.append(rec);d.loc[d.gene==r.gene,'status']='analysed'
 result=pd.DataFrame(rows).sort_values('p');result['p_bonferroni_8']=np.minimum(result.p*len(PANEL),1);result.to_csv(OUT/'mr_results.csv',index=False);pd.DataFrame(records).to_csv(OUT/'matched_CAD_records.csv',index=False);d.to_csv(OUT/'candidate_instruments.csv',index=False);(OUT/'sources.json').write_text(json.dumps({'protein_source':SUN,'protein_doi':'10.1038/s41586-018-0175-2','outcome_source':URL,'outcome_doi':'10.1038/s41588-022-01233-6','panel':PANEL,'genome_matching':'Original GRCh37 markername matches supplement coordinates; harmonised CAD position is GRCh38','sun_sha256':hashlib.sha256((DATA/'sun2018.xlsx').read_bytes()).hexdigest(),'independence':'Exposure/outcome covariance assumed zero; participant-level overlap not verified'},indent=2));print(result.to_string(index=False),flush=True)
if __name__=='__main__':main()
