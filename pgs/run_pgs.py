"""CardioBridge PGS module: calculate a published CAD score in a user-supplied genotype dosage table.

Default score: PGS012581 (PRS169_CAD), PGS Catalog.
This module does NOT ship individual-level genotypes and does NOT claim independent
clinical validation. It downloads the public scoring file and calculates scores only
when a compatible target genotype table is supplied.

Target format (TSV/TSV.GZ):
sample_id  rsID  effect_allele  other_allele  dosage
where dosage is the number (0..2) of copies of effect_allele.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,urllib.request
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs';CACHE=ROOT/'data';OUT.mkdir(parents=True,exist_ok=True);CACHE.mkdir(exist_ok=True)
PGS_ID='PGS012581'
URL=f'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{PGS_ID}/ScoringFiles/{PGS_ID}.txt.gz'
def download():
 p=CACHE/f'{PGS_ID}.txt.gz'
 if not p.exists():
  with urllib.request.urlopen(URL,timeout=120) as r:p.write_bytes(r.read())
 return p
def read_score(path):
 with gzip.open(path,'rt') as f:
  d=pd.read_csv(f,sep='\t',comment='#',low_memory=False)
 need={'effect_allele','effect_weight'}
 if not need.issubset(d):raise ValueError(f'Missing score columns: {need-set(d.columns)}')
 idcol='rsID' if 'rsID' in d.columns else None
 if idcol is None:raise ValueError('This implementation requires rsID in the scoring file.')
 d=d[d[idcol].notna()].copy();d.effect_weight=pd.to_numeric(d.effect_weight,errors='coerce');d=d[d.effect_weight.notna()]
 return d.rename(columns={idcol:'rsID'})
def score(genotypes,weights):
 g=pd.read_csv(genotypes,sep='\t',compression='infer')
 req={'sample_id','rsID','effect_allele','dosage'}
 if not req.issubset(g):raise ValueError(f'Missing genotype columns: {req-set(g.columns)}')
 g.effect_allele=g.effect_allele.str.upper();w=weights[['rsID','effect_allele','effect_weight']].copy();w.effect_allele=w.effect_allele.str.upper()
 z=g.merge(w,on=['rsID','effect_allele'],how='inner',validate='many_to_one');z.dosage=pd.to_numeric(z.dosage,errors='coerce');z=z[z.dosage.between(0,2)]
 z['contribution']=z.dosage*z.effect_weight
 n_total=w.rsid.nunique() if 'rsid' in w else w.rsID.nunique()
 result=z.groupby('sample_id').agg(raw_PGS=('contribution','sum'),variants_scored=('rsID','nunique')).reset_index();result['score_variant_fraction']=result.variants_scored/n_total
 return result,z
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--genotypes',required=False);a=ap.parse_args();p=download();w=read_score(p)
 manifest={'pgs_id':PGS_ID,'trait':'coronary artery disease','score_name':'PRS169_CAD','catalog_variant_count':169,'development_method':'Genome-wide significant SNPs; R2 < 0.01','effect_weight_type':'beta','scoring_file_url':URL,'scoring_file_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'variants_read':len(w),'validation_status':'No independent individual-level validation performed in CardioBridge.'}
 (OUT/'pgs_manifest.json').write_text(json.dumps(manifest,indent=2))
 if not a.genotypes:
  print(json.dumps(manifest,indent=2));print('No --genotypes supplied: metadata/weight QC only; no participant scores calculated.');return
 r,z=score(a.genotypes,w);r.to_csv(OUT/'participant_scores.csv',index=False);z.to_csv(OUT/'matched_score_variants.csv.gz',index=False,compression='gzip');print(r.describe(include='all').to_string())
if __name__=='__main__':main()
