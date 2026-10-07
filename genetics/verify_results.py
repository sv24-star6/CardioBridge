"""Numerical and source consistency checks for genetic analyses."""
from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from scipy.stats import norm
from run_colocalisation import posterior,abf
OUT=Path(__file__).resolve().parent/'outputs'
def main():
 checks={};b1=np.array([2.,5.,11.]);b2=np.array([7.,3.,13.]);p1,p2,p12=1e-4,1e-4,1e-5
 weights=np.array([1,p1*sum(b1),p2*sum(b2),p1*p2*sum(b1[i]*b2[j] for i,j in itertools.product(range(3),repeat=2) if i!=j),p12*sum(b1*b2)])
 assert np.allclose(posterior(np.log(b1),np.log(b2)),weights/weights.sum(),atol=1e-14);checks['coloc_matches_independent_hypothesis_enumeration']=True
 beta=np.array([0.,.1,.5]);se=np.array([.03,.05,.1]);prior=.15;direct=norm.pdf(beta,0,np.sqrt(se**2+prior**2))/norm.pdf(beta,0,se)
 assert np.allclose(np.exp(abf(beta,se,prior)),direct);checks['ABF_matches_normal_density_ratio']=True
 m=pd.read_csv(OUT/'mr_results.csv');c=pd.read_csv(OUT/'colocalisation_sensitivity.csv');q=pd.read_csv(OUT/'colocalisation_qc.csv')
 assert len(m)==6 and m.gene.is_unique and m.rsid.is_unique;assert np.allclose(m.OR,np.exp(m.CAD_beta_aligned/m.protein_beta));assert np.allclose(m.SE_delta,np.sqrt((m.CAD_se/m.protein_beta)**2+(m.CAD_beta_aligned*m.protein_se/m.protein_beta**2)**2));assert np.allclose(m.p_bonferroni_8,np.minimum(8*m.p,1));assert (m.F>10).all() and (abs(m.protein_EAF-m.CAD_EAF_aligned)<.15).all();checks['MR_ratios_uncertainty_multiplicity_and_instruments']=True
 pp=c[['PP_H0','PP_H1','PP_H2','PP_H3','PP_H4']];assert np.isfinite(pp).all().all() and np.allclose(pp.sum(axis=1),1) and (pp>=0).all().all();assert q.lead_beta_agrees.all() and q.lead_se_agrees.all() and (q.shared_SNPs>100).all();checks['coloc_probabilities_and_downloaded_sentinel_consistency']=True
 for gene in m.gene:
  z=pd.read_csv(OUT/f'coloc_variants_{gene}.csv.gz');assert z.key.is_unique and np.isclose(z.SNP_PP_given_H4.sum(),1)
 checks['variant_keys_unique_and_conditional_posteriors_normalised']=True;(OUT/'verification.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
if __name__=='__main__':main()
