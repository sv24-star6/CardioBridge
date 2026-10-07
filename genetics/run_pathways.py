"""Reactome annotations from UniProt cross-references, not enrichment testing."""
from pathlib import Path
import concurrent.futures,hashlib,json
import pandas as pd
from remote_tabix import fetch
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data';OUT=ROOT/'outputs'
def main():
 d=pd.read_csv(OUT/'candidate_instruments.csv');rows=[];prov=[]
 def get(r):
  path=DATA/f'uniprot_{r.uniprot}.txt';url=f'https://rest.uniprot.org/uniprotkb/{r.uniprot}.txt'
  if not path.exists():path.write_bytes(fetch(url))
  b=path.read_bytes();out=[]
  for line in b.decode().splitlines():
   if line.startswith('DR   Reactome;'):
    parts=[p.strip().rstrip('.') for p in line.split(';')]
    out.append(dict(gene=r.gene,uniprot=r.uniprot,reactome_id=parts[1],pathway=parts[2],MR_status=r.status))
  return out,dict(gene=r.gene,url=url,sha256=hashlib.sha256(b).hexdigest(),pathway_count=len(out))
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
  for a,p in ex.map(get,[r for _,r in d.iterrows()]):rows+=a;prov.append(p)
 pd.DataFrame(rows).to_csv(OUT/'reactome_annotations.csv',index=False);(OUT/'pathway_sources.json').write_text(json.dumps(prov,indent=2));print(pd.DataFrame(rows).to_string(index=False),flush=True)
if __name__=='__main__':main()
