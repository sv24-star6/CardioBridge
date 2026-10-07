"""Minimal read-only TBI/BGZF region retrieval with verified HTTP byte ranges."""
import gzip,struct,urllib.request,time
from pathlib import Path
def fetch(url,start=None,end=None):
 for attempt in range(3):
  try:
   req=urllib.request.Request(url,headers={} if start is None else {'Range':f'bytes={start}-{end}'})
   with urllib.request.urlopen(req,timeout=90) as r:
    if start is not None and (r.status!=206 or not r.headers.get('Content-Range','').startswith(f'bytes {start}-')): raise ValueError('Server did not honour byte range')
    return r.read()
  except Exception:
   if attempt==2: raise
   time.sleep(1)
def parse_index(content):
 b=gzip.decompress(content);off=4;assert b[:4]==b'TBI\x01'
 def take(fmt):
  nonlocal off
  v=struct.unpack_from('<'+fmt,b,off);off+=struct.calcsize('<'+fmt);return v
 n,fmt,colseq,colbeg,colend,meta,skip,l_nm=take('8i');names=b[off:off+l_nm].rstrip(b'\0').decode().split('\0');off+=l_nm;refs={}
 for name in names:
  nb,=take('i');bins={}
  for j in range(nb):
   bin_id,nc=take('Ii');bins[bin_id]=[take('QQ') for _ in range(nc)]
  ni,=take('i');linear=[take('Q')[0] for _ in range(ni)];refs[name]=(bins,linear)
 return refs,dict(format=fmt,colseq=colseq,colbeg=colbeg,colend=colend,meta=meta)
def region_bins(start,end):
 end-=1
 return [0]+list(range(1+(start>>26),2+(end>>26)))+list(range(9+(start>>23),10+(end>>23)))+list(range(73+(start>>20),74+(end>>20)))+list(range(585+(start>>17),586+(end>>17)))+list(range(4681+(start>>14),4682+(end>>14)))
def decompress_blocks(raw):
 blocks={};i=0
 while i+18<=len(raw):
  if raw[i:i+2]!=b'\x1f\x8b':break
  xlen=struct.unpack_from('<H',raw,i+10)[0];size=None;j=i+12
  while j<i+12+xlen:
   slen=struct.unpack_from('<H',raw,j+2)[0]
   if raw[j:j+2]==b'BC':size=struct.unpack_from('<H',raw,j+4)[0]+1
   j+=4+slen
  if size is None or i+size>len(raw):break
  blocks[i]=gzip.decompress(raw[i:i+size]);i+=size
 return blocks
class RemoteTabix:
 def __init__(self,url,index_path,cache):
  self.url=url;self.cache=Path(cache);self.cache.mkdir(exist_ok=True);self.refs,self.meta=parse_index(Path(index_path).read_bytes())
 def rawrange(self,start,end):
  path=self.cache/f'{start}_{end}.bin'
  if not path.exists():path.write_bytes(fetch(self.url,start,end))
  return path.read_bytes()
 def header(self):
  return b''.join(decompress_blocks(self.rawrange(0,131071)).values()).decode().splitlines()[0]
 def query(self,chrom,start,end):
  bins,linear=self.refs[str(chrom)];floor=linear[start>>14] if start>>14<len(linear) else 0;chunks=[]
  for b in region_bins(start,end):chunks.extend((a,c) for a,c in bins.get(b,[]) if c>floor)
  chunks.sort();merged=[]
  for a,b in chunks:
   if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
   else:merged.append((a,b))
  lines=set()
  for a,b in merged:
   lo=a>>16;hi=(b>>16)+65535;blocks=decompress_blocks(self.rawrange(lo,hi));pieces=[]
   for rel,data in blocks.items():
    addr=lo+rel
    if addr>b>>16:break
    left=(a&65535) if addr==a>>16 else 0;right=(b&65535) if addr==b>>16 else len(data);pieces.append(data[left:right])
   for line in b''.join(pieces).decode().splitlines():
    fields=line.split('\t')
    if len(fields)<3:continue
    try:
     ch=fields[self.meta['colseq']-1];pos=int(fields[self.meta['colbeg']-1])
     if ch==str(chrom) and start<pos<=end:lines.add(line)
    except ValueError:continue
  return sorted(lines)
