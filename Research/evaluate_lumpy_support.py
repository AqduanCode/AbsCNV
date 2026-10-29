# Recovered from executed project turn 01a12fb2-96cc-7193-b11e-1af00cba3069.
# Analysis unchanged; shell quoting decoded and an optional input-root argument added.
import csv,gzip,re,json,sys
from pathlib import Path
if len(sys.argv)!=2: raise SystemExit("Usage: evaluate_lumpy_support.py PUBLIC_WGS_DATASET")
base=Path(sys.argv[1])
vcfs={'FASR':base/'LUMPY_SV'/'GSM3336915_FASR-lumpy-sv.vcf.gz','MCF7':base/'LUMPY_SV'/'GSM3336911_MCF7-lumpy-sv.vcf.gz','TAMR':base/'LUMPY_SV'/'GSM3336913_TAMR-lumpy-sv.vcf.gz'}
abspath=base/'output'/'CNVresult_filter.csv'

def parse(mode):
    out={}; filter_counts={}
    for sample,path in vcfs.items():
        rows=[]; fc={}
        with gzip.open(path,'rt',encoding='utf-8') as h:
            for line in h:
                if line.startswith('#'): continue
                f=line.rstrip().split('\t'); filt=f[6]; fc[filt]=fc.get(filt,0)+1
                info=dict(x.split('=',1) if '=' in x else (x,True) for x in f[7].split(';'))
                typ=str(info.get('SVTYPE','')); chrom=f[0]
                if typ not in {'DEL','DUP'} or not re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])',chrom): continue
                if mode=='strict' and filt!='PASS': continue
                if mode=='nofail' and filt not in {'.','PASS'}: continue
                pos=int(f[1]); end=int(info['END']); su=int(info.get('SU',0)); length=abs(int(info.get('SVLEN',end-pos)))
                if su<5 or length<10000: continue
                rows.append({'sample':sample,'id':f[2],'chrom':chrom,'start':min(pos,end),'end':max(pos,end),'type':typ,'length':length,'SU':su,'FILTER':filt})
        out[sample]=rows; filter_counts[sample]=fc
    return out,filter_counts

def ro(a,b):
    ov=max(0,min(a['end'],b['end'])-max(a['start'],b['start']))
    return 0 if ov==0 else min(ov/(a['end']-a['start']),ov/(b['end']-b['start']))
def same(a,b):
    return a['chrom']==b['chrom'] and a['type']==b['type'] and (ro(a,b)>=.5 or (abs(a['start']-b['start'])<=1000 and abs(a['end']-b['end'])<=1000))
def overlap(c,e): return c['chrom']==e['chrom'] and c['start']<e['end'] and c['end']>e['start']

events={}
with open(abspath,encoding='utf-8-sig',newline='') as h:
    for i,row in enumerate(csv.DictReader(h),1):
        a=row['pair_sampleA_name']; b=row['pair_sampleB_name']; ca=int(float(row['copy_number_A'])); cb=int(float(row['copy_number_B']))
        if ca==cb: continue
        events.setdefault((a,b),[]).append({'id':i,'chrom':row['chr'],'start':int(row['start_pos']),'end':int(row['end_pos']),'ca':ca,'cb':cb})

def counttype(rows): return {t:sum(x['type']==t for x in rows) for t in ('DEL','DUP')}
for mode in ('nofail',):
    calls,fc=parse(mode)
    result={'mode':mode,'filter_counts':fc,'sample_counts':{s:counttype(v) for s,v in calls.items()},'pairs':[]}
    for (a,b),evs in events.items():
        ua=[v for v in calls[a] if not any(same(v,x) for x in calls[b])]
        ub=[v for v in calls[b] if not any(same(v,x) for x in calls[a])]
        matched=set(); matched_by_type={'DEL':set(),'DUP':set()}; event_any=0; event_del=0; event_dup=0
        for e in evs:
            if e['ca']>e['cb']:
                expected=[(x,a) for x in ua if x['type']=='DUP']+[(x,b) for x in ub if x['type']=='DEL']
            else:
                expected=[(x,a) for x in ua if x['type']=='DEL']+[(x,b) for x in ub if x['type']=='DUP']
            hits=[(x,s) for x,s in expected if overlap(x,e)]
            if hits: event_any+=1
            if any(x['type']=='DEL' for x,s in hits): event_del+=1
            if any(x['type']=='DUP' for x,s in hits): event_dup+=1
            for x,s in hits:
                key=(s,x['id']); matched.add(key); matched_by_type[x['type']].add(key)
        result['pairs'].append({'pair':f'{a}-{b}','abscnv_events':len(evs),'a':a,'b':b,'a_unique':counttype(ua),'b_unique':counttype(ub),'pair_unique':counttype(ua+ub),'matched_lumpy_calls':{t:len(matched_by_type[t]) for t in ('DEL','DUP')},'matched_lumpy_calls_total':len(matched),'abscnv_events_any':event_any,'abscnv_events_with_DEL_support':event_del,'abscnv_events_with_DUP_support':event_dup})
    print(json.dumps(result,ensure_ascii=False))
