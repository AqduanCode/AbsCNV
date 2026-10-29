"""Direction-aware interval metrics used by the paired benchmark."""
from __future__ import annotations
import csv,math
from pathlib import Path

def read_table(path: Path):
    delim = "\t" if path.suffix.lower() in (".tsv", ".benchmark") else ","
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))

def write_table(path: Path, rows, fields=None):
    if not rows and not fields: return
    fields = fields or list(rows[0].keys())
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        w=csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)

def load_truth(path):
    out=[]
    for i,r in enumerate(read_table(path),1):
        d=float(r["delta_B_minus_A"])
        if d:
            out.append(dict(id=r.get("event_id") or f"truth_{i}",chrom=r["chrom"],start=int(r["start"]),end=int(r["end"]),direction="gain" if d>0 else "loss"))
    return out

def load_my(path):
    out=[]
    for i,r in enumerate(read_table(path),1):
        if r.get("copy_number_A","") != "" and r.get("copy_number_B","") != "":
            d=float(r["copy_number_B"])-float(r["copy_number_A"])
        else: d=float(r["depth_diff"])
        if d:
            out.append(dict(id=f"MyMethod_{i}",chrom=r["chr"],start=int(float(r["start_pos"])),end=int(float(r["end_pos"])),direction="gain" if d>0 else "loss"))
    return out

def load_normalized(path):
    return [dict(id=r["id"],chrom=r["chrom"],start=int(r["start"]),end=int(r["end"]),direction=r["direction"].lower()) for r in read_table(path)]

def validate(records, label):
    ids=set()
    for r in records:
        if r["id"] in ids: raise ValueError(f"{label}: duplicate id {r['id']}")
        ids.add(r["id"])
        if r["direction"] not in ("gain","loss"): raise ValueError(f"{label}: invalid direction {r['direction']}")
        if r["start"] >= r["end"]: raise ValueError(f"{label}: non-positive interval {r}")

def overlap(a,b): return max(0,min(a["end"],b["end"])-max(a["start"],b["start"]))

def metric(tp,np,nt):
    p=tp/np if np else math.nan; r=tp/nt if nt else math.nan
    return p,r,2*p*r/(p+r) if p+r else 0.0

def candidates(preds, truths, threshold, mode):
    ans=[]
    for p in preds:
        plen=p["end"]-p["start"]
        for t in truths:
            if p["chrom"]!=t["chrom"] or p["direction"]!=t["direction"]: continue
            ov=overlap(p,t); tlen=t["end"]-t["start"]
            if ov<=0: continue
            pf,tf=ov/plen,ov/tlen
            ok=(pf>=threshold and tf>=threshold) if mode=="reciprocal" else (pf>=threshold or tf>=threshold)
            if ok:
                iou=ov/(plen+tlen-ov)
                ans.append((min(pf,tf),iou,ov,pf,tf,p,t))
    return sorted(ans,key=lambda x:(x[0],x[1],x[2]),reverse=True)

def event_eval(tool,preds,truths,standard,threshold,mode,multiplicity):
    cand=candidates(preds,truths,threshold,mode); used_p=set(); used_t=set(); matches=[]
    for score,iou,ov,pf,tf,p,t in cand:
        if p["id"] in used_p: continue
        if multiplicity=="one_to_one" and t["id"] in used_t: continue
        used_p.add(p["id"]); used_t.add(t["id"])
        matches.append(dict(standard=standard,tool=tool,prediction_id=p["id"],truth_id=t["id"],chrom=p["chrom"],direction=p["direction"],prediction_start=p["start"],prediction_end=p["end"],truth_start=t["start"],truth_end=t["end"],overlap_bp=ov,prediction_overlap_fraction=pf,truth_overlap_fraction=tf,iou=iou))
    rows=[]
    for direction in ("all","gain","loss"):
        np=sum(direction=="all" or p["direction"]==direction for p in preds)
        nt=sum(direction=="all" or t["direction"]==direction for t in truths)
        matched_preds={m["prediction_id"] for m in matches if direction=="all" or m["direction"]==direction}
        matched_truth={m["truth_id"] for m in matches if direction=="all" or m["direction"]==direction}
        if multiplicity=="one_to_one":
            tp=len(matched_preds); precision,recall,f1=metric(tp,np,nt)
        else:
            precision=len(matched_preds)/np if np else math.nan
            recall=len(matched_truth)/nt if nt else math.nan
            f1=2*precision*recall/(precision+recall) if precision+recall else 0.0
            tp=len(matched_preds)
        rows.append(dict(standard=standard,tool=tool,direction=direction,matched_predictions=len(matched_preds),matched_truth_events=len(matched_truth),false_positive_predictions=np-len(matched_preds),false_negative_truth_events=nt-len(matched_truth),predicted_events=np,truth_events=nt,precision=precision,recall=recall,f1=f1))
    up=[dict(standard=standard,tool=tool,**p) for p in preds if p["id"] not in used_p]
    ut=[dict(standard=standard,tool=tool,**t) for t in truths if t["id"] not in used_t]
    return rows,matches,up,ut

def merge_by_key(records):
    groups={}
    for r in records: groups.setdefault((r["chrom"],r["direction"]),[]).append((r["start"],r["end"]))
    for k,vals in groups.items():
        merged=[]
        for s,e in sorted(vals):
            if not merged or s>merged[-1][1]: merged.append([s,e])
            else: merged[-1][1]=max(merged[-1][1],e)
        groups[k]=merged
    return groups

def sum_bp(groups,d=None): return sum(e-s for (_,direction),vals in groups.items() if d is None or d==direction for s,e in vals)

def intersect_bp(a,b,d=None):
    total=0
    for key,av in a.items():
        if d is not None and key[1]!=d: continue
        bv=b.get(key,[]); i=j=0
        while i<len(av) and j<len(bv):
            total+=max(0,min(av[i][1],bv[j][1])-max(av[i][0],bv[j][0]))
            if av[i][1]<=bv[j][1]: i+=1
            else: j+=1
    return total

def base_eval(tool,preds,truths,standard):
    pm,tm=merge_by_key(preds),merge_by_key(truths); rows=[]
    for label,d in (("all",None),("gain","gain"),("loss","loss")):
        pb,tb=sum_bp(pm,d),sum_bp(tm,d); mb=intersect_bp(pm,tm,d); p,r,f=metric(mb,pb,tb)
        rows.append(dict(standard=standard,tool=tool,direction=label,matched_bp=mb,predicted_cnv_bp=pb,truth_cnv_bp=tb,false_positive_bp=pb-mb,false_negative_bp=tb-mb,precision=p,recall=r,f1=f))
    return rows
