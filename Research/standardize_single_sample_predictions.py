"""Rebuild the five unmerged evaluation tables from native caller outputs."""
from pathlib import Path
import argparse,csv,math,re,shutil

def write(path,rows,fields):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();root=a.dataset;out=a.output_dir;out.mkdir(parents=True,exist_ok=True)
    native=root/'native'
    with (root/'output/copynumber_background.csv').open(encoding='utf-8-sig') as f:
        reader=csv.DictReader(f);fields=reader.fieldnames;rows=[]
        for r in reader:
            try:cn=float(r['Sample_A'])
            except (ValueError,TypeError):continue
            if math.isfinite(cn) and cn!=2:rows.append(r)
    write(out/'MyMethod_raw.tsv',rows,fields)
    fields=['type','chrom','start','end','length','normalized_RD','e_val1','e_val2','e_val3','e_val4','q0'];rows=[]
    with (native/'CNVnator/Sample_A.cnvnator.calls.tsv').open() as f:
        for line in f:
            x=line.split()
            if not x or x[0] not in {'deletion','duplication'}:continue
            m=re.fullmatch(r'(chr[^:]+):(\d+)-(\d+)',x[1])
            if m is None:raise ValueError(x[1])
            row=[x[0],m[1],m[2],m[3],*x[2:]]
            if len(row)!=len(fields):raise ValueError('Unexpected CNVnator columns')
            rows.append(dict(zip(fields,row)))
    write(out/'CNVnator_raw.tsv',rows,fields)
    fields=['chrom','start','end','copy_number','status'];rows=[]
    with (native/'ControlFREEC/Sample_A.sorted.bam_CNVs').open() as f:
        for line in f:
            x=line.split()
            if len(x)>=5 and x[4].lower() in {'gain','loss'}:rows.append(dict(zip(fields,x[:5])))
    write(out/'ControlFREEC_raw.tsv',rows,fields)
    cony=native/'CONY/Sample_A.CONY.Single.10kb.BaselinePloidy2.CNV.tsv'
    shutil.copy2(cony,out/'CONY_raw.tsv')
    fields=['segStart','segEnd','state','length','copyEst','segMedian','zstat','chrom','posStart','posEnd'];rows=[]
    for n in range(1,23):
        with (native/'RDXplorer'/f'chr{n}.sum.filt').open(encoding='utf-8-sig') as f:
            reader=csv.DictReader(f,delimiter='\t')
            rows.extend({k:r[k] for k in fields} for r in reader if int(float(r['state'])) in {1,3})
    write(out/'RDXplorer_raw.tsv',rows,fields)
    print('Wrote five unmerged caller tables:',out)

if __name__=='__main__':main()
