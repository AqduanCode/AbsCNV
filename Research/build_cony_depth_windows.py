"""Build autosomal fixed windows with GC and non-ambiguous fractions for CONY."""
import re,sys

def emit(chrom,start0,seq,bed,meta):
    seq=seq.upper();width=len(seq)
    if width==0:return
    end0=start0+width
    acgt=sum(seq.count(base) for base in 'ACGT');gc=seq.count('G')+seq.count('C')
    bed.write(f'{chrom}\t{start0}\t{end0}\n')
    meta.write(f'{chrom}\t{start0}\t{end0}\t{width}\t{acgt/width:.8f}\t{gc/width:.8f}\n')

def main():
    if len(sys.argv)!=5:raise SystemExit('Usage: build_cony_depth_windows.py REF.fa WINDOW OUT.bed OUT.meta.tsv')
    fasta,window_text,bed_path,meta_path=sys.argv[1:];window=int(window_text)
    autosome=re.compile(r'^chr(?:[1-9]|1[0-9]|2[0-2])$')
    with open(fasta) as handle,open(bed_path,'w') as bed,open(meta_path,'w') as meta:
        meta.write('chrom\tstart0\tend0\twidth\tnonAmb\tGC\n')
        chrom=None;keep=False;buf='';start0=0
        for raw in handle:
            if raw.startswith('>'):
                if keep and buf:emit(chrom,start0,buf,bed,meta)
                chrom=raw[1:].split()[0];keep=autosome.fullmatch(chrom) is not None
                buf='';start0=0;continue
            if not keep:continue
            buf+=raw.strip()
            while len(buf)>=window:
                emit(chrom,start0,buf[:window],bed,meta);buf=buf[window:];start0+=window
        if keep and buf:emit(chrom,start0,buf,bed,meta)

if __name__=='__main__':main()
