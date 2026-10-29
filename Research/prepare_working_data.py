"""Copy archived datasets into a new working directory and decompress text inputs.

Original gzip files (including VCFs and plotting snapshots) are retained.
Only .csv.gz/.tsv.gz/.txt.gz receive additional uncompressed working copies.
"""
from pathlib import Path
import argparse,gzip,shutil,csv,os

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();source=args.source.resolve();output=args.output.resolve()
    if not source.is_dir():parser.error('Source datasets directory does not exist')
    if output.exists() or output.is_relative_to(source):parser.error('Use a new output directory outside the source archive')
    files=[p for p in source.rglob('*') if p.is_file()]
    names={p.relative_to(source).as_posix() for p in files}
    for name in names:
        if name.endswith(('.csv.gz','.tsv.gz','.txt.gz')) and name[:-3] in names:
            parser.error('Archive has both compressed and uncompressed versions: '+name)
    output.mkdir(parents=True)
    for p in files:
        dest=output/p.relative_to(source);dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,dest)
        if dest.name.endswith(('.csv.gz','.tsv.gz','.txt.gz')):
            with gzip.open(dest,'rb') as fi,dest.with_suffix('').open('xb') as fo:
                shutil.copyfileobj(fi,fo,4194304)
    aliases=output/'input_aliases.csv'
    if aliases.is_file():
        with aliases.open(encoding='utf-8-sig',newline='') as handle:
            for row in csv.DictReader(handle):
                source_file=(output/row['source']).resolve()
                target=(output/row['target']).resolve()
                if not source_file.is_relative_to(output) or not target.is_relative_to(output):
                    raise ValueError('Input alias escapes the working directory')
                if not source_file.is_file() or target.exists():raise ValueError('Invalid input alias: '+str(target))
                target.parent.mkdir(parents=True,exist_ok=True)
                try:os.link(source_file,target)
                except OSError:shutil.copy2(source_file,target)
    print('Prepared working datasets:',output)

if __name__=='__main__':main()
