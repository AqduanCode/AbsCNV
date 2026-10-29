"""Export the final LUMPY support counts as the manuscript summary CSV.

Calls evaluate_lumpy_support.py unchanged and tabulates its nofail result using
the original summary formulas. No new overlap or filtering rule is introduced.
"""
import argparse,csv,json,subprocess,sys
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    raw=subprocess.check_output([sys.executable,str(Path(__file__).with_name('evaluate_lumpy_support.py')),str(a.dataset)],text=True,encoding='utf-8')
    result=next(r for r in map(json.loads,raw.splitlines()) if r['mode']=='nofail')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f)
        w.writerow(['pair','AbsCNV_events','pair_specific_filtered_LUMPY_calls','supported_AbsCNV_events','supported_fraction','supported_percent'])
        for r in result['pairs']:
            n=r['abscnv_events'];k=r['abscnv_events_any'];total=sum(r['pair_unique'].values())
            w.writerow([r['pair'],n,total,k,k/n,100*k/n])
    print(a.output)

if __name__=='__main__':main()
