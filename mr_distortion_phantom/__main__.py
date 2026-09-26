import argparse,csv,json
from pathlib import Path
import numpy as np
from .core import evaluate

def reproduce(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    source=Path(__file__).parent/'data';groups={}
    with (source/'landmarks.csv').open() as f:
        for row in csv.DictReader(f):groups.setdefault(row['protocol']+'_'+row['state'],[]).append(row)
    expected=json.loads((source/'expected_metrics.json').read_text());results={}
    for key,rows in groups.items():
        ct=np.array([[float(r['ct_'+a+'_mm']) for a in 'LPS'] for r in rows]);mr=np.array([[float(r['mr_'+a+'_mm']) for a in 'LPS'] for r in rows]);fit=np.array([r['role']=='registration' for r in rows])
        stats,v,R,t=evaluate(ct,mr,fit)
        for k in stats:
            if abs(stats[k]-expected[key][k])>1e-9:raise AssertionError('Published numeric regression failed: '+key+' '+k)
        results[key]=stats
    (output/'metrics.json').write_text(json.dumps(results,indent=2))
    print(f'Recomputed {len(results)} series from measured landmarks; all metrics match at 1e-9 mm.')
    return results

def main():
    p=argparse.ArgumentParser(description='Independent phantom research analysis')
    sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('reproduce');r.add_argument('--output',default='results')
    r=sub.add_parser('figures');r.add_argument('--output',default='.')
    a=sub.add_parser('analyze');a.add_argument('--ct',required=True);a.add_argument('--original',required=True);a.add_argument('--corrected',required=True);a.add_argument('--ct-crop',nargs=6,type=int);a.add_argument('--output',required=True);a.add_argument('--phantom-only',action='store_true',required=True)
    args=p.parse_args()
    if args.command=='reproduce':reproduce(args.output)
    elif args.command=='figures':
        from .figures import build
        build(Path(args.output))
    else:
        from .dicom import analyze
        analyze(args)
if __name__=='__main__':main()
