import argparse,json
from pathlib import Path
from . import search,verify_certificate
from .workflow import analyze

def main():
    p=argparse.ArgumentParser(description='Exact graph Gram witnesses; failure is inconclusive.')
    sub=p.add_subparsers(dest='command',required=True)
    v=sub.add_parser('verify');v.add_argument('file')
    for name in ('search','analyze'):
        a=sub.add_parser(name);a.add_argument('graph6');a.add_argument('--rank',required=True,type=int)
        a.add_argument('--trials',type=int,default=100);a.add_argument('--seed',type=int,default=0)
        a.add_argument('--backend',choices=['auto','sage','sympy'],default='auto');a.add_argument('--output')
        if name=='search':a.add_argument('--negative',type=int,default=0)
    args=p.parse_args()
    if args.command=='verify':
        obj=json.loads(Path(args.file).read_text())
        if isinstance(obj,dict) and 'certificate' in obj:
            if obj['certificate'] is None:raise ValueError('no certificate in exhausted search result')
            obj=obj['certificate']
        cs=obj if isinstance(obj,list) else [obj]
        print(json.dumps([verify_certificate(c) for c in cs],indent=2));return
    kw=dict(trials=args.trials,seed=args.seed,backend=args.backend)
    if args.command=='analyze':result=analyze(args.graph6,args.rank,**kw)
    else:
        if args.negative:kw.update(order_policy='eligible',coefficients=(-1,0,1))
        result=search(args.graph6,args.rank,signature=(args.rank-args.negative,args.negative),**kw)
    text=json.dumps(result,indent=2)
    if args.output:Path(args.output).write_text(text+'\n')
    else:print(text)
if __name__=='__main__':main()
