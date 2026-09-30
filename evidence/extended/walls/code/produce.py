from pathlib import Path
import argparse,time
import walls,components
from core import locate_parent,dump
p=argparse.ArgumentParser();p.add_argument('--parent',type=Path);p.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'results');args=p.parse_args()
parent,tmp=locate_parent(args.parent);args.out.mkdir(exist_ok=True,parents=True);t=time.time()
try:
 w=walls.run(parent,args.out);c=components.run(parent,args.out);b=components.binary_example(args.out)
 dump(args.out/'production.json',{'status':'PASS','mixed_wall_count':len(w),'component_counts':{str(x['k']):x['components'] for x in c},'binary_rc_mds':b['rc_mds_count']})
 print('PRODUCE_PASS seconds',round(time.time()-t,3),flush=True)
finally:
 if tmp:tmp.cleanup()
