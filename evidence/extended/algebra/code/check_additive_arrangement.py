"""Exact three-magnitude Freiman-2 arrangement in 0<x<y<z.

This counts the additive-floor arrangement, not terminal chambers.
"""
from fractions import Fraction as F
from itertools import product, combinations
from functools import reduce
from math import gcd
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]

def require(ok, name):
    if not ok: raise AssertionError(name)

def main():
    normals=set()
    for c in product(range(-4,5),repeat=3):
        norm=sum(abs(x) for x in c)
        if not norm or reduce(gcd,map(abs,c))!=1: continue
        if norm>(4 if norm%2==0 else 2): continue
        if next(x for x in c if x)<0: continue
        suffix=[sum(c[j:]) for j in range(3)]
        if min(suffix)<0<max(suffix): normals.add(c)
    normals=sorted(normals)
    expected=[(0,3,-1),(1,-2,1),(1,2,-1),(2,1,-1),(3,-1,0),(3,0,-1)]
    require(normals==expected,'complete primitive wall enumeration')
    points={}
    for a,b in combinations(normals,2):
        det=a[1]*b[2]-a[2]*b[1]
        if not det: continue
        u=F(-a[0]*b[2]+a[2]*b[0],det)
        v=F(-a[1]*b[0]+a[0]*b[1],det)
        if u>1 and v>u:
            points[u,v]=sum(c[0]+c[1]*u+c[2]*v==0 for c in normals)
    require(points=={(F(2),F(3)):2,(F(3),F(5)):3,(F(3),F(7)):2,(F(3),F(9)):2},'intersection lattice')
    chambers=1+len(normals)+sum(m-1 for m in points.values())
    edges=len(normals)+sum(points.values())
    cells={}
    for u in (F(3,2),F(5,2),F(4)):
        heights=sorted({u,*(-F(c[0]+c[1]*u,c[2]) for c in normals if c[2] and -F(c[0]+c[1]*u,c[2])>u)})
        heights.append(heights[-1]+2)
        for a,b in zip(heights,heights[1:]):
            v=(a+b)/2
            values=[c[0]+c[1]*u+c[2]*v for c in normals]
            require(all(values),'sample off all walls')
            signs=tuple(1 if x>0 else -1 for x in values)
            cells[signs]={'signs':signs,'u':str(u),'v':str(v)}
    require(chambers==12 and edges==15 and len(cells)==12,'independent chamber witnesses')
    out={'status':'PASS','scope':'rational additive floor, not terminal arrangement','check_count':3,
         'normals':normals,'interior_intersections':[{'u':str(u),'v':str(v),'multiplicity':m} for (u,v),m in sorted(points.items())],
         'open_chambers':chambers,'one_dimensional_faces':edges,'zero_dimensional_faces':len(points),
         'representative_sign_cells':[cells[c] for c in sorted(cells)]}
    (ROOT/'results/three_magnitude_arrangement.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS additive arrangement: 6 walls, 12 chambers, 15 edges, 4 vertices')
if __name__=='__main__':main()
