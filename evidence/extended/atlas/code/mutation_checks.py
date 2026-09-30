"""Typed semantic mutations, including coherently rehashed wrong certificates."""
from pathlib import Path
import copy,json,struct
import verify_atlas as v
ROOT=Path(__file__).resolve().parents[1]

def main():
    A=json.loads((ROOT/'results/atlas_k8.json').read_text());m=v.get_models(8);v.verify_catalog(A,m)
    row=A['rows'][100];c=v.compute(m,row['ratio']);cid=row['id'];bits=(ROOT/f'results/certificates/{cid}.bits').read_bytes();levels=(ROOT/f'results/certificates/{cid}.levels').read_bytes()
    v.verify_certificate(row,m,c,bits,levels);records=[]
    def reject(name,expected,fn):
        try:fn()
        except v.VerificationError as ex:
            v.chk(str(ex)==expected,'MUTATION_WRONG_REASON_'+name);records.append({'mutation':name,'output':str(ex),'status':'REJECTED'})
        else:raise RuntimeError('MUTATION_ACCEPTED_'+name)
    def badcat(name,expected,alter):
        B=copy.deepcopy(A);alter(B);reject(name,expected,lambda:v.verify_catalog(B,m))
    badcat('omit_one_critical_root','ROOT_COMPLETENESS_ORDER',lambda b:b['roots'].pop())
    badcat('swap_first_two_roots','ROOT_COMPLETENESS_ORDER',lambda b:b['roots'].__setitem__(slice(0,2),b['roots'][1::-1]))
    badcat('replace_by_admissible_field_nonroot','ROOT_COMPLETENESS_ORDER',lambda b:b['roots'][0].__setitem__('value',[6,0,5]))
    badcat('omit_point_stratum','CELL_COVERAGE',lambda b:b['rows'].pop(1))
    badcat('sample_on_domain_boundary','INTERVAL_SAMPLE',lambda b:b['rows'][0].__setitem__('ratio',[1,0,1]))
    badcat('inactive_midpoint_witness','ACTIVE_WITNESS',lambda b:b['roots'][0].__setitem__('active_midpoint_word_code',0))
    idx=next(i for i,r in enumerate(A['roots']) if 'lower_lambda' in r['families'])
    badcat('erase_lambda_resonance_tag','LAMBDA_ROOT_TAG',lambda b:b['roots'][idx]['families'].remove('lower_lambda'))
    badcat('alter_embedding_multiplicity','EMBEDDING_SPECTRUM',lambda b:b['embedding_spectrum'][0].__setitem__('accidental_words',33))
    badcat('field_gate_as_realizability','EMBEDDING_SPECTRUM',lambda b:b['embedding_spectrum'][0].__setitem__('value',[2,0,1]))
    badcat('structural_embedding_witness','EMBEDDING_WITNESS',lambda b:b['embedding_spectrum'][0].__setitem__('witness_code',0))
    def cert(name,expected,alter):
        r=copy.deepcopy(row);b=bytearray(bits);ls=list(struct.unpack('<'+'H'*(len(levels)//2),levels));alter(r,b,ls)
        lb=struct.pack('<'+'H'*len(ls),*ls)
        # Rehash mutated bytes deliberately, so the semantic predicate must reject.
        r['bitset_sha256']=v.bh(bytes(b));r['level_sha256']=v.bh(lb)
        reject(name,expected,lambda:v.verify_certificate(r,m,c,bytes(b),lb))
    cert('coherently_rehashed_terminal_bit','SELECTOR_MISMATCH',lambda r,b,l:b.__setitem__(0,b[0]^1))
    cert('lowered_longest_path','LEVEL_BOUNDS',lambda r,b,l:r.__setitem__('lower_longest_edges',140))
    edge=next((w for i,w in enumerate(m['ww']) if c['term'][i]==0 and v.encode(w[1:])!=v.encode(w[:-1])))
    head=v.encode(edge[1:]);cert('coherently_rehashed_edge_level','LEVEL_EDGE_INEQUALITY',lambda r,b,l:l.__setitem__(head,0))
    cert('repeat_path_vertex','PATH_SIZE_OR_REPEAT',lambda r,b,l:r['longest_path_lower_vertex_codes'].__setitem__(1,r['longest_path_lower_vertex_codes'][0]))
    cert('wrong_zero_count','ZERO_EMBEDDING_COUNT',lambda r,b,l:r.__setitem__('zero_embedding_words',257))
    cert('wrong_minimum_cardinality','CARDINALITY',lambda r,b,l:r.__setitem__('selected_count',8229))
    cert('wrong_positive_block_length','BLOCK_LENGTHS',lambda r,b,l:r['block_lengths'].__setitem__('1',28))
    cert('confuse_upper_and_lower_edge_units','LINE_GRAPH_UNITS',lambda r,b,l:r.__setitem__('upper_longest_edges',141))
    v.chk(len(records)==18,'MUTATION_COUNT');result={'status':'PASS','mutations_rejected':len(records),'mutations':records}
    (ROOT/'results/mutation_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    for r in records:print(r['mutation']+': '+r['output'])
    print('PASS: 18/18 semantic mutations rejected, including coherent rehashes')
if __name__=='__main__':main()
