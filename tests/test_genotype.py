"""Tests on synthetic examples. Do not infer study phenotype-genotype agreement."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from phenotype_gene_crosscheck import crosscheck
from gene_frequency_by_esbl import summarize

class CountsTest(unittest.TestCase):
    def test_requires_coverage_to_infer_negative_gene_calls(self):
        cohort=['1','2','3','4']
        p=[{'isolate_id':'1','laboratory_esbl_phenotype':'Positive'},
           {'isolate_id':'2','laboratory_esbl_phenotype':'Positive'},
           {'isolate_id':'3','laboratory_esbl_phenotype':'Negative'},
           {'isolate_id':'4','laboratory_esbl_phenotype':'Negative'}]
        calls=[{'isolate_id':'1','gene':'blaCTX-M-1','status':'verified','call_class':'acquired'},
          {'isolate_id':'3','gene':'tet(A)','status':'verified','call_class':'acquired'}]
        coverage=[{'isolate_id':'1','annotation_complete':'Yes'},
                  {'isolate_id':'2','annotation_complete':'Yes'},
                  {'isolate_id':'3','annotation_complete':'Yes'},
                  {'isolate_id':'4','annotation_complete':'No'}]
        rows=list(crosscheck(cohort,p,calls,coverage))
        positives=[x for x in rows if x['phenotypic_group']=='ESBL-positive' and x['blaCTXM_carriage']=='CTX-M-present']
        self.assertEqual(positives[0]['isolates_n'],1)
        self.assertEqual(positives[0]['group_denominator_n'],2)
        negatives=[x for x in rows if x['phenotypic_group']=='ESBL-negative' and x['blaCTXM_carriage']=='CTX-M-not-detected']
        self.assertEqual(negatives[0]['group_denominator_n'],1) # isolate 4 incomplete
        freq=list(summarize(cohort,p,calls,coverage))
        tet=[x for x in freq if x['gene']=='tet(A)'][0]
        self.assertEqual(tet['ESBL_negative_n'],1)
        self.assertEqual(tet['ESBL_negative_N'],1)

if __name__=='__main__':unittest.main()
