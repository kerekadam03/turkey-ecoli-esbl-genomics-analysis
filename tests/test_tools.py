"""Synthetic input behavioural unit tests. They do NOT validate paper results."""
from __future__ import annotations
import csv
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import fasta_qc,analysis_esbl,reconcile_registry,resistome_summary,mge_context_summary,sampling_summary,audit_pangenome_phylogeny
from common import wilson_interval

class ToolsTest(unittest.TestCase):
    def test_fasta_qc_length_n50_and_duplicates(self):
        seqs=[('ctg1','ACGT'*75),('ctg2','N'*199)]
        r=fasta_qc.summarize(seqs,'12_spades_megahit.gam.fasta','demo.zip')
        self.assertEqual(r['assembly_size_bp'],499)
        self.assertEqual(r['n50_bp'],300)
        self.assertEqual(r['contigs_lt_200'],1)
        self.assertEqual(r['gc_percent_acgt'],'50.0000')
        self.assertEqual(r['isolate_id'],'12')
        r=fasta_qc.summarize([('x','A'*200),('x','G'*200)],'test.fna')
        self.assertEqual(r['duplicated_ids'],1)

    def test_zip_qc(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.zip'
            with zipfile.ZipFile(p,'w') as z:z.writestr('part/58.fasta','>x\n'+'A'*250+'\n>y\n'+'C'*400+'\n')
            rows=list(fasta_qc.scan_paths([p]))
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['isolate_id'],'58')
            self.assertEqual(rows[0]['n50_bp'],400)

    def test_mic_synergy_independently_of_lab_call(self):
        row={'isolate_id':'X','cefotaxime_mic_ug_ml':'8','cefotaxime_clavulanate_mic_ug_ml':'1',
           'initial_esbl_screen_positive':'No','laboratory_esbl_phenotype':'Negative'}
        r=analysis_esbl.calculate(row)
        self.assertEqual(r['cefotaxime_synergy_ge8'],'Yes')
        self.assertEqual(r['screen_plus_synergy'],'Not confirmed using supplied CTX screen and synergy data')
        self.assertEqual(r['laboratory_esbl_phenotype'],'Negative')
        self.assertEqual(r['fold_reduction'],'8.00000')
        row['cefotaxime_mic_ug_ml']='>=16'
        self.assertEqual(analysis_esbl.calculate(row)['cefotaxime_synergy_ge8'],'Not evaluable')

    def test_unique_registry_requires_explicit_cohort(self):
        bios=[{'accession':'SAMN1','sample_name':'Ecoli_Turkey_1','bioproject_accession':'PRJNA1475617','organism':'Escherichia coli'},
              {'accession':'SAMN2','sample_name':'Ecoli_Turkey_2','bioproject_accession':'PRJNA1475617','organism':'Escherichia coli'}]
        genome=[{'biosample_accession':'SAMN1','genome_acc':'JCBX','status':'Processed','filename':'1.fasta'}]
        rows,_=reconcile_registry.build(bios,genome)
        self.assertEqual(rows[0]['manuscript_cohort'],'Not verified')
        with self.assertRaises(ValueError):reconcile_registry.build(bios,genome,['1','1'],expected_cohort=2)
        explicit,_=reconcile_registry.build(bios,genome,['1'],expected_cohort=1)
        self.assertEqual(explicit[0]['manuscript_cohort'],'Yes')
        self.assertEqual(explicit[1]['manuscript_cohort'],'No')

    def test_resistome_distinguishes_categories(self):
        calls=[{'isolate_id':'1','gene':'blaCTX-M-1','call_class':'acquired','status':'verified','plasmid_predicted':'Yes','mge_proximal':'No'},
               {'isolate_id':'1','gene':'mdfA','call_class':'intrinsic','status':'verified','plasmid_predicted':'No','mge_proximal':'No'}]
        out=resistome_summary.summarize(calls,['1','2'])
        self.assertEqual(out[0]['n_acquired'],1)
        self.assertEqual(out[0]['n_intrinsic'],1)
        self.assertEqual(out[0]['n_distinct_verified_gene_labels'],2)
        self.assertEqual(out[1]['n_distinct_verified_gene_labels'],0)
        self.assertIn('does NOT imply',out[1]['interpretation'])

    def test_separate_MGE_carriers_from_hits(self):
        rows=[{'isolate_id':'1','gene':'blaCTX-M-1','contig_id':'abc','plasmid_predicted':'Yes','mge_element':'IS26','mge_hit_id':'hit1','distance_bp':'5'},
              {'isolate_id':'1','gene':'blaCTX-M-1','contig_id':'abc','plasmid_predicted':'Yes','mge_element':'IS26','mge_hit_id':'hit2','distance_bp':'200'},
              {'isolate_id':'2','gene':'blaCTX-M-27','contig_id':'zz','plasmid_predicted':'Unknown','mge_element':'','mge_hit_id':'','distance_bp':''}]
        out=mge_context_summary.summarize(rows)
        self.assertEqual(len(out),2)
        self.assertEqual(out[0]['mge_hits_n'],2)
        self.assertEqual(sum(x['has_mge_proximity']=='Yes' for x in out),1)

    def test_sampling_intervals_and_grouping(self):
        m=[{'isolate_id':'A','sample_id':'s1','sampling_site':'Cloacal','farm_id':'f1'},
           {'isolate_id':'B','sample_id':'s2','sampling_site':'Cloacal','farm_id':'f1'},
           {'isolate_id':'C','sample_id':'s3','sampling_site':'Tracheal','farm_id':'f2'}]
        r=sampling_summary.compute(m,{'A':True,'B':False,'C':None})
        clo=[x for x in r if x['stratum']=='sampling_site' and x['category']=='Cloacal'][0]
        self.assertEqual(clo['known_esbl_status_n'],2)
        self.assertEqual(clo['esbl_positive_n'],1)
        l,u=wilson_interval(1,2)
        self.assertLess(l,.5);self.assertGreater(u,.5)

    def test_panaroo_threshold_boundary_and_alignment(self):
        calls=[{'No. isolates':'100'},{'No. isolates':'95'},{'No. isolates':'50'},{'No. isolates':'10'}]
        result=audit_pangenome_phylogeny.gene_categories(calls,100)
        self.assertEqual(result['strict_core_99_to_100'],1)
        self.assertEqual(result['soft_core_95_to_lt99'],1)
        self.assertEqual(result['shell_15_to_lt95'],1)
        self.assertEqual(result['cloud_lt15'],1)
        with tempfile.TemporaryDirectory() as d:
            a=Path(d)/'a.fa';a.write_text('>a\nACGT\n>b\nACGT\n')
            align=audit_pangenome_phylogeny.alignment_dimensions(a)
            self.assertEqual(align['alignment_bp'],4)
            self.assertTrue(align['lengths_equal'])

if __name__=='__main__':unittest.main()
