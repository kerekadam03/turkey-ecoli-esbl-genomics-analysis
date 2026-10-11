#!/usr/bin/env python3
"""Regenerate descriptive figures from verified flat-file inputs, not hardcoded counts.

Figures are new reference implementations; they are not the original plots used
in the manuscript. Requires matplotlib. Each figure is standalone and strictly
checks the required columns. The status of breakpoint/ECOFF interpretations
must be supplied explicitly in the input TSV.
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path
from collections import defaultdict,Counter
from common import table,tri_bool


def figure1(data,out):
    """Stacked S/I/R by drug. Does not plot ECOFF or cefotaxime-clavulanate as S/I/R."""
    import matplotlib.pyplot as plt
    c=table(data,['drug','category','n','N','interpretation_type'])
    tmp=defaultdict(dict);Ns={}
    categories=['Susceptible','Intermediate','Resistant']
    for r in c:
        if r['interpretation_type'].strip().lower()!='clinical':continue
        if 'clavulanate' in r['drug'].lower() and 'cefotaxime' in r['drug'].lower():continue
        drug=r['drug'];cat=r['category']
        if cat not in categories:raise ValueError('Nonclinical or unexpected category in Figure 1: '+cat)
        n=int(r['n']);N=int(r['N']);
        if not 0<=n<=N:raise ValueError('Invalid n/N')
        if cat in tmp[drug]:raise ValueError('Duplicate drug/category')
        tmp[drug][cat]=n
        if drug in Ns and Ns[drug]!=N:raise ValueError('Inconsistent denominator')
        Ns[drug]=N
    if not tmp:raise ValueError('No clinical drugs to plot')
    for drug,count in tmp.items():
        if set(count)!=set(categories) or sum(count.values())!=Ns[drug]:
            raise ValueError(f'Incomplete categories or incorrect sum for {drug}')
    agents=list(tmp);labels=agents
    vals={cat:[100*tmp[d][cat]/Ns[d] for d in agents] for cat in categories}
    palette={'Susceptible':'#0C985B','Intermediate':'#F0B323','Resistant':'#CB3128'}
    fig,ax=plt.subplots(figsize=(max(10,len(agents)*0.82),6))
    bottoms=[0.0]*len(agents)
    for cat in categories:
        v=vals[cat]
        ax.bar(range(len(agents)),v,bottom=bottoms,color=palette[cat],label=cat,edgecolor='white',linewidth=.35)
        for i,percent in enumerate(v):
            if percent>=5.5:
                ax.text(i,bottoms[i]+percent/2,f'{percent:.1f}%',color='white',
                     fontweight='bold',fontsize=8,ha='center',va='center')
            elif 0<percent<2 and bottoms[i]+percent>=95:
                ax.text(i,101,f'{percent:.1f}%',ha='center',va='bottom',fontsize=8,clip_on=False)
        bottoms=[a+b for a,b in zip(bottoms,v)]
    ax.set_xticks(range(len(agents)));ax.set_xticklabels(labels,rotation=45,ha='right')
    ax.set_ylim(0,106);ax.set_ylabel('Proportion of isolates (%)')
    ax.set_xlabel('Antimicrobial agent',labelpad=12)
    ax.grid(axis='y',color='#DDE4EA',alpha=.7)
    ax.set_axisbelow(True)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.35),ncol=3,frameon=False)
    fig.subplots_adjust(bottom=.42,top=.93,left=.08,right=.98)
    fig.savefig(out,dpi=600)
    plt.close(fig)


def figure2(data,out):
    import matplotlib.pyplot as plt
    c=table(data,['isolate_id','log2_fold_reduction'])
    values=[float(x['log2_fold_reduction']) for x in c if x['log2_fold_reduction']]
    if not values:raise ValueError('No paired MIC-derived fold ratios')
    fig,ax=plt.subplots(figsize=(8,5));ax.hist(values,bins=25,color='#2C6DA6',edgecolor='white')
    ax.axvline(3,color='#B5553F',linestyle='--',label='8-fold cefotaxime MIC reduction')
    ax.set(xlabel='log2(cefotaxime MIC / cefotaxime-clavulanate MIC)',ylabel='Number of isolates')
    ax.legend(frameon=False);fig.tight_layout();fig.savefig(out,dpi=600);plt.close(fig)


def figure3(data,out):
    """Prevalence of named gene labels from a verified per-gene n/N table.

    The input must distinguish acquired genes from intrinsic / efflux determinants;
    category assignment is not guessed by this plotting code.
    """
    import matplotlib.pyplot as plt
    rows=table(data,['gene','carriers_n','denominator_n','call_class'])
    if not rows:raise ValueError('No gene counts to plot')
    counts=[]
    for r in rows:
        n=int(r['carriers_n']);N=int(r['denominator_n'])
        if N<=0 or n<0 or n>N:raise ValueError(f'Invalid n/N for {r["gene"]}')
        counts.append((r['gene'],100*n/N,r['call_class']))
    if len({x[0] for x in counts})!=len(counts):raise ValueError('Repeated gene name')
    counts.sort(key=lambda x:x[1],reverse=True)
    fig,ax=plt.subplots(figsize=(9,max(5,len(counts)*.35)))
    ypos=list(range(len(counts)))
    ax.barh(ypos,[x[1] for x in counts],color='#2765A3')
    ax.set_yticks(ypos,[x[0] for x in counts]);ax.invert_yaxis()
    ax.set_xlabel('Percentage of included study isolates (%)')
    ax.set_xlim(0,max(100,max(x[1] for x in counts)*1.07))
    fig.tight_layout();fig.savefig(out,dpi=600);plt.close(fig)


def figure4(data,out):
    import matplotlib.pyplot as plt
    rows=table(data,['isolate_id','blaCTXM_gene_labels','blactxm_on_plasmid_predicted_contig','has_mge_proximity'])
    groups=defaultdict(list)
    for r in rows:
        # one allele per isolate is expected for original manuscript; detect otherwise
        genes=[s.strip() for s in r['blaCTXM_gene_labels'].split(';') if s.strip()]
        if len(genes)!=1:raise ValueError('Cannot plot more than one blaCTX-M allele per isolate without explicit counting rules')
        groups[genes[0]].append(r)
    if not groups:raise ValueError('No blaCTX-M-positive isolates')
    labels=list(groups)
    fields=['Carriers','Plasmid-predicted contig','MGE-proximal']
    values=[[len(groups[k]),sum(tri_bool(x['blactxm_on_plasmid_predicted_contig']) is True for x in groups[k]),
                 sum(x['has_mge_proximity']=='Yes' for x in groups[k])] for k in labels]
    fig,ax=plt.subplots(figsize=(8,5));import numpy as np
    xx=np.arange(len(labels));width=.25
    for j,f in enumerate(fields):ax.bar(xx+(j-1)*width,[v[j] for v in values],width,label=f)
    ax.set_xticks(xx,labels);ax.set_ylabel('Number of isolates');ax.legend(frameon=False)
    fig.tight_layout();fig.savefig(out,dpi=600);plt.close(fig)


def figure6(data,out):
    import matplotlib.pyplot as plt
    vals=defaultdict(list)
    for r in table(data,['isolate_id','laboratory_esbl_phenotype','n_distinct_verified_gene_labels']):
        if r['laboratory_esbl_phenotype'] not in ('Positive','Negative'):continue
        if not r['n_distinct_verified_gene_labels']:continue
        vals[r['laboratory_esbl_phenotype']].append(float(r['n_distinct_verified_gene_labels']))
    if len(vals)!=2:raise ValueError('ESBL positive and negative groups required')
    fig,ax=plt.subplots(figsize=(6,5));
    ax.boxplot([vals['Negative'],vals['Positive']],tick_labels=['ESBL−','ESBL+'])
    ax.set_ylabel('Number of distinct verified gene labels per isolate')
    fig.tight_layout();fig.savefig(out,dpi=600);plt.close(fig)


def figure7(data,out):
    import matplotlib.pyplot as plt
    rows=[r for r in table(data,['stratum','category','known_esbl_status_n','esbl_positive_n']) if r['stratum']=='sampling_site']
    if not rows:raise ValueError('No sampling-site data')
    labels=[r['category'] for r in rows];n=[int(r['esbl_positive_n']) for r in rows];N=[int(r['known_esbl_status_n']) for r in rows]
    if any(d<=0 for d in N):raise ValueError('No denominator for site comparison')
    ys=[100*a/b for a,b in zip(n,N)]
    fig,ax=plt.subplots(figsize=(7,5));bars=ax.bar(labels,ys,color=['#2468AB','#0D9181'][:len(labels)])
    for j,b in enumerate(bars):ax.text(b.get_x()+b.get_width()/2,b.get_height()+.5,f'{n[j]}/{N[j]} ({ys[j]:.1f}%)',ha='center')
    ax.set_ylim(0,max(30,max(ys)*1.3));ax.set_ylabel('ESBL-positive isolates (%)')
    fig.tight_layout();fig.savefig(out,dpi=600);plt.close(fig)

FUNCS={'1':figure1,'2':figure2,'3':figure3,'4':figure4,'6':figure6,'7':figure7}

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--figure',choices=list(FUNCS),required=True)
    p.add_argument('--input',required=True)
    p.add_argument('--output',required=True)
    args=p.parse_args(argv)
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    FUNCS[args.figure](args.input,args.output)
    print('New reference figure generated:',args.output)

if __name__=='__main__':main()
