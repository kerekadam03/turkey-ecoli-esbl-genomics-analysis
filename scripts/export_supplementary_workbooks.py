#!/usr/bin/env python3
"""Create formatted XLSX workbooks from checked, normalized S1--S7 TSV tables.

New reconstruction. Requires optional ``artifact_tool``; this does not recover
missing original assays. The origin of each value must be documented separately.
No empty field is filled by inference.  A matching *_provenance.tsv is advised.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from common import table


def create(input_path, output):
    try:
        from artifact_tool import Workbook,SpreadsheetFile
    except ImportError as exc:
        raise RuntimeError('This optional exporter requires artifact_tool') from exc
    src=Path(input_path)
    import csv
    with src.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.reader(f,delimiter='\t')
        lines=list(reader)
    if not lines or not lines[0]:raise ValueError('TSV must contain headings')
    cols=len(lines[0]);rows=len(lines)
    if any(len(r)!=cols for r in lines):raise ValueError('Jagged TSV; refusing to silently truncate rows')
    wb=Workbook.create()
    sheet=wb.worksheets.add(src.stem[:30])
    sheet.get_range_by_indexes(0,0,rows,cols).values=lines
    sheet.get_range_by_indexes(0,0,1,cols).format={
      'fill':'#183B57','font':{'bold':True,'color':'#FFFFFF','size':11},
      'row_height':42,'wrap_text':True,'vertical_alignment':'center'}
    sheet.get_range_by_indexes(0,0,rows,cols).format.column_width=22
    for i in range(min(cols,4)):
        sheet.get_range_by_indexes(0,i,rows,1).format.column_width=24
    sheet.freeze_panes.freeze_rows(1)
    SpreadsheetFile.export_xlsx(wb).save(str(output))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',nargs='+',required=True,help='Verified S1.tsv to S7.tsv')
    p.add_argument('--outdir',required=True)
    a=p.parse_args(argv)
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    for src in a.inputs:
        dest=out/(Path(src).stem+'.xlsx')
        create(src,dest)
        print('Exported',dest)

if __name__=='__main__':main()
