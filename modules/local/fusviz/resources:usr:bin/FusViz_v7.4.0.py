#!/usr/bin/env python3

"""
MSK/CTI's FUSionVIsualiZer (FusViz) is an in-house tool to visualize fusions in RNAseq context specially tailored for the MSK-TARGET panel.
"""

__author__       = "Juan Blanco Heredia"
__email__        = "blancoj@mskcc.org"
__version__      = "7.4.0"
__status__       = "Production"

import os
import re
import sys
import pysam
import argparse
import colorsys
import subprocess
import numpy as np
import pandas as pd
import urllib.parse
from typing import List
import matplotlib as mpl
from pycirclize import Circos
from pyranges import PyRanges
import matplotlib.pyplot as plt
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.colors as mcolors
from matplotlib.colors import to_rgba
import matplotlib.patches as mpatches
from matplotlib import patheffects as pe
from pycirclize.utils import ColorCycler
mpl.rcParams['pdf.use14corefonts'] = True
plt.rcParams['axes.unicode_minus'] = False
from matplotlib.patches import FancyBboxPatch
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_pdf import PdfPages

def remove_chr(x: str) -> str:
    if pd.isna(x) or x == "":
        return ""
    elif x == 'chrM':
        return "MT"
    return re.sub(r"^(chr|CHR)", "", str(x))

def parse_args():
    parser = argparse.ArgumentParser( description="MSK/CTI's FUSionVIsualiZer (FusViz) is an in-house tool to visualize fusions in RNAseq context specially tailored for the MSK-TARGET panel." )
    parser.add_argument('--color1'                   , default='#5A73B3', help='Primary color defaults to #5A73B3, other options are #FF6B6B, #FFEA0A, #5974B2, #FF7300, #54F7CC, #FA7550')
    parser.add_argument('--color2'                   , default="#60187D", help='Secondary color defaults to #60187D, other options are #8900CC, #4ECDC4, #9650F5, #9683F5, #00D9FF, #339242')
    parser.add_argument('--output'                   , required=True, help='Output PDF file')
    parser.add_argument('--fusions'                  , required=True, help='Path to fusions file (e.g., path/to/Sample.final.cff)')
    parser.add_argument('--minMapQ'                  , type=int, default=20, help='Minimum MAPQ for reads considered in junctions/bridge (default: 20)')
    parser.add_argument('--version'                  , action='version', version='FusViz 7.4.0')
    parser.add_argument('--fontSize'                 , type=float, default=1, help='Font size (default: 1)')
    parser.add_argument('--pdfWidth'                 , type=float, default=11.692, help='PDF width (default: 11.692)')
    parser.add_argument('--cytobands'                , default='', help='Path to cytobands file (optional)')
    parser.add_argument('--pdfHeight'                , type=float, default=8.267, help='PDF height (default: 8.267)')
    parser.add_argument('--alignments'               , default='', help='Path to BAM file (optional)')
    parser.add_argument('--fontFamily'               , default='Helvetica', help='Font family (default: Helvetica)')
    parser.add_argument('--annotation'               , required=True, help='Path to annotation GTF file')
    parser.add_argument('--chromosomes'              , default='', help='Path to chromosomes file (optional)')
    parser.add_argument('--fusionsFormat'            , default='auto', choices=['auto', 'tsv', 'vcf', 'cff'], help='Input format for --fusions (auto-detect by default).')
    parser.add_argument('--sashimiWindow'            , type=int, default=500, help='Bp window each side of the breakpoints to count junctions (default: 5000)')
    parser.add_argument('--proteinDomains'           , default='', help='Path to protein domains GFF3 file (optional)')
    parser.add_argument('--swapDomainColors'         , action='store_true', help='Swap protein domain colors: use adjust_lightness(color2) for left domains and adjust_lightness(color1) for right domains')
    parser.add_argument('--transcriptSelection'      , default='provided', choices=['coverage','provided','canonical'], help='Transcript selection method (default: provided)')
    parser.add_argument('--optimizeDomainColors'     , type=lambda x: x.lower() in ['true','1','t','yes'], default=False, help='Optimize domain colors (default: False)')
    parser.add_argument('--mergeDomainsOverlappingBy', type=float, default=0.9, help='Merge domains overlapping by (default: 0.9)')
    args = parser.parse_args()
    for f in [args.fusions, args.annotation]:
        if not os.path.isfile(f):
            print(f"[ERROR] Required file not found: {f}", file=sys.stderr)
            sys.exit(1)
    for f in [args.alignments, args.cytobands, args.proteinDomains]:
        if f and not os.path.isfile(f):
            print(f"[WARNING] Optional file not found: {f}", file=sys.stderr)
    return args

def parse_vcf_to_df(vcf_path: str) -> pd.DataFrame:
    rows = []
    vf = pysam.VariantFile(vcf_path)
    sample = next(iter(vf.header.samples)) if len(vf.header.samples) > 0 else None
    def _get(info, key, default=''):
        try:
            v = info.get(key, default)
            if isinstance(v, tuple):
                v = ','.join(str(x) for x in v)
            return '.' if v in (None, 'nan', 'NaN') else str(v)
        except Exception:
            return '.'
    for rec in vf:
        info  = rec.info
        chra  = _get(info, 'CHRA', rec.chrom)
        chrb  = _get(info, 'CHRB', chra)
        posa  = _get(info, 'POSA', rec.pos)
        posb  = _get(info, 'POSB', rec.pos)
        genea = _get(info, 'GENEA', '.')
        geneb = _get(info, 'GENEB', '.')
        orientation = _get(info, 'ORIENTATION', 'nan,nan')
        s1, s2 = orientation.split(',') if ',' in orientation else ('nan', 'nan')
        s1 = s1 if s1 in ['+', '-'] else '+'
        s2 = s2 if s2 in ['+', '-'] else '+'
        frame = _get(info, 'FRAME_STATUS', 'unclear').lower()
        if frame in ('inframe', 'in-frame'):
            reading_frame = 'in-frame'
        elif frame in ('frameshift', 'out-of-frame', 'out_of_frame'):
            reading_frame = 'out-of-frame'
        else:
            reading_frame = 'unclear'
        tx_a = _get(info, 'TRANSCRIPT_ID_A', '.')
        tx_b = _get(info, 'TRANSCRIPT_ID_B', '.')
        dv = rv = ffpm = '.'
        if sample and rec.format:
            call = rec.samples[sample]
            dv = str(call.get('DV', '.')) if 'DV' in rec.format else '.'
            rv = str(call.get('RV', '.')) if 'RV' in rec.format else '.'
            ffpm = str(call.get('FFPM', '.')) if 'FFPM' in rec.format else '.'
        rows.append({
            'gene5_chr'                : chra,
            'gene3_chr'                : chrb,
            'gene5_breakpoint'         : posa,
            'gene3_breakpoint'         : posb,
            'gene5_strand'             : s1,
            'gene3_strand'             : s2,
            'gene5_renamed_symbol'     : genea,
            'gene3_renamed_symbol'     : geneb,
            'gene5_tool_annotation'    : 'exon',
            'gene3_tool_annotation'    : 'exon',
            'gene5_transcript_id'      : tx_a,
            'gene3_transcript_id'      : tx_b,
            'Fusion_effect'            : reading_frame,
            'max_split_cnt'            : rv,
            'SVTYPE'                   : _get(info, 'SVTYPE', rec.info.get('SVTYPE', 'BND')),
            'DV'                       : dv,
            'RV'                       : rv,
            'FFPM'                     : ffpm,
            'somatic_flags'            : '.',
        })
    df = pd.DataFrame(rows, dtype=str)
    return df

def parse_tsv_to_df(tsv_path: str) -> pd.DataFrame:
    df_raw = pd.read_csv(tsv_path, sep='\t', dtype=str, low_memory=False)
    first_col = df_raw.columns[0]
    last_col = df_raw.columns[-1]
    def _gene_symbol(val: str) -> str:
        if pd.isna(val):
            return '.'
        v = str(val)
        return v.split('^', 1)[0] if '^' in v else v
    def _bp_triplet(val: str):
        if pd.isna(val):
            return ('.', '0', '+')
        v = str(val)
        parts = v.split(':')
        if len(parts) >= 3:
            chrom, pos, strand = parts[0], parts[1], parts[2]
            return (chrom, pos, strand if strand in ['+', '-'] else '+')
        if len(parts) >= 2:
            chrom, pos = parts[0], parts[1]
            return (chrom, pos, '+')
        return ('.', '0', '+')
    def _frame(val: str) -> str:
        if pd.isna(val):
            return 'unclear'
        v = str(val).strip().lower()
        if 'in-frame' in v or v == 'inframe':
            return 'in-frame'
        if 'frame' in v:
            return 'out-of-frame'
        return 'unclear'
    def _to_int_str(x):
        if pd.isna(x) or str(x).strip() == '':
            return '0'
        try:
            return str(int(float(str(x))))
        except:
            return '0'
    def _max2(a, b):
        ai = int(_to_int_str(a))
        bi = int(_to_int_str(b))
        return str(max(ai, bi))
    rows = []
    if first_col == '#gene1':
        for _, r in df_raw.iterrows():
            if 'action' in r and str(r.get('action', '')).strip().lower() == 'drop':
                continue
            g5 = _gene_symbol(r.get('#gene1', '.'))
            g3 = _gene_symbol(r.get('gene2', '.'))
            l_chr, l_pos, _ = _bp_triplet(r.get('breakpoint1', '.'))
            r_chr, r_pos, _ = _bp_triplet(r.get('breakpoint2', '.'))
            s1 = r.get('strand1(gene/fusion)', '+')
            s2 = r.get('strand2(gene/fusion)', '+')
            l_strand = str(s1).split('/', 1)[0] if not pd.isna(s1) else '+'
            r_strand = str(s2).split('/', 1)[0] if not pd.isna(s2) else '+'
            split1 = r.get('split_reads1', '0')
            split2 = r.get('split_reads2', '0')
            disc = r.get('discordant_mates', '0')
            rows.append({
                'gene5_chr'            : l_chr or '.',
                'gene3_chr'            : r_chr or '.',
                'gene5_breakpoint'     : l_pos or '0',
                'gene3_breakpoint'     : r_pos or '0',
                'gene5_strand'         : l_strand if l_strand in ['+','-'] else '+',
                'gene3_strand'         : r_strand if r_strand in ['+','-'] else '+',
                'gene5_renamed_symbol' : g5 or '.',
                'gene3_renamed_symbol' : g3 or '.',
                'gene5_tool_annotation': 'exon',
                'gene3_tool_annotation': 'exon',
                'gene5_transcript_id'  : r.get('transcript_id1', '.') or '.',
                'gene3_transcript_id'  : r.get('transcript_id2', '.') or '.',
                'cds_left_range'       : '.',
                'cds_right_range'      : '.',
                'Fusion_effect'        : _frame(r.get('reading_frame', '')),
                'max_split_cnt'        : _max2(split1, split2),
                'SVTYPE'               : 'BND',
                'DV'                   : _max2(split1, split2),
                'RV'                   : _to_int_str(disc),
                'FFPM'                 : '.',
                'somatic_flags'        : '.',
            })
    elif last_col == 'somatic_flags':
        for _, r in df_raw.iterrows():
            if 'action' in r and str(r.get('action', '')).strip().lower() == 'drop':
                continue
            fusion = r.get('fusion', '.') or '.'
            if fusion == '.' or '::' not in fusion:
                g5 = '.'
                g3 = '.'
            else:
                parts = str(fusion).split('::', 1)
                g5 = _gene_symbol(parts[0])
                g3 = _gene_symbol(parts[1])
            bp = r.get('breakpoint', '.')
            if pd.isna(bp) or '|' not in str(bp):
                l_chr, l_pos, l_strand = ('.','0','+')
                r_chr, r_pos, r_strand = ('.','0','+')
            else:
                left_bp, right_bp = str(bp).split('|', 1)
                l_chr, l_pos, l_strand = _bp_triplet(left_bp)
                r_chr, r_pos, r_strand = _bp_triplet(right_bp)
            dv_candidates = [
                r.get('JunctionReadCount', None),
                r.get('junction_reads', None),
                r.get('split_reads', None),
                r.get('SR', None),
                r.get('total_support', None),
                r.get('support', None),
            ]
            rv_candidates = [
                r.get('SpanningFragCount', None),
                r.get('spanning', None),
                r.get('PE', None),
            ]
            som_tgs = r.get('somatic_flags', None)
            if som_tgs is None or pd.isna(som_tgs):
                sf = '0'
            elif isinstance(som_tgs, (list, tuple)):
                sf = next((x for x in som_tgs if x is not None and str(x).strip()!=''), '0')
            else:
                sf = str(som_tgs).strip() if str(som_tgs).strip() != '' else '0'
            dv = next((x for x in dv_candidates if x is not None and str(x).strip()!=''), '0')
            rv = next((x for x in rv_candidates if x is not None and str(x).strip()!=''), '.')
            rows.append({
                'gene5_chr'            : l_chr,
                'gene3_chr'            : r_chr,
                'gene5_breakpoint'     : l_pos,
                'gene3_breakpoint'     : r_pos,
                'gene5_strand'         : l_strand,
                'gene3_strand'         : r_strand,
                'gene5_renamed_symbol' : g5,
                'gene3_renamed_symbol' : g3,
                'gene5_tool_annotation': 'exon',
                'gene3_tool_annotation': 'exon',
                'gene5_transcript_id'  : r.get('tx5', '.') or '.',
                'gene3_transcript_id'  : r.get('tx3', '.') or '.',
                'cds_left_range'       : '.',
                'cds_right_range'      : '.',
                'Fusion_effect'        : _frame(r.get('Fusion_effect', '')),
                'max_split_cnt'        : _to_int_str(dv),
                'SVTYPE'               : 'BND',
                'DV'                   : _to_int_str(dv),
                'RV'                   : (_to_int_str(rv) if rv != '.' else '.'),
                'FFPM'                 : '.',
                'somatic_flags'        : sf,
            })
    elif first_col == 'SV_BC':
        def _norm_effect(val: str) -> str:
            if val is None or pd.isna(val): return 'unclear'
            v = str(val).strip().lower()
            if v in ('in-frame','inframe','frame-preserving') or 'in-frame' in v or 'inframe' in v or 'in frame' in v:
                return 'in-frame'
            if v in ('out-of-frame','outframe') or ('frame' in v and 'in' not in v):
                return 'out-of-frame'
            if v in 'antisense' or 'antisense' in v:
                return 'antisense'
            if v in ('unclear','unknown','na','.'):
                return 'unclear'
            return 'unclear'
        for _, r in df_raw.iterrows():
            l_chr_raw = r.get('CHROM', '')
            if pd.isna(l_chr_raw) or not str(l_chr_raw).strip() or str(l_chr_raw).strip() == '.':
                l_chr_raw = r.get('chr1', '')
            if pd.isna(l_chr_raw) or not str(l_chr_raw).strip() or str(l_chr_raw).strip() == '.':
                l_chr = '.'
            else:
                l_chr = remove_chr(str(l_chr_raw).strip())
            r_chr_raw = r.get('CHR2', '')
            if pd.isna(r_chr_raw) or not str(r_chr_raw).strip() or str(r_chr_raw).strip() == '.':
                r_chr_raw = r.get('chr2', '')
            if pd.isna(r_chr_raw) or not str(r_chr_raw).strip() or str(r_chr_raw).strip() == '.':
                r_chr = '.'
            else:
                r_chr = remove_chr(str(r_chr_raw).strip())
            l_pos = _to_int_str(r.get('POS', '0'))
            r_pos = _to_int_str(r.get('END', '0'))
            strands = r.get('STRANDS', None)
            if pd.isna(strands) or not strands or len(str(strands)) < 2:
                l_strand, r_strand = '+', '+'
            else:
                s = str(strands)
                l_strand = s[0] if s[0] in ['+','-'] else '+'
                r_strand = s[1] if s[1] in ['+','-'] else '+'
            sr = r.get('SR', '0')
            pe = r.get('PE', '0')
            fusion_effect_raw = r.get('fusion', None)
            fusion_effect_norm = _norm_effect(fusion_effect_raw)
            rows.append({
                'gene5_chr'            : l_chr,
                'gene3_chr'            : r_chr,
                'gene5_breakpoint'     : l_pos,
                'gene3_breakpoint'     : r_pos,
                'gene5_strand'         : l_strand,
                'gene3_strand'         : r_strand,
                'gene5_renamed_symbol' : r.get('gene1', '.'),
                'gene3_renamed_symbol' : r.get('gene2', '.'),
                'gene5_tool_annotation': 'exon',
                'gene3_tool_annotation': 'exon',
                'gene5_transcript_id'  : '.',
                'gene3_transcript_id'  : '.',
                'cds_left_range'       : '.',
                'cds_right_range'      : '.',
                'Fusion_effect'        : fusion_effect_norm,
                'max_split_cnt'        : _to_int_str(sr),
                'SVTYPE'               : r.get('SVTYPE', 'BND') or 'BND',
                'DV'                   : _to_int_str(sr),
                'RV'                   : _to_int_str(pe),
                'FFPM'                 : '.',
                'somatic_flags'        : '.',
            })
    elif first_col == 'sample':
        for _, r in df_raw.iterrows():
            fusion = r.get('fusion', '.') or '.'
            if fusion == '.' or '::' not in fusion:
                g5 = '.'
                g3 = '.'
            else:
                parts = str(fusion).split('::', 1)
                g5 = _gene_symbol(parts[0])
                g3 = _gene_symbol(parts[1])
            bp = r.get('breakpoint', '.')
            if pd.isna(bp) or '|' not in str(bp):
                l_chr, l_pos, l_strand = ('.','0','+')
                r_chr, r_pos, r_strand = ('.','0','+')
            else:
                left_bp, right_bp = str(bp).split('|', 1)
                l_chr, l_pos, l_strand = _bp_triplet(left_bp)
                r_chr, r_pos, r_strand = _bp_triplet(right_bp)
            dv_candidates = [
                r.get('JunctionReadCount', None),
                r.get('junction_reads', None),
                r.get('split_reads', None),
                r.get('SR', None),
                r.get('total_support', None),
                r.get('support', None),
            ]
            rv_candidates = [
                r.get('SpanningFragCount', None),
                r.get('spanning', None),
                r.get('PE', None),
            ]
            dv = next((x for x in dv_candidates if x is not None and str(x).strip()!=''), '0')
            rv = next((x for x in rv_candidates if x is not None and str(x).strip()!=''), '.')
            rows.append({
                'gene5_chr'            : l_chr,
                'gene3_chr'            : r_chr,
                'gene5_breakpoint'     : l_pos,
                'gene3_breakpoint'     : r_pos,
                'gene5_strand'         : l_strand,
                'gene3_strand'         : r_strand,
                'gene5_renamed_symbol' : g5,
                'gene3_renamed_symbol' : g3,
                'gene5_tool_annotation': 'exon',
                'gene3_tool_annotation': 'exon',
                'gene5_transcript_id'  : r.get('tx5', '.') or '.',
                'gene3_transcript_id'  : r.get('tx3', '.') or '.',
                'cds_left_range'       : '.',
                'cds_right_range'      : '.',
                'Fusion_effect'        : _frame(r.get('Fusion_effect', '')),
                'max_split_cnt'        : _to_int_str(dv),
                'SVTYPE'               : 'BND',
                'DV'                   : _to_int_str(dv),
                'RV'                   : (_to_int_str(rv) if rv != '.' else '.'),
                'FFPM'                 : '.',
                'somatic_flags'        : '.',
            })
    else:
        for _, r in df_raw.iterrows():
            l_chr, l_pos, l_strand = _bp_triplet(r.get('LeftBreakpoint', '.'))
            r_chr, r_pos, r_strand = _bp_triplet(r.get('RightBreakpoint', '.'))
            gene5 = _gene_symbol(r.get('LeftGene', '.'))
            gene3 = _gene_symbol(r.get('RightGene', '.'))
            rows.append({
                'gene5_chr'            : l_chr,
                'gene3_chr'            : r_chr,
                'gene5_breakpoint'     : l_pos,
                'gene3_breakpoint'     : r_pos,
                'gene5_strand'         : l_strand,
                'gene3_strand'         : r_strand,
                'gene5_renamed_symbol' : gene5,
                'gene3_renamed_symbol' : gene3,
                'gene5_tool_annotation': 'exon',
                'gene3_tool_annotation': 'exon',
                'gene5_transcript_id'  : r.get('CDS_LEFT_ID', '.') or '.',
                'gene3_transcript_id'  : r.get('CDS_RIGHT_ID', '.') or '.',
                'cds_left_range'       : (str(r.get('CDS_LEFT_RANGE', '.')).strip()  or '.'),
                'cds_right_range'      : (str(r.get('CDS_RIGHT_RANGE', '.')).strip() or '.'),
                'Fusion_effect'        : _frame(r.get('PROT_FUSION_TYPE', '.')),
                'max_split_cnt'        : str(r.get('JunctionReadCount', '.')),
                'SVTYPE'               : 'BND',
                'DV'                   : str(r.get('JunctionReadCount', '.')),
                'RV'                   : str(r.get('SpanningFragCount', '.')),
                'FFPM'                 : '.',
                'somatic_flags'        : '.',
            })
    return pd.DataFrame(rows, dtype=str)

def parse_txt_to_df(txt_path: str) -> pd.DataFrame:
    df_raw = pd.read_csv(txt_path, sep='\t', dtype=str, low_memory=False)
    first_col = df_raw.columns[0]
    def _gene_symbol(val: str) -> str:
        if pd.isna(val):
            return '.'
        v = str(val)
        return v.split('^', 1)[0] if '^' in v else v
    def _strand_from_flag(val: str) -> str:
        if pd.isna(val):
            return '+'
        v = str(val).strip()
        if v in ['+','-']:
            return v
        if v == '1':
            return '-'
        if v == '0':
            return '+'
        return '+'
    def _fusion_effect(val: str) -> str:
        if val is None or pd.isna(val):
            return 'unclear'
        v = str(val).strip().lower()
        if 'antisense' in v:
            return 'antisense'
        if 'in-frame' in v or 'inframe' in v or 'in frame' in v:
            return 'in-frame'
        if 'frame' in v and 'in' not in v:
            return 'out-of-frame'
        return 'unclear'
    def _to_int_str(x):
        if pd.isna(x) or str(x).strip() == '':
            return '0'
        try:
            return str(int(float(str(x))))
        except:
            return '0'
    rows = []
    if first_col == 'TumorId':
        for _, r in df_raw.iterrows():
            l_chr = r.get('Chr1', '.') or '.'
            r_chr = r.get('Chr2', '.') or '.'
            l_pos = r.get('Pos1', '0') or '0'
            r_pos = r.get('Pos2', '0') or '0'
            l_strand = _strand_from_flag(r.get('Str1', '+'))
            r_strand = _strand_from_flag(r.get('Str2', '+'))
            g5 = _gene_symbol(r.get('Gene1', '.'))
            g3 = _gene_symbol(r.get('Gene2', '.'))
            tx5 = r.get('Transcript1', '.') or '.'
            tx3 = r.get('Transcript2', '.') or '.'
            fusion_raw = r.get('Fusion', None)
            fusion_effect = _fusion_effect(fusion_raw)
            fusion_effect = 'unclear'
            support = _to_int_str(r.get('TotalReadSupport', '0'))
            rows.append({
                'gene5_chr'            : l_chr,
                'gene3_chr'            : r_chr,
                'gene5_breakpoint'     : l_pos,
                'gene3_breakpoint'     : r_pos,
                'gene5_strand'         : l_strand,
                'gene3_strand'         : r_strand,
                'gene5_renamed_symbol' : g5,
                'gene3_renamed_symbol' : g3,
                'gene5_tool_annotation': 'exon',
                'gene3_tool_annotation': 'exon',
                'gene5_transcript_id'  : tx5.split('.')[0],
                'gene3_transcript_id'  : tx3.split('.')[0],
                'cds_left_range'       : '.',
                'cds_right_range'      : '.',
                'Fusion_effect'        : fusion_effect,
                'max_split_cnt'        : support,
                'SVTYPE'               : 'BND',
                'DV'                   : support,
                'RV'                   : '0',
                'FFPM'                 : '.',
                'somatic_flags'        : '.',
            })
        return pd.DataFrame(rows, dtype=str)
    else:
        for _, r in df_raw.iterrows():
            l_chr = r.get('chr1', '.') or '.'
            r_chr = r.get('chr2', '.') or '.'
            l_pos = r.get('pos1', '0') or '0'
            r_pos = r.get('pos2', '0') or '0'
            l_strand = _strand_from_flag(r.get('str1', '+'))
            r_strand = _strand_from_flag(r.get('str2', '+'))
            g5 = _gene_symbol(r.get('gene1', '.'))
            g3 = _gene_symbol(r.get('gene2', '.'))
            tx5 = r.get('transcript1', '.') or '.'
            tx3 = r.get('transcript2', '.') or '.'
            fusion_raw = r.get('Fusion', None)
            fusion_effect = _fusion_effect(fusion_raw)
            rows.append({
                'gene5_chr'            : l_chr,
                'gene3_chr'            : r_chr,
                'gene5_breakpoint'     : l_pos,
                'gene3_breakpoint'     : r_pos,
                'gene5_strand'         : l_strand,
                'gene3_strand'         : r_strand,
                'gene5_renamed_symbol' : g5,
                'gene3_renamed_symbol' : g3,
                'gene5_tool_annotation': 'exon',
                'gene3_tool_annotation': 'exon',
                'gene5_transcript_id'  : tx5,
                'gene3_transcript_id'  : tx3,
                'cds_left_range'       : '.',
                'cds_right_range'      : '.',
                'Fusion_effect'        : fusion_effect,
                'max_split_cnt'        : '0',
                'SVTYPE'               : 'BND',
                'DV'                   : '0',
                'RV'                   : '0',
                'FFPM'                 : '.',
                'somatic_flags'        : '.',
            })
        return pd.DataFrame(rows, dtype=str)

def parse_cff_to_df(cff_path: str) -> pd.DataFrame:
    try:
        df_raw = pd.read_csv(cff_path, sep='\t', dtype=str, low_memory=False)
    except Exception as e:
        print(f"[ERROR] Failed to read CFF file {cff_path}: {e}", file=sys.stderr)
        return pd.DataFrame()
    required_columns = [
        'gene5_chr','gene5_breakpoint','gene5_strand',
        'gene3_chr','gene3_breakpoint','gene3_strand',
        'gene5_renamed_symbol','gene3_renamed_symbol',
        'is_inframe'
    ]
    missing_cols = [col for col in required_columns if col not in df_raw.columns]
    if missing_cols:
        print(f"[ERROR] Missing required columns in CFF file: {missing_cols}", file=sys.stderr)
        return pd.DataFrame()
    def _norm_effect(val: str) -> str:
        if val is None or pd.isna(val): return 'unclear'
        v = str(val).strip().lower()
        if v in ('in-frame','inframe','frame-preserving') or 'in-frame' in v or 'inframe' in v:
            return 'in-frame'
        if v in ('out-of-frame','outframe') or ('frame' in v and 'in' not in v):
            return 'out-of-frame'
        if v in ('unclear','unknown','na','.'):
            return 'unclear'
        return 'unclear'
    def _infer_effect_from_tool(a: str, b: str) -> str:
        a_eff = _norm_effect(a)
        b_eff = _norm_effect(b)
        if a_eff == 'in-frame' and b_eff == 'in-frame':
            return 'in-frame'
        if a_eff == 'out-of-frame' or b_eff == 'out-of-frame':
            return 'out-of-frame'
        if (a_eff == 'in-frame' and b_eff == 'unclear') or (b_eff == 'in-frame' and a_eff == 'unclear'):
            return 'in-frame'
        return 'unclear'
    def calculate_cds_range_from_breakpoint(gtf_df, contig, gene_symbol, transcript_id, breakpoint, strand, direction):
        if gtf_df is None or gtf_df.empty:
            return '.'
        gene_cds = gtf_df[
            (gtf_df['type'] == 'CDS') &
            (gtf_df['contig'] == contig) &
            (gtf_df['geneName'] == gene_symbol)
        ].copy()
        if (not transcript_id) or transcript_id == '.':
            if gene_cds.empty:
                return '.'
            transcript_id = gene_cds['transcript'].iloc[0]
        txid_core = transcript_id.split('.')[0]
        tx_cds = gene_cds[gene_cds['transcript'].str.contains(txid_core, na=False)].copy()
        if tx_cds.empty:
            return '.'
        tx_cds['start'] = pd.to_numeric(tx_cds['start'], errors='coerce').astype('Int64')
        tx_cds['end']   = pd.to_numeric(tx_cds['end'],   errors='coerce').astype('Int64')
        tx_cds = tx_cds.dropna(subset=['start','end']).astype({'start':int,'end':int}).sort_values('start')
        retained_length = 0
        for _, cds in tx_cds.iterrows():
            cds_start, cds_end = cds['start'], cds['end']
            cds_len = max(0, cds_end - cds_start)
            if strand == '+':
                if direction == 'downstream' and cds_start <= breakpoint:
                    retained_length += min(cds_len, max(0, breakpoint - cds_start + 1)) if cds_end > breakpoint else cds_len
                elif direction == 'upstream' and cds_end >= breakpoint:
                    if cds_start <= breakpoint:
                        retained_length += max(0, cds_end - breakpoint)
                    else:
                        retained_length += cds_len
            else:
                if direction == 'upstream' and cds_end >= breakpoint:
                    retained_length += min(cds_len, max(0, cds_end - breakpoint)) if cds_start < breakpoint else cds_len
                elif direction == 'downstream' and cds_start <= breakpoint:
                    if cds_end >= breakpoint:
                        retained_length += max(0, breakpoint - cds_start + 1)
                    else:
                        retained_length += cds_len
        return f"1-{retained_length}" if retained_length > 0 else '.'
    try:
        annotation = annotation[annotation['type'].isin(['CDS'])].copy()
        annotation['contig'] = annotation['contig'].apply(lambda x: re.sub(r"^(chr|CHR)", "", str(x)))
        annotation['geneName'] = annotation['attributes'].apply(
            lambda x: re.search(r'gene_name "([^"]+)"', str(x)).group(1) if re.search(r'gene_name "([^"]+)"', str(x)) else ''
        )
        annotation['transcript'] = annotation['attributes'].apply(
            lambda x: re.search(r'transcript_id "([^"]+)"', str(x)).group(1) if re.search(r'transcript_id "([^"]+)"', str(x)) else ''
        )
    except Exception:
        annotation = None
    rows = []
    for _, r in df_raw.iterrows():
        gene5_chr = str(r['gene5_chr']).replace('chr', '')
        gene3_chr = str(r['gene3_chr']).replace('chr', '')
        gene5_bp  = int(str(r['gene5_breakpoint']))
        gene3_bp  = int(str(r['gene3_breakpoint']))
        g5s = str(r['gene5_strand'])
        g3s = str(r['gene3_strand'])
        gene5_strand = g5s if g5s in ['+','-'] else '+'
        gene3_strand = g3s if g3s in ['+','-'] else '+'
        gene5_symbol = str(r['gene5_renamed_symbol'])
        gene3_symbol = str(r['gene3_renamed_symbol'])
        g5_tx_raw = r.get('gene5_transcript_id', '.')
        g3_tx_raw = r.get('gene3_transcript_id', '.')
        gene5_transcript = '.' if (g5_tx_raw is None or pd.isna(g5_tx_raw) or str(g5_tx_raw).strip()=='.') else str(g5_tx_raw)
        gene3_transcript = '.' if (g3_tx_raw is None or pd.isna(g3_tx_raw) or str(g3_tx_raw).strip()=='.') else str(g3_tx_raw)
        gene5_direction = 'downstream' if gene5_strand == '+' else 'upstream'
        gene3_direction = 'upstream' if gene3_strand == '-' else 'downstream'
        if annotation is not None:
            cds_left_range  = calculate_cds_range_from_breakpoint(annotation, gene5_chr, gene5_symbol, gene5_transcript, gene5_bp, gene5_strand, gene5_direction)
            cds_right_range = calculate_cds_range_from_breakpoint(annotation, gene3_chr, gene3_symbol, gene3_transcript, gene3_bp, gene3_strand, gene3_direction)
        else:
            cds_left_range = '.'
            cds_right_range = '.'
        fusion_effect_raw = 'in-frame' if r.get('is_inframe', '') else None
        fusion_effect_norm = _norm_effect(fusion_effect_raw)
        if fusion_effect_norm == 'unclear':
            fusion_effect_norm = _infer_effect_from_tool(
                str(r.get('gene5_tool_annotation','')),
                str(r.get('gene3_tool_annotation',''))
            )
        max_split_cnt = (str(r.get('max_split_cnt', '.')).strip()
                         if pd.notna(r.get('max_split_cnt')) else '.')
        row = {
            'gene5_chr'            : gene5_chr,
            'gene3_chr'            : gene3_chr,
            'gene5_breakpoint'     : str(gene5_bp),
            'gene3_breakpoint'     : str(gene3_bp),
            'gene5_strand'         : gene5_strand,
            'gene3_strand'         : gene3_strand,
            'gene5_renamed_symbol' : gene5_symbol,
            'gene3_renamed_symbol' : gene3_symbol,
            'gene5_tool_annotation': str(r.get('gene5_tool_annotation', 'exon')),
            'gene3_tool_annotation': str(r.get('gene3_tool_annotation', 'exon')),
            'gene5_transcript_id'  : gene5_transcript,
            'gene3_transcript_id'  : gene3_transcript,
            'cds_left_range'       : cds_left_range,
            'cds_right_range'      : cds_right_range,
            'Fusion_effect'        : fusion_effect_norm,
            'max_split_cnt'        : max_split_cnt,
            'SVTYPE'               : 'BND',
            'DV'                   : max_split_cnt,
            'RV'                   : (str(r.get('max_span_cnt', '.')).strip()
                                      if pd.notna(r.get('max_span_cnt')) else '.'),
            'FFPM'                 : '.',
            'somatic_flags'        : '.',
        }
        rows.append(row)
    return pd.DataFrame(rows, dtype=str)

def normalize_fusions(fusions):
    fusions                     = fusions.copy()
    fusions['contig1']          = fusions['gene5_chr'].apply(remove_chr)
    fusions['contig2']          = fusions['gene3_chr'].apply(remove_chr)
    fusions['Breakpoint1']      = fusions['gene5_breakpoint'].astype(int)
    fusions['Breakpoint2']      = fusions['gene3_breakpoint'].astype(int)
    fusions['strand1']          = fusions['gene5_strand']
    fusions['strand2']          = fusions['gene3_strand']
    fusions['direction1']       = fusions['strand1'].map({'+': 'downstream', '-': 'upstream'})
    fusions['direction2']       = fusions['strand2'].map({'+': 'upstream', '-': 'downstream'})
    fusions['gene_id1']         = fusions['gene5_renamed_symbol']
    fusions['gene_id2']         = fusions['gene3_renamed_symbol']
    fusions['FUSTYPE']          = '.'
    fusions['site1']            = fusions['gene5_tool_annotation'] if 'gene5_tool_annotation' in fusions.columns else 'exon'
    fusions['site2']            = fusions['gene3_tool_annotation'] if 'gene3_tool_annotation' in fusions.columns else 'exon'
    fusions['transcript_id1']   = fusions['gene5_transcript_id'].str.replace(r'\.\d+$', '', regex=True)
    fusions['transcript_id2']   = fusions['gene3_transcript_id'].str.replace(r'\.\d+$', '', regex=True)
    fusions['reading_frame']    = fusions['Fusion_effect'].apply(parse_fusion_frame) if 'Fusion_effect' in fusions.columns else 'unclear'
    fusions['split_reads']      = fusions['max_split_cnt'] if 'max_split_cnt' in fusions.columns else '.'
    fusions['display_contig1']  = fusions['gene5_chr']
    fusions['display_contig2']  = fusions['gene3_chr']
    fusions['confidence']       = 'high'
    def classify_fusion_type(row):
        if row['FUSTYPE']   == 'INV':
            return 'Inversion'
        elif row['FUSTYPE'] == 'DUP':
            return 'Duplication'
        elif row['FUSTYPE'] == 'TRA':
            return 'Translocation'
        elif row['FUSTYPE'] == 'DEL':
            return 'Deletion'
        elif row['contig1'] != row['contig2']:
            return 'Translocation'
        else:
            s1, s2 = row['strand1'], row['strand2']
            b1, b2 = int(row['Breakpoint1']), int(row['Breakpoint2'])
            if s1 == s2:
                if s1 == '+':
                    return 'Deletion' if b1 < b2 else 'Duplication'
                else:
                    return 'Deletion' if b1 > b2 else 'Duplication'
            else:
                return 'Inversion'
    fusions['type'] = fusions.apply(classify_fusion_type, axis=1)
    return fusions

def parse_gtf_attribute(attribute, attr_str):
    match = re.search(rf'{attribute} "([^"]+)"', attr_str)
    if match:
        return match.group(1)
    return ''

def draw_ideogram(ax, cytobands, contig, breakpoint, gene, xpos, width):
    pos = breakpoint
    bands = cytobands[cytobands['contig'] == contig]
    if bands.empty:
        return
    chrom_start = bands['start'].astype(int).min()
    chrom_end = bands['end'].astype(int).max()
    def norm(x):
        return xpos + width * (x - chrom_start) / (chrom_end - chrom_start)
    body_height = 0.1
    acen_bands = bands[bands['giemsa'] == 'acen']
    if not acen_bands.empty:
        acen_start_bp = acen_bands['start'].astype(int).min()
        acen_end_bp = acen_bands['end'].astype(int).max()
        acen_start = norm(acen_start_bp)
        acen_end = norm(acen_end_bp)
    else:
        acen_start = acen_end = (norm(chrom_start) + norm(chrom_end)) / 2
        acen_start_bp = acen_end_bp = (chrom_start + chrom_end) // 2
    total_width = norm(chrom_end) - norm(chrom_start)
    grad_res_x = 300
    grad_res_y = 75
    grad_y = np.linspace(0, 1, grad_res_y)
    grad_y = np.clip(1.0 - 1.2 * np.abs(grad_y - 0.75) * 2, 0.0, 1.0)
    grad = np.tile(grad_y[:, None], (1, grad_res_x))
    ax.imshow(
        grad,
        extent=[norm(chrom_start), norm(chrom_end), 0, body_height],
        aspect='auto', cmap='Greys', alpha=1.0, zorder=0
    )
    for _, band in bands.iterrows():
        left = norm(int(band['start']))
        right = norm(int(band['end']))
        if band['giemsa'] == 'acen':
            continue
        color = (
            '#FFFFFF' if band['giemsa'] == 'gneg' else
            '#000000' if band['giemsa'] == 'gpos100' else
            '#888888' if 'gpos' in band['giemsa'] else
            '#CCCCCC'
        )
        ax.add_patch(plt.Rectangle(
            (left, 0), right - left, body_height,
            color=color, ec='none', lw=0, zorder=1, alpha=0.7
        ))
    ax.add_patch(FancyBboxPatch(
        (norm(chrom_start) - 0.035, 0), total_width + 0.055, body_height,
        boxstyle=f"round,pad=0,rounding_size={body_height/3}",
        ec='#000000', fc='none', lw=2, zorder=3
    ))
    if not acen_bands.empty:
        ax.add_patch(plt.Rectangle(
            (acen_start - 0.007, 0), acen_end - (acen_start - 0.014), body_height,
            color='#FFFFFF', ec='#FFFFFF', lw=3, zorder=4
        ))
        n_lines = 200
        left_x1 = acen_start - 0.01
        left_x2 = ((acen_start + acen_end) / 2) - 0.005
        left_y1 = 0
        left_y2 = body_height
        for i in range(n_lines):
            y = left_y1 + (left_y2 - left_y1) * i / (n_lines - 1)
            y_norm = y / body_height
            intensity = np.clip(1.0 - 1.2 * abs(y_norm - 0.75) * 2, 0.0, 1.0)
            red_intensity = 0.4 + 0.6 * intensity
            color = (red_intensity, 0.1, 0.1)
            line_width = abs(left_x2 - left_x1) * (body_height / 2 - abs(y - body_height / 2)) / (body_height / 2)
            if line_width > 0:
                ax.plot([left_x1, left_x1 + line_width], [y, y],
                       color=color, linewidth=2.25, alpha=1, zorder=5)
        right_x1 = acen_end + 0.01
        right_x2 = (acen_start + acen_end) / 2
        right_y1 = 0
        right_y2 = body_height
        for i in range(n_lines):
            y = right_y1 + (right_y2 - right_y1) * i / (n_lines - 1)
            y_norm = y / body_height
            intensity = np.clip(1.0 - 1.2 * abs(y_norm - 0.75) * 2, 0.0, 1.0)
            red_intensity = 0.4 + 0.6 * intensity
            color = (red_intensity, 0.1, 0.1)
            line_width = abs(right_x2 - right_x1) * (body_height / 2 - abs(y - body_height / 2)) / (body_height / 2)
            if line_width > 0:
                ax.plot([right_x1 - line_width, right_x1], [y, y],
                       color=color, linewidth=2.25, alpha=1, zorder=5)
    bp_x = norm(breakpoint)
    ax.axvline(bp_x, 0.1, body_height * 3, color='#F400A1', lw=2, linestyle=':', dash_capstyle='round', zorder=6)
    ax.text(xpos + width / 2, body_height + 0.20, f'{gene}\n', ha='center', va='bottom', fontsize=22, fontweight='bold', fontstyle='italic')
    ax.text(xpos + width / 2, body_height + 0.18, f'Chromosome {contig}', ha='center', va='bottom', fontsize=12)
    ax.text(xpos + width / 2, body_height + 0.10, f'Position {pos}', ha='center', va='bottom', fontsize=10)
    band_name = bands[
        (bands['start'].astype(int) <= breakpoint) &
        (bands['end'].astype(int) >= breakpoint)
    ]['name']
    if not band_name.empty:
        ax.text(bp_x, body_height + 0.02, band_name.values[0],
                ha='center', va='bottom', fontsize=10)
    ax.set_ylim(-0.05, body_height + 0.18)
    ax.set_xlim(xpos - 0.07 * width, xpos + width + 0.07 * width)
    ax.axis('off')

def plot_chromosome_ideograms_to_axis(fusions_norm, cytobands, ax):
    ax.axis('off')
    if cytobands is None or fusions_norm.empty:
        ax.text(0.5, 0.5, "No cytobands data available", ha='center', va='center', fontsize=14)
        return
    for idx, row in fusions_norm.iterrows():
        left_x = -0.15
        right_x = 1.75
        width = 0.5
        draw_ideogram(ax, cytobands, row['contig1'], row['Breakpoint1'], row['gene_id1'], xpos=left_x, width=width)
        draw_ideogram(ax, cytobands, row['contig2'], row['Breakpoint2'], row['gene_id2'], xpos=right_x, width=width)
    ax.set_xlim(-0.2, 2.3)
    ax.set_ylim(-0.05, 0.5)

def get_per_exon_coverage(bam_path: str, contig: str, exons: pd.DataFrame) -> List[np.ndarray]:
    coverage_segments = []
    for _, row in exons.iterrows():
        start = int(row['start'])
        end = int(row['end'])
        exon_len = end - start
        if exon_len <= 0:
            coverage_segments.append(np.zeros(1))
            continue
        region = f"{contig}:{start}-{end}"
        try:
            result = subprocess.run(
                ["samtools", "depth", "-r", region, bam_path],
                capture_output=True, text=True, check=True
            )
            lines = result.stdout.strip().split('\n')
            coverage = np.zeros(exon_len)
            for line in lines:
                parts = line.split('\t')
                if len(parts) >= 3:
                    pos = int(parts[1])
                    cov = int(parts[2])
                    if start <= pos < end:
                        coverage[pos - start] = cov
            coverage_segments.append(coverage)
        except subprocess.CalledProcessError:
            coverage_segments.append(np.zeros(exon_len))
    return coverage_segments

def _try_provided_transcript(exons, transcript_id):
    if not transcript_id or transcript_id == '.':
        return None, None
    normalized_tx_id = str(transcript_id).strip()
    normalized_tx_id = re.sub(r'\.\d+$', '', normalized_tx_id)
    if normalized_tx_id in exons['transcript'].values:
        tx_exons = exons[exons['transcript'] == normalized_tx_id].copy()
        if not tx_exons.empty:
            return tx_exons, normalized_tx_id
    return None, None

def _try_canonical_transcript(exons):
    appris_levels = ['level 1', 'level 2', 'level 3', 'level 4']
    appris_types = [
        'appris_principal',
        'appris_candidate',
        'appris_principal_1',
        'appris_principal_2',
        'appris_principal_3',
        'appris_principal_4',
        'appris_principal_5'
    ]
    attrs_mask = exons['attributes'].fillna('').astype(str)
    for level in appris_levels:
        for appris in appris_types:
            level_mask = attrs_mask.str.contains(level, na=False, case=False)
            appris_mask = attrs_mask.str.contains(appris, na=False, case=False)
            level_exons = exons[level_mask & appris_mask].copy()
            if not level_exons.empty:
                tx_ids = level_exons['transcript'].unique()
                if len(tx_ids) > 0:
                    selected_tx = tx_ids[0]
                    tx_exons = exons[exons['transcript'] == selected_tx].copy()
                    return tx_exons, selected_tx
    return None, None

def _calculate_transcript_coverage(bam_path, contig, exons_df):
    try:
        coverage_segments = get_per_exon_coverage(bam_path, contig, exons_df)
        total_coverage = sum(np.sum(seg) for seg in coverage_segments)
        return total_coverage
    except Exception:
        return None

def _try_coverage_transcript(exons, contig, alignments_path):
    if not alignments_path or not os.path.isfile(alignments_path):
        return None, None
    unique_transcripts = exons['transcript'].dropna().unique()
    if len(unique_transcripts) == 0:
        return None, None
    best_tx = None
    best_coverage = -1
    for tx_id in unique_transcripts:
        tx_exons = exons[exons['transcript'] == tx_id].copy()
        if tx_exons.empty:
            continue
        coverage = _calculate_transcript_coverage(alignments_path, contig, tx_exons)
        if coverage is not None and coverage > best_coverage:
            best_coverage = coverage
            best_tx = tx_id
    if best_tx is not None:
        tx_exons = exons[exons['transcript'] == best_tx].copy()
        return tx_exons, best_tx
    return None, None

def _fallback_most_common_transcript(exons):
    if exons.empty:
        return None, None
    tx_counts = exons['transcript'].value_counts()
    if tx_counts.empty:
        return None, None
    best_tx = tx_counts.idxmax()
    tx_exons = exons[exons['transcript'] == best_tx].copy()
    return tx_exons, best_tx

def find_exons(annotation, contig, gene_id, transcript_id, transcript_selection='provided', alignments_path=None):
    ann = annotation.copy()
    def parse_gtf_attribute(key, attr):
        m = re.search(rf'{key}\s+"([^"]+)"', str(attr))
        return m.group(1) if m else ''
    def remove_chr(x):
        return str(x).replace('chr', '')
    if 'contig' in ann.columns:
        ann['contig_orig'] = ann['contig'].copy()
        ann['contig']   = ann['contig'].apply(remove_chr)
    if 'geneName' not in ann.columns or ann['geneName'].isna().any():
        ann['geneName'] = ann['attributes'].apply(lambda x: parse_gtf_attribute('gene_name', x))
    if 'gene_id' not in ann.columns or ann['gene_id'].isna().any():
        ann['gene_id']  = ann['attributes'].apply(lambda x: parse_gtf_attribute('gene_id', x))
    exons = pd.DataFrame()
    contig_clean = remove_chr(contig)
    contig_variants = [contig_clean, 'chr' + contig_clean] if not contig_clean.startswith('chr') else [contig_clean, contig_clean.replace('chr', '')]
    if transcript_id and transcript_id != '.':
        normalized_tx_id = str(transcript_id).strip()
        normalized_tx_id = re.sub(r'\.\d+$', '', normalized_tx_id)
        if 'transcript' not in ann.columns:
            ann['transcript'] = ann['attributes'].apply(lambda x: parse_gtf_attribute('transcript_id', x))
            ann['transcript'] = ann['transcript'].str.replace(r'\.\d+$', '', regex=True)
        tx_rows = ann[(ann['transcript'] == normalized_tx_id) & (ann['contig'] == contig_clean)].copy()
        if tx_rows.empty:
            for test_contig in contig_variants:
                if test_contig != contig_clean:
                    tx_rows = ann[(ann['transcript'] == normalized_tx_id) & (ann['contig'] == test_contig)].copy()
                    if not tx_rows.empty:
                        break
        if tx_rows.empty:
            attrs_str = ann['attributes'].fillna('').astype(str)
            tx_pattern = rf'transcript_id\s+"{re.escape(normalized_tx_id)}(?:\.[\d]+)?"'
            tx_match = attrs_str.str.contains(tx_pattern, na=False, regex=True)
            if not tx_match.any():
                tx_match = attrs_str.str.contains(re.escape(normalized_tx_id), na=False, regex=False)
            if tx_match.any():
                tx_rows = ann[tx_match & (ann['contig'] == contig_clean)].copy()
                if tx_rows.empty:
                    for test_contig in contig_variants:
                        if test_contig != contig_clean:
                            tx_rows = ann[tx_match & (ann['contig'] == test_contig)].copy()
                            if not tx_rows.empty:
                                break
        if not tx_rows.empty:
            tx_rows = tx_rows.copy()
            tx_rows['geneName'] = tx_rows['attributes'].apply(lambda x: parse_gtf_attribute('gene_name', x))
            actual_gene_names = tx_rows['geneName'].dropna().unique()
            if len(actual_gene_names) > 0:
                actual_gene_name = actual_gene_names[0]
                exons = ann[(ann['geneName'] == actual_gene_name) & (ann['contig'] == contig_clean)].copy()
                if exons.empty:
                    for test_contig in contig_variants:
                        if test_contig != contig_clean:
                            exons = ann[(ann['geneName'] == actual_gene_name) & (ann['contig'] == test_contig)].copy()
                            if not exons.empty:
                                break
    if exons.empty:
        exons = ann[(ann['geneName'] == gene_id) & (ann['contig'] == contig_clean)].copy()
        if exons.empty:
            for test_contig in contig_variants:
                if test_contig != contig_clean:
                    exons = ann[(ann['geneName'] == gene_id) & (ann['contig'] == test_contig)].copy()
                    if not exons.empty:
                        break
    if exons.empty:
        gene_id_lower = str(gene_id).lower()
        attrs_str = ann['attributes'].fillna('').astype(str).str.lower()
        gene_name_match = attrs_str.str.contains(rf'gene_name\s+"[^"]*{re.escape(gene_id_lower)}[^"]*"', na=False, regex=True)
        general_match = attrs_str.str.contains(re.escape(gene_id_lower), na=False, regex=False)
        gene_match = gene_name_match | general_match
        matched_rows = ann[gene_match & (ann['contig'] == contig_clean)].copy()
        if matched_rows.empty:
            for test_contig in contig_variants:
                if test_contig != contig_clean:
                    matched_rows = ann[gene_match & (ann['contig'] == test_contig)].copy()
                    if not matched_rows.empty:
                        break
        if not matched_rows.empty:
            matched_rows['geneName'] = matched_rows['attributes'].apply(lambda x: parse_gtf_attribute('gene_name', x))
            actual_gene_names = matched_rows['geneName'].dropna().unique()
            if len(actual_gene_names) > 0:
                actual_gene_name = actual_gene_names[0]
                exons = ann[(ann['geneName'] == actual_gene_name) & (ann['contig'] == contig_clean)].copy()
                if exons.empty:
                    for test_contig in contig_variants:
                        if test_contig != contig_clean:
                            exons = ann[(ann['geneName'] == actual_gene_name) & (ann['contig'] == test_contig)].copy()
                            if not exons.empty:
                                break
    if exons.empty and 'gene_id' in ann.columns:
        exons = ann[(ann['gene_id'] == gene_id) & (ann['contig'] == contig_clean)].copy()
        if exons.empty:
            for test_contig in contig_variants:
                if test_contig != contig_clean:
                    exons = ann[(ann['gene_id'] == gene_id) & (ann['contig'] == test_contig)].copy()
                    if not exons.empty:
                        break
    if exons.empty:
        return pd.DataFrame(), None
    if 'type' in exons.columns:
        exon_filtered = exons[exons['type'] == 'exon'].copy()
        if not exon_filtered.empty:
            exons = exon_filtered
        elif exons.empty:
            return pd.DataFrame(), None
    if 'transcript' not in exons.columns:
        def parse_gtf_attribute(key, attr):
            m = re.search(rf'{key}\s+"([^"]+)"', str(attr))
            return m.group(1) if m else ''
        exons['transcript'] = exons['attributes'].apply(lambda x: parse_gtf_attribute('transcript_id', x))
        exons['transcript'] = exons['transcript'].str.replace(r'\.\d+$', '', regex=True)
    exons['start'] = exons['start'].astype(int)
    exons['end'] = exons['end'].astype(int)
    exons = exons[exons['transcript'].notna() & (exons['transcript'] != '')].copy()
    if exons.empty:
        return pd.DataFrame(), None
    result = None, None
    if transcript_selection == 'provided':
        result = _try_provided_transcript(exons, transcript_id)
        if result[0] is None:
            result = _try_canonical_transcript(exons)
            if result[0] is None:
                result = _try_coverage_transcript(exons, contig_clean, alignments_path)
    elif transcript_selection == 'canonical':
        result = _try_canonical_transcript(exons)
        if result[0] is None:
            result = _try_coverage_transcript(exons, contig_clean, alignments_path)
    elif transcript_selection == 'coverage':
        result = _try_coverage_transcript(exons, contig_clean, alignments_path)
        if result[0] is None:
            result = _try_canonical_transcript(exons)
    if result[0] is None:
        result = _fallback_most_common_transcript(exons)
    if result[0] is None or result[0].empty:
        return pd.DataFrame(), None
    return result

def _merge_intervals_for_panelB(segs):
    if not segs:
        return []
    segs = sorted((int(s), int(e)) for s, e in segs if s is not None and e is not None and e > s)
    merged = [segs[0]]
    for s, e in segs[1:]:
        ms, me = merged[-1]
        if s <= me:
            merged[-1] = (ms, max(me, e))
        else:
            merged.append((s, e))
    return merged

def _subtract_intervals_for_panelB(whole, parts):
    S, E = int(whole[0]), int(whole[1])
    if S >= E:
        return []
    if not parts:
        return [(S, E)]
    outs, cur = [], S
    for s, e in parts:
        if e <= cur:
            continue
        if s > cur:
            outs.append((cur, min(s, E)))
        cur = max(cur, e)
        if cur >= E:
            break
    if cur < E:
        outs.append((cur, E))
    return outs

def _draw_greys_gradient(ax, x0, w, y0, h, alpha=1, zorder=1):
    import numpy as _np
    if w <= 0 or h <= 0:
        return
    res_x = max(int(w * 1000), 10)
    res_y = 150
    gy = _np.linspace(0, 1, res_y)
    gy = _np.clip(1.0 - 2.0 * _np.abs(gy - 0.8) * 1.5, 0.0, 1.0)
    grad = _np.tile(gy[:, None], (1, res_x))
    ax.imshow(grad, extent=[x0, x0 + w, y0, y0 + h],
              aspect='auto', cmap='Greys', alpha=alpha, zorder=zorder)

def draw_exon_with_cds_utr(ax, start_x, width_x, g_start, g_end, cds_intervals,
                                  color, y_bottom, y_top, utr_center_frac=0.5,
                                  utr_alpha=0.55, cds_alpha=0.55, zorder_base=1):
    exon_len_g = max(1, int(g_end) - int(g_start))
    def gx(pos):
        return start_x + (int(pos) - int(g_start)) * width_x / exon_len_g
    cds_list = []
    for cs, ce in (cds_intervals or []):
        s = max(int(g_start), int(cs))
        e = min(int(g_end), int(ce))
        if s < e:
            cds_list.append((s, e))
    cds_merged = _merge_intervals_for_panelB(cds_list)
    utr_segments = _subtract_intervals_for_panelB((int(g_start), int(g_end)), cds_merged)
    exon_h = y_top - y_bottom
    utr_h = max(0.0, exon_h * utr_center_frac)
    utr_y0 = y_bottom + (exon_h - utr_h) / 2.0
    for us, ue in utr_segments:
        ux0 = gx(us)
        uw = gx(ue) - ux0
        if uw <= 0:
            continue
        _draw_greys_gradient(ax, ux0, uw, utr_y0, utr_h, alpha=1, zorder=zorder_base)
        ax.add_patch(plt.Rectangle((ux0, utr_y0), uw, utr_h, color=color, ec=None, alpha=utr_alpha, zorder=zorder_base+1))
    for cs, ce in cds_merged:
        cx0 = gx(cs)
        cw = gx(ce) - cx0
        if cw <= 0:
            continue
        _draw_greys_gradient(ax, cx0, cw, y_bottom, exon_h, alpha=1, zorder=zorder_base+2)
        ax.add_patch(plt.Rectangle((cx0, y_bottom), cw, exon_h, color=color, ec=None, alpha=cds_alpha, zorder=zorder_base+3))

def plot_gene(ax, exons, cov_segments, gene_name, color, breakpoint, x_offset, width, transcript_id, cds_intervals=None):
    y_exon_bottom = 0.0
    y_exon_top = 0.15
    coords = []
    curr_x = 0
    total_exon_len = sum(int(row['end']) - int(row['start']) for _, row in exons.iterrows())
    exons = exons.copy()
    exons = exons.sort_values(by='start', ascending=True).reset_index(drop=True)
    for i, row in exons.iterrows():
        exon_width = int(row['end']) - int(row['start'])
        coords.append((curr_x, exon_width, int(row['start']), int(row['end'])))
        curr_x += exon_width
        if i != len(exons) - 1:
            curr_x += int(0.02 * total_exon_len)
    total_len2 = curr_x
    scale2 = width / total_len2 if total_len2 > 0 else 1.0
    coords_scaled = [(x_offset + start * scale2, width2 * scale2, g_start, g_end) for (start, width2, g_start, g_end) in coords]
    max_cov = max((np.max(c) for c in cov_segments), default=1)
    coverage_bottom = y_exon_top + 0.01
    coverage_top = coverage_bottom + 0.2
    ax.text(x_offset - 0.004, coverage_top, f'{int(max_cov)}', ha='right', va='bottom', fontsize=9, color=color)
    ax.text(x_offset - 0.003, (coverage_top + coverage_bottom)/2, 'Coverage',
            ha='right', va='center', fontsize=9, color=color, rotation=90)
    for (start, width2, g_start, g_end), cov in zip(coords_scaled, cov_segments):
        exon_len = len(cov)
        bin_width = width2 / exon_len if exon_len > 0 else width2
        for i in range(exon_len):
            x = start + i * bin_width
            h = (cov[i] / max_cov) * 0.2 if max_cov > 0 else 0
            ax.add_patch(plt.Rectangle((x, y_exon_top + 0.01), bin_width, h, color=color, alpha=1, zorder=1))
    for i, ((start, width2, g_start, g_end), (_, exon_row)) in enumerate(zip(coords_scaled, exons.iterrows())):
        exon_height = y_exon_top - y_exon_bottom
        draw_exon_with_cds_utr(
            ax,
            start_x=start,
            width_x=width2,
            g_start=g_start,
            g_end=g_end,
            cds_intervals=cds_intervals,
            color=color,
            y_bottom=y_exon_bottom,
            y_top=y_exon_top,
            utr_center_frac=0.5,
            utr_alpha=0.55,
            cds_alpha=0.55,
            zorder_base=1
        )
        if 'exon_number' in exon_row:
            exon_num = exon_row['exon_number']
        elif 'attributes' in exon_row:
            match = re.search(r'exon_number\s+"?(\d+)"?', str(exon_row['attributes']))
            exon_num = match.group(1) if match else str(i + 1)
        else:
            exon_num = str(i + 1)
        try:
            exon_num_int = int(exon_num)
        except Exception:
            exon_num_int = i + 1
        exon_center_x = start + width2 / 2
        exon_center_y = y_exon_bottom + exon_height / 2
        dy = 0.1 * exon_height
        y_label = exon_center_y + dy if (exon_num_int % 2 == 1) else exon_center_y - dy
        txt = ax.text(
            exon_center_x, y_label, exon_num,
            ha='center', va='center', fontsize=11, fontweight='bold',
            color='#FFFFFF', zorder=10
        )
        txt.set_path_effects([pe.withStroke(linewidth=0.25, foreground='#000000')])
    if breakpoint is not None:
        breakpoint_x = None
        for (start, width2, g_start, g_end) in coords_scaled:
            if g_start <= breakpoint <= g_end:
                relative_pos = (breakpoint - g_start) / (g_end - g_start)
                breakpoint_x = start + relative_pos * width2
                break
        if breakpoint_x is None:
            for i in range(len(coords_scaled) - 1):
                curr_exon = coords_scaled[i]
                next_exon = coords_scaled[i + 1]
                curr_g_end = curr_exon[3]
                next_g_start = next_exon[2]
                if curr_g_end < breakpoint < next_g_start:
                    curr_visual_end = curr_exon[0] + curr_exon[1]
                    next_visual_start = next_exon[0]
                    breakpoint_x = (curr_visual_end + next_visual_start) / 2
                    break
        if breakpoint_x is not None:
            ax.plot([breakpoint_x, breakpoint_x], [y_exon_bottom, coverage_top],
                   color='#F400A1', lw=1.5, linestyle=':', dash_capstyle='round', zorder=2)
    line, = ax.plot(
        [x_offset, x_offset],
        [coverage_bottom, coverage_top],
        lw=1,
        color='#2A2A2A',
        alpha=1.0,
        solid_capstyle='butt',
        zorder=10,
        clip_on=False,
    )
    line.set_path_effects([pe.Stroke(linewidth=1, foreground='#808080'), pe.Normal()])

def reorder_cov_to_match_plot_gene(exons: pd.DataFrame, cov_segs: List[np.ndarray]) -> List[np.ndarray]:
    if len(cov_segs) == 0 or exons.empty:
        return cov_segs
    keys_unsorted = [(int(r['start']), int(r['end'])) for _, r in exons.iterrows()]
    cov_map = {k: cov for k, cov in zip(keys_unsorted, cov_segs)}
    exons_sorted = exons.sort_values(by='start', ascending=True)
    cov_sorted = []
    for _, r in exons_sorted.iterrows():
        k = (int(r['start']), int(r['end']))
        cov_sorted.append(cov_map.get(k, np.zeros(int(r['end']) - int(r['start']))))
    return cov_sorted

def retained_from_cds_range(coding_exons, cds_intervals, strand, cds_range):
    return coding_exons.copy()

def determine_retained_from_breakpoint(exons_df, breakpoint, direction, strand):
    if exons_df.empty:
        return pd.DataFrame(columns=['start', 'end', 'strand'])
    exons = exons_df.copy()
    exons['start'] = exons['start'].astype(int)
    exons['end'] = exons['end'].astype(int)
    retained = []
    for _, exon in exons.iterrows():
        exon_start = int(exon['start'])
        exon_end = int(exon['end'])
        if direction == 'downstream':
            if exon_start <= breakpoint:
                retained_start = exon_start
                retained_end = min(exon_end, breakpoint)
                if retained_start < retained_end:
                    retained.append({
                        'start': retained_start,
                        'end': retained_end,
                        'strand': strand
                    })
        else:
            if exon_end >= breakpoint:
                retained_start = max(exon_start, breakpoint)
                retained_end = exon_end
                if retained_start < retained_end:
                    retained.append({
                        'start': retained_start,
                        'end': retained_end,
                        'strand': strand
                    })
    return pd.DataFrame(retained)

def extract_cds_intervals(gtf_df, contig, gene_id, transcript_id):
    import pandas as pd
    if gtf_df is None or gtf_df.empty or not transcript_id:
        return []
    feat_col = 'type' if 'type' in gtf_df.columns else ('feature' if 'feature' in gtf_df.columns else None)
    if feat_col is None or 'contig' not in gtf_df.columns or 'transcript' not in gtf_df.columns:
        return []
    tx = str(transcript_id).split('.', 1)[0]
    tcol = gtf_df['transcript'].astype(str).str.replace(r'\.\d+$', '', regex=True)
    mask = (gtf_df[feat_col] == 'CDS') & (gtf_df['contig'] == contig) & (tcol == tx)
    if gene_id:
        g = str(gene_id).split('.', 1)[0]
        gmasks = []
        for col in ('gene_id','gene','geneName','gene_name','gene_symbol','geneId','geneid'):
            if col in gtf_df.columns:
                gmasks.append(gtf_df[col].astype(str).str.replace(r'\.\d+$','',regex=True) == g)
        if gmasks:
            mask = mask & pd.concat(gmasks, axis=1).any(axis=1)
    cds = gtf_df.loc[mask, ['start','end']].copy()
    if cds.empty:
        return []
    cds['start'] = pd.to_numeric(cds['start'], errors='coerce')
    cds['end']   = pd.to_numeric(cds['end'],   errors='coerce')
    cds = cds.dropna(subset=['start','end']).astype({'start': int, 'end': int}).sort_values('start')
    return [(int(s), int(e)) for s, e in cds[['start','end']].itertuples(index=False)]

def get_unified_transcript_data(annotation, contig, gene_id, direction, breakPoint, transcript_id, transcript_selection='provided', alignments_path=None):
    original_tx_id = transcript_id
    cds = pd.DataFrame()
    for attempt in range(2):
        exons, selected_tx = find_exons(annotation, contig, gene_id, transcript_id, transcript_selection, alignments_path)
        if exons is None or exons.empty:
            exons = pd.DataFrame(columns=['contig','start','end','strand','attributes','type'])
        if not exons.empty and ('start' not in exons.columns or 'end' not in exons.columns):
            raise RuntimeError(f"[ERROR] Exons missing required columns (start/end) for {gene_id}")
        if exons.empty or selected_tx is None or pd.isna(selected_tx) or str(selected_tx).strip() == '':
            if attempt == 0 and original_tx_id and original_tx_id != '.':
                transcript_id = None
                continue
            gene_found = not annotation[annotation['geneName'] == gene_id].empty if 'geneName' in annotation.columns else False
            contig_found = not annotation[annotation['contig'] == contig].empty if 'contig' in annotation.columns else False
            error_msg = f"[ERROR] No exon rows for the selected transcript(s). Gene: {gene_id}, Contig: {contig}"
            if not gene_found:
                error_msg += f" (Gene '{gene_id}' not found in annotation)"
            elif not contig_found:
                error_msg += f" (Contig '{contig}' not found in annotation)"
            raise RuntimeError(error_msg)
        exons['start'] = exons['start'].astype(int)
        exons['end']   = exons['end'].astype(int)
        if 'attributes' not in exons.columns:
            exons['attributes'] = ''
        if 'exon_number' not in exons.columns:
            exons['exon_number'] = exons['attributes'].apply(
                lambda x: (re.search(r'exon_number\s+"?(\d+)"?', str(x)).group(1)
                        if re.search(r'exon_number\s+"?(\d+)"?', str(x)) else ''))
        def parse_gtf_attribute(key, attr):
            m = re.search(rf'{key}\s+"([^"]+)"', str(attr))
            return m.group(1) if m else ''
        def remove_chr(x):
            return str(x).replace('chr', '')
        gtf = annotation.copy()
        gtf['contig']     = gtf['contig'].apply(remove_chr)
        gtf['geneName']   = gtf['attributes'].apply(lambda x: parse_gtf_attribute('gene_name', x))
        gtf['transcript'] = gtf['attributes'].apply(lambda x: parse_gtf_attribute('transcript_id', x))
        gtf['transcript'] = gtf['transcript'].str.replace(r'\.\d+$', '', regex=True)
        def _cds_for_tx(gtf_df, gene, tx):
            if tx is None or pd.isna(tx) or tx == '':
                return pd.DataFrame()
            try:
                cds = gtf_df[(gtf_df['type']=='CDS') &
                             (gtf_df['geneName']==gene) &
                             (gtf_df['transcript']==tx)]
                return cds
            except Exception:
                return pd.DataFrame()
        cds = _cds_for_tx(gtf, gene_id, selected_tx)
        if cds.empty and not exons.empty:
            break
        if cds.empty:
            if attempt == 0 and original_tx_id and original_tx_id != '.':
                transcript_id = None
                continue
        break
    def _cds_intervals(df):
        if df is None or df.empty:
            return []
        return [(int(s), int(e)) for s, e in zip(df['start'], df['end'])]
    has_cds = cds is not None and not cds.empty
    if has_cds:
        cds_intervals = _cds_intervals(cds)
    else:
        cds_intervals = []
    def _intersect_exons(exons_df, intervals):
        if not intervals or exons_df.empty:
            if exons_df.empty:
                return pd.DataFrame(columns=['start','end','strand'])
            return exons_df[['start','end','strand']].copy()
        rows = []
        for _, ex in exons_df.iterrows():
            es, ee = int(ex['start']), int(ex['end'])
            for (cs, ce) in intervals:
                s = max(es, int(cs)); e = min(ee, int(ce))
                if s < e:
                    rows.append({'start': s, 'end': e, 'strand': ex.get('strand','+')})
        return pd.DataFrame(rows, columns=['start','end','strand'])
    coding_exons = _intersect_exons(exons, cds_intervals)
    strand_val = '+'
    if not exons.empty and 'strand' in exons.columns:
        strand_val = exons['strand'].iloc[0] if not pd.isna(exons['strand'].iloc[0]) else '+'
    retained_exons = determine_retained_from_breakpoint(coding_exons, breakPoint, direction, strand_val)
    def _annotate_retained_with_exon_numbers(retained, exons_full):
        if retained is None or retained.empty or exons_full is None or exons_full.empty:
            return retained if retained is not None else pd.DataFrame(columns=['start','end'])
        ex_full = exons_full.copy().sort_values('start', ascending=True).reset_index(drop=True)
        if 'exon_number' not in ex_full.columns:
            if 'attributes' in ex_full.columns:
                ex_full['exon_number'] = ex_full['attributes'].apply(
                    lambda x: (re.search(r'exon_number\s+"?(\d+)"?', str(x)).group(1)
                            if re.search(r'exon_number\s+"?(\d+)"?', str(x)) else ''))
            else:
                ex_full['exon_number'] = ''
        ex_full['exon_number'] = ex_full['exon_number'].replace('', np.nan)
        out = []
        for _, r in retained.iterrows():
            rs, re = int(r['start']), int(r['end'])
            hit = ex_full[(ex_full['start'] <= rs) & (ex_full['end'] >= re)]
            exon_num = hit['exon_number'].iloc[0] if not hit.empty else np.nan
            rr = r.copy()
            rr['exon_number'] = '' if pd.isna(exon_num) else str(exon_num)
            out.append(rr)
        return pd.DataFrame(out) if out else retained
    retained_exons = _annotate_retained_with_exon_numbers(retained_exons, exons)
    if selected_tx is None or pd.isna(selected_tx) or str(selected_tx).strip() == '':
        raise RuntimeError(f"[ERROR] No transcript found for {gene_id} on {contig}.")
    if exons is None or exons.empty:
        raise RuntimeError(f"[ERROR] No exons found for {gene_id} transcript {selected_tx}.")
    if coding_exons is None:
        coding_exons = pd.DataFrame(columns=['start','end','strand'])
    if cds_intervals is None:
        cds_intervals = []
    if retained_exons is None:
        retained_exons = pd.DataFrame(columns=['start','end'])
    selected_tx = str(selected_tx).strip()
    return selected_tx, exons, coding_exons, cds_intervals, retained_exons

def plot_coverage_gene_structure_to_axis(fusions_norm, args, ax):
    ax.axis('off')
    try:
        annotation = pd.read_csv(args.annotation, sep='\t', comment='#', header=None,
                                names=['contig','src','type','start','end','score','strand','frame','attributes'], dtype=str)
    except Exception as e:
        ax.text(0.5, 0.5, f"Error reading annotation: {e}", ha='center', va='center', fontsize=14)
        return
    annotation = fix_missing_type_column(annotation)
    annotation['contig'] = annotation['contig'].apply(remove_chr)
    annotation['geneName'] = annotation['attributes'].apply(lambda x: parse_gtf_attribute('gene_name', x))
    annotation['transcript'] = annotation['attributes'].apply(lambda x: parse_gtf_attribute('transcript_id', x))
    annotation['transcript'] = annotation['transcript'].str.replace(r'\.\d+$', '', regex=True)
    for idx, row in fusions_norm.iterrows():
        tx1, exons1, coding_exons1, cds_intervals1, retained_exons1 = get_unified_transcript_data(
            annotation, row['contig1'], row['gene_id1'], row['direction1'],
            row['Breakpoint1'], row['transcript_id1'], args.transcriptSelection, args.alignments)
        tx2, exons2, coding_exons2, cds_intervals2, retained_exons2 = get_unified_transcript_data(
            annotation, row['contig2'], row['gene_id2'], row['direction2'],
            row['Breakpoint2'], row['transcript_id2'], args.transcriptSelection, args.alignments)
        if exons1.empty or exons2.empty:
            continue
        exon_len1 = sum(int(r['end']) - int(r['start']) for _, r in exons1.iterrows())
        exon_len2 = sum(int(r['end']) - int(r['start']) for _, r in exons2.iterrows())
        fixed_intron_width = 0.01
        g1_introns = (len(exons1) - 1) * fixed_intron_width
        g2_introns = (len(exons2) - 1) * fixed_intron_width
        g1_span = exon_len1 + g1_introns
        g2_span = exon_len2 + g2_introns
        total_span = g1_span + g2_span
        gap_frac = 0.025
        gap_span = gap_frac * total_span
        total_panel_span = g1_span + gap_span + g2_span
        g1_width = g1_span / total_panel_span
        g2_width = g2_span / total_panel_span
        gap_width = gap_span / total_panel_span
        if args.alignments:
            cov_segs1_raw = get_per_exon_coverage(args.alignments, row['contig1'], exons1)
            cov_segs2_raw = get_per_exon_coverage(args.alignments, row['contig2'], exons2)
            cov_segs1 = reorder_cov_to_match_plot_gene(exons1, cov_segs1_raw)
            cov_segs2 = reorder_cov_to_match_plot_gene(exons2, cov_segs2_raw)
        else:
            cov_segs1 = [np.zeros(int(r['end']) - int(r['start'])) for _, r in exons1.iterrows()]
            cov_segs2 = [np.zeros(int(r['end']) - int(r['start'])) for _, r in exons2.iterrows()]
        plot_gene(ax, exons1, cov_segs1, row['gene_id1'], args.color1, row['Breakpoint1'],
                  0, g1_width, tx1, cds_intervals=cds_intervals1)
        plot_gene(ax, exons2, cov_segs2, row['gene_id2'], args.color2, row['Breakpoint2'],
                  g1_width + gap_width, g2_width, tx2, cds_intervals=cds_intervals2)
        L_layout, L_map, _ = _build_linear_exon_layout(exons1, 0.0, g1_width)
        R_layout, R_map, _ = _build_linear_exon_layout(exons2, g1_width + gap_width, g2_width)
        y_underline = -0.01
        def _draw_retained_underlines(retained_df: pd.DataFrame, mapper, color: str):
            if retained_df.empty:
                return
            spans = []
            BODY_LW = 3
            HEAD_LEN_PT = 0.1
            HEAD_W_PT = 0.05
            for _, seg in retained_df.iterrows():
                try:
                    xs = mapper(int(seg['start']))
                    xe = mapper(int(seg['end']))
                except Exception:
                    continue
                if xs is None or xe is None or xs == xe:
                    continue
                spans.append((min(xs, xe), max(xs, xe)))
            if not spans:
                return
            x_left = min(s[0] for s in spans)
            x_right = max(s[1] for s in spans)
            strand = retained_df.iloc[0].get('strand', '+')
            if strand == '+':
                x0, x1 = x_left, x_right
            else:
                x0, x1 = x_right, x_left
            length = abs(x1 - x0)
            if length <= 1e-6:
                return
            ms = max(12, min(48, 240 * length))
            arrow = mpatches.FancyArrowPatch(
                (x0, y_underline - 0.005), (x1, y_underline - 0.005),
                arrowstyle=f'-|>,head_length={HEAD_LEN_PT},head_width={HEAD_W_PT}',
                mutation_scale=ms,
                linewidth=BODY_LW,
                color=color,
                zorder=6,
                shrinkA=0, shrinkB=0,
                capstyle='butt', joinstyle='miter'
            )
            ax.add_patch(arrow)
        try:
            strand1 = exons1['strand'].iloc[0] if 'strand' in exons1.columns and not exons1.empty else '+'
            strand2 = exons2['strand'].iloc[0] if 'strand' in exons2.columns and not exons2.empty else '+'
            retained_full1 = determine_retained_from_breakpoint(
                exons1[['start', 'end', 'strand']] if {'start','end','strand'}.issubset(exons1.columns) else exons1,
                row['Breakpoint1'], row['direction1'], strand1
            )
            retained_full2 = determine_retained_from_breakpoint(
                exons2[['start', 'end', 'strand']] if {'start','end','strand'}.issubset(exons2.columns) else exons2,
                row['Breakpoint2'], row['direction2'], strand2
            )
        except Exception:
            retained_full1 = retained_exons1
            retained_full2 = retained_exons2
        _draw_retained_underlines(retained_full1, L_map, args.color1)
        _draw_retained_underlines(retained_full2, R_map, args.color2)
        x_bp_left = L_map(int(row['Breakpoint1']))
        x_bp_right = R_map(int(row['Breakpoint2']))
        if x_bp_left is not None and x_bp_right is not None:
            try:
                bridge_support = collect_fusion_bridge_strength(
                    args.alignments,
                    row['contig1'], int(row['Breakpoint1']),
                    row['contig2'], int(row['Breakpoint2']),
                    args.sashimiWindow, args.minMapQ
                ) if args.alignments else 0
            except Exception:
                bridge_support = 0
            apex_y = 0.5
            midx = (x_bp_left + x_bp_right) / 2.0
            segs = 60
            brdy = 0.36
            ptsx, ptsy = [], []
            for t in np.linspace(0, 1, segs):
                bx = (1 - t) ** 2 * x_bp_left + 2 * (1 - t) * t * midx + t ** 2 * x_bp_right
                by = (1 - t) ** 2 * brdy + 2 * (1 - t) * t * apex_y + t ** 2 * brdy
                ptsx.append(bx)
                ptsy.append(by)
            points = np.column_stack([ptsx, ptsy])
            segments = np.stack([points[:-1], points[1:]], axis=1)
            t = np.linspace(0, 1, segs)
            support_scale = min(bridge_support, 30) / 30.0
            lw_apex = 1.0 + 3.0 * support_scale
            lw_ends = 1.0
            fade = np.abs(t - 0.5) / 0.5
            lw_profile = lw_apex - (lw_apex - lw_ends) * fade
            lc = LineCollection( segments, linewidths=lw_profile[:-1], color='#F400A1', capstyle='round', joinstyle='round', alpha=1, zorder=2 )
            ax.add_collection(lc)
            ax.text(midx, apex_y - 0.045, 'Junction', ha='center', va='center', fontsize=10, fontweight='bold', color='#F400A1', zorder=8)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.15, 0.5)

def fix_missing_type_column(annotation):
    if 'type' not in annotation.columns:
        if 'feature' in annotation.columns:
            annotation['type'] = annotation['feature']
        elif 'attributes' in annotation.columns:
            annotation['type'] = annotation['attributes'].apply(
                lambda x: 'exon' if 'exon' in str(x).lower() else 'CDS' if 'cds' in str(x).lower() else 'gene'
            )
        else:
            annotation['type'] = 'exon'
    return annotation

def remove_introns_from_domains(coding_exons, retained_domains):
    if coding_exons.empty or retained_domains.empty:
        return retained_domains
    ex = coding_exons[['start', 'end', 'strand']].copy()
    ex['start'] = pd.to_numeric(ex['start'], errors='coerce').astype('Int64')
    ex['end']   = pd.to_numeric(ex['end'],   errors='coerce').astype('Int64')
    ex = ex.dropna(subset=['start','end']).astype({'start': int, 'end': int})
    if ex.empty:
        return retained_domains
    strand = ex['strand'].iloc[0] if 'strand' in ex.columns else '+'
    ex = ex.sort_values('start', ascending=(strand == '+')).reset_index(drop=True)
    offsets = []
    running = 0
    for _, r in ex.iterrows():
        s, e = int(r['start']), int(r['end'])
        L = max(0, e - s)
        offsets.append((s, e, running))
        running += L
    if running == 0:
        return retained_domains
    def map_pos(pos):
        for s, e, off in offsets:
            if s <= pos <= e:
                return off + (pos - s) if strand == '+' else off + (e - pos)
        if pos < offsets[0][0]:
            return 0
        if pos > offsets[-1][1]:
            return running
        for i in range(len(offsets) - 1):
            s0, e0, off0 = offsets[i]
            s1, e1, off1 = offsets[i + 1]
            if e0 < pos < s1:
                return off0 + (e0 - s0)
        return 0
    rd = retained_domains.copy()
    if 'Start' not in rd.columns or 'End' not in rd.columns:
        if 'start' in rd.columns and 'End' not in rd.columns: rd['End'] = rd['start']
        if 'end' in rd.columns and 'Start' not in rd.columns: rd['Start'] = rd['end']
    rd['Start'] = pd.to_numeric(rd['Start'], errors='coerce').fillna(0).astype(int)
    rd['End']   = pd.to_numeric(rd['End'],   errors='coerce').fillna(0).astype(int)
    rd['Start'] = rd['Start'].map(map_pos)
    rd['End']   = rd['End'].map(map_pos)
    starts = rd['Start'].values
    ends   = rd['End'].values
    swapped = starts > ends
    if swapped.any():
        rd.loc[swapped, 'Start'], rd.loc[swapped, 'End'] = ends[swapped], starts[swapped]
    return rd

def merge_similar_domains(domains, merge_domains_overlapping_by):
    if domains.empty:
        return domains
    try:
        domains = domains.copy()
        if 'start' in domains.columns and 'Start' not in domains.columns:
            domains['Start'] = domains['start']
        if 'end' in domains.columns and 'End' not in domains.columns:
            domains['End'] = domains['end']
        if 'Start' not in domains.columns or 'End' not in domains.columns:
            print(f'[ERROR] Missing Start/End columns: {list(domains.columns)}')
            return domains
        domains['Start'] = pd.to_numeric(domains['Start'], errors='coerce')
        domains['End']   = pd.to_numeric(domains['End'],   errors='coerce')
        domains = domains.dropna(subset=['Start', 'End'])
        if domains.empty:
            return domains
        if 'proteinDomainID' in domains.columns:
            merged_by_id = []
            for uid, grp in domains.groupby('proteinDomainID', sort=False):
                if grp.empty:
                    continue
                row = grp.iloc[0].copy()
                row['Start'] = grp['Start'].min()
                row['End']   = grp['End'].max()
                merged_by_id.append(row)
            domains = pd.DataFrame(merged_by_id)
            if domains.empty:
                return domains
        domains['width'] = domains['End'] - domains['Start']
        domains = domains.sort_values('width', ascending=False)
        merged = pd.DataFrame(columns=domains.columns)
        for _, domain in domains.iterrows():
            overlap = False
            for _, m in merged.iterrows():
                try:
                    overlap_frac = (min(domain['End'], m['End']) - max(domain['Start'], m['Start'])) / (domain['End'] - domain['Start'] + 1)
                    if overlap_frac > merge_domains_overlapping_by:
                        overlap = True
                        break
                except (ZeroDivisionError, TypeError):
                    continue
            if not overlap:
                merged = pd.concat([merged, domain.to_frame().T], ignore_index=True)
        if not merged.empty:
            merged = merged.drop(columns='width', errors='ignore')
        return merged
    except Exception as e:
        print(f'[ERROR] Exception in safe_merge_similar_domains: {e}')
        print(f'[ERROR] Domains shape: {domains.shape}')
        print(f'[ERROR] Domains columns: {list(domains.columns)}')
        return domains

def assign_domain_colors(domains, optimize_domain_colors=True):
    if domains.empty:
        return domains
    domains = domains.copy()
    unique_ids = domains['proteinDomainID'].unique()
    if optimize_domain_colors:
        colors = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())
        color_map = {uid: colors[i % len(colors)] for i, uid in enumerate(unique_ids)}
        domains['color'] = domains['proteinDomainID'].map(color_map)
    else:
        domains['color'] = domains['color'].fillna('#CCCCCC')
    return domains

def normalize_lengths(coding_exons1, coding_exons2, retained_domains1, retained_domains2):
    try:
        coding_exons1 = coding_exons1.copy() if not coding_exons1.empty else pd.DataFrame()
        coding_exons2 = coding_exons2.copy() if not coding_exons2.empty else pd.DataFrame()
        retained_domains1 = retained_domains1.copy() if not retained_domains1.empty else pd.DataFrame()
        retained_domains2 = retained_domains2.copy() if not retained_domains2.empty else pd.DataFrame()
        for df, name in [(coding_exons1, 'exons1'), (coding_exons2, 'exons2')]:
            if not df.empty:
                if 'start' not in df.columns or 'end' not in df.columns:
                    print(f"[ERROR] Missing start/end columns in {name}: {list(df.columns)}")
                    return coding_exons1, coding_exons2, retained_domains1, retained_domains2, 0, 0
                df['start'] = pd.to_numeric(df['start'], errors='coerce')
                df['end'] = pd.to_numeric(df['end'], errors='coerce')
                df['length'] = df['end'] - df['start'] + 1
                df['length'] = df['length'].fillna(0).clip(lower=0)
        for df, name in [(retained_domains1, 'domains1'), (retained_domains2, 'domains2')]:
            if not df.empty:
                if 'start' in df.columns and 'Start' not in df.columns:
                    df['Start'] = df['start']
                if 'end' in df.columns and 'End' not in df.columns:
                    df['End'] = df['end']
                if 'Start' not in df.columns or 'End' not in df.columns:
                    print(f"[ERROR] Missing Start/End columns in {name}: {list(df.columns)}")
                    available_cols = list(df.columns)
                    start_candidates = [col for col in available_cols if 'start' in col.lower()]
                    end_candidates = [col for col in available_cols if 'end' in col.lower()]
                    if start_candidates and end_candidates:
                        df['Start'] = pd.to_numeric(df[start_candidates[0]], errors='coerce')
                        df['End'] = pd.to_numeric(df[end_candidates[0]], errors='coerce')
                        print(f"[INFO] Created Start/End from {start_candidates[0]}/{end_candidates[0]}")
                    else:
                        print(f"[WARNING] Cannot create Start/End columns for {name}, skipping normalization")
                        continue
                df['Start'] = pd.to_numeric(df['Start'], errors='coerce')
                df['End'] = pd.to_numeric(df['End'], errors='coerce')
        coding_length1 = coding_exons1['length'].sum() if not coding_exons1.empty and 'length' in coding_exons1.columns else 0
        coding_length2 = coding_exons2['length'].sum() if not coding_exons2.empty and 'length' in coding_exons2.columns else 0
        total_length = coding_length1 + coding_length2
        if total_length == 0:
            print("[WARNING] Total coding length is 0, skipping normalization")
            return coding_exons1, coding_exons2, retained_domains1, retained_domains2, 0, 0
        if not coding_exons1.empty and 'length' in coding_exons1.columns:
            coding_exons1['length'] = coding_exons1['length'] / total_length
        if not coding_exons2.empty and 'length' in coding_exons2.columns:
            coding_exons2['length'] = coding_exons2['length'] / total_length
        if not retained_domains1.empty and 'Start' in retained_domains1.columns and 'End' in retained_domains1.columns:
            retained_domains1['Start'] = retained_domains1['Start'] / total_length
            retained_domains1['End'] = retained_domains1['End'] / total_length
        if not retained_domains2.empty and 'Start' in retained_domains2.columns and 'End' in retained_domains2.columns:
            retained_domains2['Start'] = retained_domains2['Start'] / total_length
            retained_domains2['End'] = retained_domains2['End'] / total_length
        return coding_exons1, coding_exons2, retained_domains1, retained_domains2, coding_length1, coding_length2
    except Exception as e:
        print(f"[ERROR] Exception in normalize_lengths: {e}")
        print(f"[ERROR] Exons1 shape: {coding_exons1.shape if not coding_exons1.empty else 'empty'}")
        print(f"[ERROR] Exons2 shape: {coding_exons2.shape if not coding_exons2.empty else 'empty'}")
        print(f"[ERROR] Domains1 shape: {retained_domains1.shape if not retained_domains1.empty else 'empty'}")
        print(f"[ERROR] Domains2 shape: {retained_domains2.shape if not retained_domains2.empty else 'empty'}")
        return coding_exons1, coding_exons2, retained_domains1, retained_domains2, 0, 0

def compute_max_overlap(domains):
    if domains.empty:
        return 1
    scale = 10_000_000
    starts = (domains['Start'] * scale).astype(int)
    ends = (domains['End'] * scale).astype(int)
    points = []
    for s, e in zip(starts, ends):
        points.append((s, 1))
        points.append((e, -1))
    points.sort()
    current = max_overlap = 0
    for _, delta in points:
        current += delta
        max_overlap = max(max_overlap, current)
    return max(1, max_overlap)

def nest_domains(domains):
    if domains.empty:
        return domains
    domains = domains.copy()
    if 'Start' not in domains.columns or 'End' not in domains.columns:
        return domains
    domains['Start'] = pd.to_numeric(domains['Start'], errors='coerce')
    domains['End']   = pd.to_numeric(domains['End'],   errors='coerce')
    domains = domains.dropna(subset=['Start', 'End'])
    if domains.empty:
        return domains
    domains = domains.reset_index(drop=True)
    n = len(domains)
    adj = [[] for _ in range(n)]
    for i in range(n):
        si, ei = domains.loc[i, 'Start'], domains.loc[i, 'End']
        for j in range(i + 1, n):
            sj, ej = domains.loc[j, 'Start'], domains.loc[j, 'End']
            if min(ei, ej) > max(si, sj):
                adj[i].append(j)
                adj[j].append(i)
    comp_idx = [-1] * n
    comps = []
    cid = 0
    for i in range(n):
        if comp_idx[i] != -1:
            continue
        stack = [i]
        comp_idx[i] = cid
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if comp_idx[v] == -1:
                    comp_idx[v] = cid
                    stack.append(v)
        comps.append(comp)
        cid += 1
    ys = [0.0] * n
    hs = [1.0] * n
    for comp in comps:
        if len(comp) == 1:
            idx = comp[0]
            ys[idx] = 0.0
            hs[idx] = 1.0
        else:
            k = len(comp)
            h = 1.0 / k
            widths = (domains.loc[comp, 'End'] - domains.loc[comp, 'Start']).to_dict()
            for level, idx in enumerate(sorted(comp, key=lambda i: -widths[i])):
                ys[idx] = level * h
                hs[idx] = h
    domains['y'] = ys
    domains['height'] = hs
    return domains

def _ensure_exon_numbers(df, strand):
    import pandas as pd, re
    if df is None or df.empty:
        return df
    if 'exon_number' in df.columns:
        df['exon_number'] = pd.to_numeric(df['exon_number'], errors='coerce')
    elif 'attributes' in df.columns:
        df['exon_number'] = df['attributes'].astype(str).str.extract(r'exon_number\s+"?(\d+)"?')[0]
        df['exon_number'] = pd.to_numeric(df['exon_number'], errors='coerce')
    order_idx = df.sort_values('start', ascending=(strand == '+')).index
    fill = pd.Series(range(1, len(order_idx) + 1), index=order_idx)
    df.loc[df['exon_number'].isna() if 'exon_number' in df.columns else order_idx, 'exon_number'] = fill
    df['exon_number'] = df['exon_number'].astype(int)
    return df

def _annotate_retained_with_exon_numbers(retained_df, full_exons):
    import pandas as pd, re
    if retained_df is None or retained_df.empty or full_exons is None or full_exons.empty:
        return retained_df.assign(exon_number=pd.Series(dtype='int'))
    if 'exon_number' not in full_exons.columns:
        if 'attributes' in full_exons.columns:
            full_exons = full_exons.copy()
            full_exons['exon_number'] = full_exons['attributes'].apply(
                lambda x: (re.search(r'exon_number\s+"?(\d+)"?', str(x)).group(1)
                        if re.search(r'exon_number\s+"?(\d+)"?', str(x)) else ''))
        else:
            full_exons = full_exons.copy()
            full_exons['exon_number'] = ''
    ex = full_exons[['start','end','exon_number']].copy()
    ex['start'] = ex['start'].astype(int); ex['end'] = ex['end'].astype(int)
    nums = []
    for _, r in retained_df.iterrows():
        s, e = int(r['start']), int(r['end'])
        hit = ex[(ex['start'] <= s) & (ex['end'] >= e)]
        nums.append(int(hit['exon_number'].iloc[0]) if not hit.empty else None)
    out = retained_df.copy()
    out['exon_number'] = nums
    return out

def filter_domains_by_exons(domains, exons):
    if exons is None or exons.df.empty:
        return pd.DataFrame()
    intersected = domains.intersect(exons)
    return intersected.df.copy()

def order_exons_for_visualization(exons, gene_strand):
    if exons.empty:
        return exons
    if gene_strand == '-':
        return exons.sort_values('start', ascending=False).reset_index(drop=True)
    else:
        return exons.sort_values('start', ascending=True).reset_index(drop=True)

def draw_coverage_bars(ax, cov_left, cov_right, junction_x, exon_y, exon_height, color_left, color_right, args):
    if not args.alignments:
        return
    max_split_cov_left = max(np.max(cov_left) if len(cov_left) > 0 else 0, 1)
    max_split_cov_right = max(np.max(cov_right) if len(cov_right) > 0 else 0, 1)
    cov_y = exon_y + exon_height / 2 + 0.01
    cov_height = 0.15
    bar_width = 0.004
    for i in range(len(cov_left)):
        if cov_left[i] > 0:
            h = (cov_left[i] / max_split_cov_left) * cov_height
            x = junction_x - (i + 1) * bar_width
            ax.add_patch(plt.Rectangle((x, cov_y), bar_width, h, color=color_left, alpha=1, zorder=2))
    for i in range(len(cov_right)):
        if cov_right[i] > 0:
            h = (cov_right[i] / max_split_cov_right) * cov_height
            x = junction_x + i * bar_width
            ax.add_patch(plt.Rectangle((x, cov_y), bar_width, h, color=color_right, alpha=1, zorder=2))
    ax.plot([junction_x - 0.05, junction_x - 0.05], [cov_y, cov_y + cov_height + 0.05], lw=0.6, color='#808080', zorder=3)
    ax.text(junction_x - 0.05, cov_y + cov_height - 0.055, f'{int(max_split_cov_left)}', ha='center', va='bottom', fontsize=9, color=color_left)
    ax.text(junction_x + 0.05, cov_y + cov_height - 0.055, f'{int(max_split_cov_right)}', ha='center', va='bottom', fontsize=9, color=color_right)
    txt = ax.text(junction_x, cov_y + cov_height + 0.01, 'Split Read Coverage', ha='center', va='bottom', fontsize=10, fontweight='bold', color="#F400A1")
    txt.set_path_effects([pe.withStroke(linewidth=0.25, foreground='#FFFFFF')])

def adjust_lightness(hex_color):
    if hex_color.upper() == '#FFFF80':
        return '#B6A400'
    if hex_color.upper() == '#80FFFF':
        return '#00A4B3'
    hex_color = hex_color.lstrip('#')
    r = int(hex_color[0:2], 16) / 255.0
    g = int(hex_color[2:4], 16) / 255.0
    b = int(hex_color[4:6], 16) / 255.0
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l = max(0, min(1, l * 0.55))
    r2, g2, b2 = colorsys.hls_to_rgb(h, l, s)
    return "#{:02X}{:02X}{:02X}".format(int(r2*255), int(g2*255), int(b2*255))

def plot_protein_structure_to_axis(fusions_norm, annotation, protein_domains, args, ax):
    ax.axis('off')
    annotation = fix_missing_type_column(annotation)
    annotation['contig'] = annotation['contig'].apply(remove_chr)
    annotation['geneName'] = annotation['attributes'].apply(lambda x: parse_gtf_attribute('gene_name', x))
    annotation['transcript'] = annotation['attributes'].apply(lambda x: parse_gtf_attribute('transcript_id', x))
    annotation['transcript'] = annotation['transcript'].str.replace(r'\.\d+$', '', regex=True)
    for idx, row in fusions_norm.iterrows():
        fusion_data = {
            'FUSION_BC': row.get('FUSION_BC', ''),
            'gene1': row.get('gene_id1', 'Gene1'),
            'gene2': row.get('gene_id2', 'Gene2'),
            'Breakpoint1': int(row.get('Breakpoint1', 0)),
            'Breakpoint2': int(row.get('Breakpoint2', 0)),
            'contig1': row.get('contig1', ''),
            'contig2': row.get('contig2', ''),
            'site1': row.get('site1', ''),
            'site2': row.get('site2', ''),
            'strand1': row.get('strand1', '+'),
            'strand2': row.get('strand2', '+'),
            'reading_frame': row.get('reading_frame', 'unclear'),
            'type': row.get('type', 'TRA'),
            'direction1': row['direction1'],
            'direction2': row['direction2'],
            'suport': row['max_split_cnt'],
            'somatic_flags': row.get('somatic_flags', '.')
        }
        split_cov1 = np.array([])
        split_cov2 = np.array([])
        if args.alignments:
            try:
                split_cov1 = count_split_reads(args.alignments, row['contig1'], fusion_data['Breakpoint1'], 10)
                split_cov2 = count_split_reads(args.alignments, row['contig2'], fusion_data['Breakpoint2'], 10)
            except Exception:
                split_cov1 = np.array([])
                split_cov2 = np.array([])
        tx1, exons1, coding_exons1, cds_intervals1, retained_exons1 = get_unified_transcript_data(
            annotation, fusion_data['contig1'], fusion_data['gene1'], fusion_data['direction1'],
            fusion_data['Breakpoint1'], row['transcript_id1'], args.transcriptSelection, args.alignments)
        tx2, exons2, coding_exons2, cds_intervals2, retained_exons2 = get_unified_transcript_data(
            annotation, fusion_data['contig2'], fusion_data['gene2'], fusion_data['direction2'],
            fusion_data['Breakpoint2'], row['transcript_id2'], args.transcriptSelection, args.alignments)
        if exons1.empty and exons2.empty:
            ax.text(0.5, 0.5, "No exons retained in fusion.", ha='center', va='center', fontsize=11)
            ax.axis('off')
            continue
        is_non_coding = (not cds_intervals1 and not cds_intervals2)
        gene1_strand = exons1['strand'].iloc[0] if not exons1.empty and 'strand' in exons1.columns else '+'
        gene2_strand = exons2['strand'].iloc[0] if not exons2.empty and 'strand' in exons2.columns else '+'
        exons1 = _ensure_exon_numbers(exons1, gene1_strand)
        exons2 = _ensure_exon_numbers(exons2, gene2_strand)
        def _intersect_exons(exons_df, intervals):
            if not intervals or exons_df.empty:
                if exons_df.empty:
                    return pd.DataFrame(columns=['start','end','strand'])
                return exons_df[['start','end','strand']].copy()
            rows = []
            for _, ex in exons_df.iterrows():
                es, ee = int(ex['start']), int(ex['end'])
                for (cs, ce) in intervals:
                    s = max(es, int(cs))
                    e = min(ee, int(ce))
                    if s < e:
                        rows.append({'start': s, 'end': e, 'strand': ex.get('strand', '+')})
            return pd.DataFrame(rows, columns=['start','end','strand'])
        coding_exons1 = _intersect_exons(exons1, cds_intervals1)
        coding_exons2 = _intersect_exons(exons2, cds_intervals2)
        retained_exons1 = determine_retained_from_breakpoint(coding_exons1, fusion_data['Breakpoint1'], row['direction1'], gene1_strand)
        retained_exons2 = determine_retained_from_breakpoint(coding_exons2, fusion_data['Breakpoint2'], row['direction2'], gene2_strand)
        if not retained_exons1.empty and 'contig' not in retained_exons1.columns:
            retained_exons1['contig'] = fusion_data['contig1']
        if not retained_exons2.empty and 'contig' not in retained_exons2.columns:
            retained_exons2['contig'] = fusion_data['contig2']
        retained_exons1 = _annotate_retained_with_exon_numbers(retained_exons1, exons1)
        retained_exons2 = _annotate_retained_with_exon_numbers(retained_exons2, exons2)
        exons_left, exons_right = retained_exons1, retained_exons2
        strand_left, strand_right = gene1_strand, gene2_strand
        color_left, color_right = args.color1, args.color2
        strand_left_for_plot = strand_left
        strand_right_for_plot = strand_right
        exons_left  = order_exons_for_visualization(exons_left, strand_left_for_plot)
        exons_right = order_exons_for_visualization(exons_right, strand_right_for_plot)
        if not exons_left.empty:
            exons_left['length'] = exons_left['end'] - exons_left['start'] + 1
        if not exons_right.empty:
            exons_right['length'] = exons_right['end'] - exons_right['start'] + 1
        pr_exons_left  = PyRanges(exons_left.rename( columns={'contig': 'Chromosome', 'start': 'Start', 'end': 'End', 'strand': 'Strand'})) if not exons_left.empty else None
        pr_exons_right = PyRanges(exons_right.rename(columns={'contig': 'Chromosome', 'start': 'Start', 'end': 'End', 'strand': 'Strand'})) if not exons_right.empty else None
        pr_domains     = PyRanges(protein_domains.rename(columns={'contig': 'Chromosome', 'start': 'Start', 'end': 'End', 'strand': 'Strand'}))
        domains_left  = filter_domains_by_exons(pr_domains, pr_exons_left) if pr_exons_left is not None else pd.DataFrame()
        domains_right = filter_domains_by_exons(pr_domains, pr_exons_right) if pr_exons_right is not None else pd.DataFrame()
        domains_left  = remove_introns_from_domains(exons_left, domains_left)
        domains_right = remove_introns_from_domains(exons_right, domains_right)
        domains_left  = merge_similar_domains(domains_left,  getattr(args, 'mergeDomainsOverlappingBy', 0.5))
        domains_right = merge_similar_domains(domains_right, getattr(args, 'mergeDomainsOverlappingBy', 0.5))
        domains_left  = assign_domain_colors( domains_left,  getattr(args, 'optimizeDomainColors', True))
        domains_right = assign_domain_colors( domains_right, getattr(args, 'optimizeDomainColors', True))
        if getattr(args, 'swapDomainColors', False):
            if not domains_left.empty:
                domains_left['color']  = adjust_lightness(args.color2)
            if not domains_right.empty:
                domains_right['color'] = adjust_lightness(args.color1)
        exons_left, exons_right, domains_left, domains_right, coding_length_left, coding_length_right = normalize_lengths(
            exons_left, exons_right, domains_left, domains_right
        )
        domains_left  = nest_domains(domains_left)
        domains_right = nest_domains(domains_right)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.2)
        ax.axis('off')
        exon_height = 0.25
        exons_y = 0.5
        padding = 0.03
        usable_height = exon_height - 2 * padding
        if not domains_left.empty:
            domains_left['y'] = exons_y - exon_height / 2 + padding + usable_height * domains_left['y']
            domains_left['height'] = domains_left['height'] * ( usable_height * 0.7 )
        if not domains_right.empty:
            domains_right['y'] = exons_y - exon_height / 2 + padding + usable_height * domains_right['y']
            domains_right['height'] = domains_right['height'] * ( usable_height * 0.7 )
        if coding_length_left > 0 and not exons_left.empty and 'length' in exons_left.columns:
            rect_start = 0
            rect_width = exons_left['length'].sum()
            rect_y = exons_y - exon_height / 2
            grad_res_x = max(int(rect_width * 1000), 10)
            grad_res_y = 150
            grad_y = np.linspace(0, 1, grad_res_y)
            grad_y = np.clip(1.0 - 2.0 * np.abs(grad_y - 0.8) * 1.5, 0.0, 1.0)
            grad = np.tile(grad_y[:, None], (1, grad_res_x))
            ax.imshow(grad, extent=[rect_start, rect_start + rect_width, rect_y, rect_y + exon_height],
                     aspect='auto', cmap='Greys', alpha=1, zorder=1)
            ax.add_patch(plt.Rectangle((rect_start, rect_y), rect_width, exon_height,
                                       color=color_left, alpha=0.55, ec=None, lw=None, zorder=2))
        if coding_length_right > 0 and not exons_right.empty and 'length' in exons_right.columns:
            offset = exons_left['length'].sum() if not exons_left.empty and 'length' in exons_left.columns else 0
            rect_start = offset
            rect_width = exons_right['length'].sum()
            rect_y = exons_y - exon_height / 2
            grad_res_x = max(int(rect_width * 1000), 10)
            grad_res_y = 150
            grad_y = np.linspace(0, 1, grad_res_y)
            grad_y = np.clip(1.0 - 2.0 * np.abs(grad_y - 0.8) * 1.5, 0.0, 1.0)
            grad = np.tile(grad_y[:, None], (1, grad_res_x))
            ax.imshow(grad, extent=[rect_start, rect_start + rect_width, rect_y, rect_y + exon_height],
                     aspect='auto', cmap='Greys', alpha=1, zorder=1)
            ax.add_patch(plt.Rectangle((rect_start, rect_y), rect_width, exon_height,
                                       color=color_right, alpha=0.55, ec=None, lw=None, zorder=2))
        if not exons_left.empty:
            exon_start = 0
            for i, (_, exon_row) in enumerate(exons_left.iterrows()):
                exon_num = exon_row['exon_number'] if 'exon_number' in exon_row and exon_row['exon_number'] else str(i + 1)
                exon_width = exon_row['length'] if 'length' in exon_row else 0
                exon_center_x = exon_start + exon_width / 2
                txt = ax.text(exon_center_x, exons_y + exon_height - 0.21, exon_num,
                             ha='center', va='bottom', fontsize=11, fontweight='bold', color='#FFFFFF', zorder=6)
                txt.set_path_effects([pe.withStroke(linewidth=0.25, foreground='#000000')])
                exon_start += exon_width
        if not exons_right.empty:
            offset = exons_left['length'].sum() if not exons_left.empty and 'length' in exons_left.columns else 0
            exon_start = offset
            for i, (_, exon_row) in enumerate(exons_right.iterrows()):
                exon_num = exon_row['exon_number'] if 'exon_number' in exon_row and exon_row['exon_number'] else str(i + 1)
                exon_width = exon_row['length'] if 'length' in exon_row else 0
                exon_center_x = exon_start + exon_width / 2
                txt = ax.text(exon_center_x, exons_y + exon_height - 0.21, exon_num,
                             ha='center', va='bottom', fontsize=11, fontweight='bold', color='#FFFFFF', zorder=6)
                txt.set_path_effects([pe.withStroke(linewidth=0.25, foreground='#000000')])
                exon_start += exon_width
        length_values = []
        if not exons_left.empty and 'length' in exons_left.columns:
            length_values.extend(exons_left['length'].values)
        if not exons_right.empty and 'length' in exons_right.columns:
            length_values.extend(exons_right['length'].values)
        if length_values:
            exon_boundaries = np.cumsum(length_values)
            for exon_boundary in exon_boundaries[:-1]:
                ax.plot(
                    [exon_boundary, exon_boundary],
                    [exons_y - exon_height / 2, exons_y + exon_height / 2],
                    color='#FFFFFF',
                    linestyle=':',
                    linewidth=0.5,
                    alpha=0.95,
                    zorder=6
                )
        offset_left = 0
        if not domains_left.empty:
            seen = set()
            placed = []
            for domain in domains_left.itertuples():
                x0 = offset_left + domain.Start
                x1 = offset_left + domain.End
                width = x1 - x0
                y = domain.y - 0.02
                h = domain.height
                ax.add_patch(plt.Rectangle((x0, y), width, h, color=domain.color, alpha=0.55, zorder=5))
                _draw_greys_gradient(ax, x0, width, y, h, alpha=1, zorder=4)
                name = domain.proteinDomainName
                if name in seen:
                    continue
                seen.add(name)
                label_x = (x0 + x1) / 2
                va = 'top'
                base_y = y - 0.10
                label_y = base_y
                step = 0.01
                for _ in range(30):
                    conflict = False
                    for (lx, ly, lw, lh) in placed:
                        if abs(label_x - lx) < (lw + 0.05) and abs(label_y - ly) < (lh + 0.03):
                            conflict = True
                            label_y -= step
                            break
                    if not conflict:
                        break
                placed.append((label_x, label_y, len(name) * 0.006, 0.025))
                arrow_start_y = y + 0.01
                arrow_end_y = label_y - 0.01
                label_color = adjust_lightness(domain.color)
                ax.annotate(
                    '',
                    xy=(label_x, arrow_end_y),
                    xytext=(label_x, arrow_start_y),
                    arrowprops=dict(arrowstyle='->', color=label_color, lw=1, alpha=0.7),
                    zorder=5
                )
                ax.text(label_x, label_y, name, ha='center', va=va, fontsize=10, color=label_color, zorder=4)
        offset_right = exons_left['length'].sum() if not exons_left.empty and 'length' in exons_left.columns else 0
        if not domains_right.empty:
            seen = set()
            placed = []
            for domain in domains_right.itertuples():
                x0 = offset_right + domain.Start
                x1 = offset_right + domain.End
                width = x1 - x0
                y = domain.y - 0.02
                h = domain.height
                ax.add_patch(plt.Rectangle((x0, y), width, h, color=domain.color, alpha=0.55, zorder=5))
                _draw_greys_gradient(ax, x0, width, y, h, alpha=1, zorder=4)
                name = domain.proteinDomainName
                if name in seen:
                    continue
                seen.add(name)
                label_x = (x0 + x1) / 2
                va = 'top'
                base_y = y - 0.15
                label_y = base_y
                step = 0.01
                for _ in range(30):
                    conflict = False
                    for (lx, ly, lw, lh) in placed:
                        if abs(label_x - lx) < (lw + 0.05) and abs(label_y - ly) < (lh + 0.03):
                            conflict = True
                            label_y -= step
                            break
                    if not conflict:
                        break
                placed.append((label_x, label_y, len(name) * 0.006, 0.025))
                arrow_start_y = y + 0.01
                arrow_end_y = label_y - 0.01
                label_color = adjust_lightness(domain.color)
                ax.annotate(
                    '',
                    xy=(label_x, arrow_end_y),
                    xytext=(label_x, arrow_start_y),
                    arrowprops=dict(arrowstyle='->', color=label_color, lw=1, alpha=0.7),
                    zorder=5
                )
                ax.text(label_x, label_y, name, ha='center', va=va, fontsize=10, color=label_color, zorder=4)
        info_y = 0.15
        def _local_exons_intersect_intervals(exons_df, intervals):
            if not intervals or exons_df.empty:
                if exons_df.empty:
                    return pd.DataFrame(columns=['start','end','strand'])
                return exons_df[['start','end','strand']].copy()
            rows = []
            for _, ex in exons_df.iterrows():
                es, ee = int(ex['start']), int(ex['end'])
                for (cs, ce) in intervals:
                    s = max(es, int(cs)); e = min(ee, int(ce))
                    if s < e:
                        rows.append({'start': s, 'end': e, 'strand': ex.get('strand', '+')})
            return pd.DataFrame(rows, columns=['start','end','strand'])
        coding_exons1_local = _local_exons_intersect_intervals(exons1, cds_intervals1)
        coding_exons2_local = _local_exons_intersect_intervals(exons2, cds_intervals2)
        retained_exons1_local = retained_from_cds_range(coding_exons1_local, cds_intervals1, row['strand1'], row.get('cds_left_range', '.'))
        retained_exons2_local = retained_from_cds_range(coding_exons2_local, cds_intervals2, row['strand2'], row.get('cds_right_range', '.'))
        lost_exons1 = coding_exons1_local[~coding_exons1_local.apply(lambda x: any((int(x['start']) >= int(r['start']) and int(x['end']) <= int(r['end'])) for _, r in retained_exons1_local.iterrows()), axis=1)] if not retained_exons1_local.empty else coding_exons1_local
        lost_exons2 = coding_exons2_local[~coding_exons2_local.apply(lambda x: any((int(x['start']) >= int(r['start']) and int(x['end']) <= int(r['end'])) for _, r in retained_exons2_local.iterrows()), axis=1)] if not retained_exons2_local.empty else coding_exons2_local
        if not protein_domains.empty:
            retained_domains1_full, lost_domains1_full = get_domains_for_regions(protein_domains, fusion_data['gene1'], retained_exons1_local, lost_exons1)
            retained_domains2_full, lost_domains2_full = get_domains_for_regions(protein_domains, fusion_data['gene2'], retained_exons2_local, lost_exons2)
        else:
            retained_domains1_full = lost_domains1_full = pd.DataFrame()
            retained_domains2_full = lost_domains2_full = pd.DataFrame()
        functional_status, confidence, explanation, status_color, detailed_explanations = determine_functional_status(retained_domains1_full, retained_domains2_full, lost_domains1_full, lost_domains2_full, fusion_data, exons_left if not exons_left.empty else pd.DataFrame(), exons_right if not exons_right.empty else pd.DataFrame())
        summary_lines = get_functional_prediction_summary(functional_status, confidence, explanation, detailed_explanations,
                                                          gene1=fusion_data['gene1'], gene2=fusion_data['gene2'],
                                                          tx1=tx1, tx2=tx2, is_non_coding=is_non_coding)
        total_suport = fusion_data['suport']
        max_split_cov_left  = max(np.max(split_cov1) if len(split_cov1) > 0 else 0, 1)
        max_split_cov_right = max(np.max(split_cov2) if len(split_cov2) > 0 else 0, 1)
        gene_left  = row['gene_id1']
        gene_right = row['gene_id2']
        somatic_flags = row['somatic_flags']
        reading_frame = {
            'in-frame'    : "In-Frame fusion",
            'antisense'   : "Antisense fusion",
            'stop-codon'  : "Stop codon before junction",
            'out-of-frame': "Out-of-Frame fusion"
        }.get(fusion_data['reading_frame'], "Reading frame unclear")
        ax.text( 0.0, info_y - 0.23,  f'{reading_frame}:', ha='left', va='center', fontsize=10, fontweight='bold', color='#2A2A2A')
        ax.text( 0.0, info_y - 0.31, rf'• Split Read Coverage for $\mathit{{{gene_left}}}$:  {max_split_cov_left}',  ha='left', va='center', fontsize=9, color='#2A2A2A' )
        ax.text( 0.0, info_y - 0.39, rf'• Split Read Coverage for $\mathit{{{gene_right}}}$: {max_split_cov_right}', ha='left', va='center', fontsize=9, color='#2A2A2A' )
        ax.text( 0.0, info_y - 0.47,  f'• Total Supporting Reads: {total_suport}', ha='left', va='center', fontsize=10, color='#2A2A2A')
        if somatic_flags == '.':
            ax.text( 0.25, info_y - 0.23,  f'Fusion Not Previously Reported', ha='left', va='center', fontsize=10, fontweight='bold', color='#2A2A2A')
        elif somatic_flags != '0':
            parts = [p.strip() for p in somatic_flags.split(',') if p.strip()]
            if parts:
                grouped = []
                for i in range(0, len(parts), 3):
                    trio = ','.join(parts[i:i+3])
                    grouped.append(trio)
                formatted_flags = ',\n'.join(grouped)
                ax.text( 0.25, info_y - 0.23,  'Fusion Previously Reported:', ha='left', va='center', fontsize=10, fontweight='bold', color='#2A2A2A')
                for i, line in enumerate(grouped):
                    ax.text( 0.26, info_y - 0.31 - i*0.08,  line, ha='left', va='center', fontsize=10, color='#2A2A2A')
        else:
            ax.text( 0.25, info_y - 0.23,  f'Fusion Not Previously Reported', ha='left', va='center', fontsize=10, fontweight='bold', color='#2A2A2A')
        if len(summary_lines) > 3:
            y_offset = info_y - 0.23
            functional_notes_idx = next((i for i, line in enumerate(summary_lines) if line == "Functional Notes:"), 3)
            lines_to_display = summary_lines[functional_notes_idx:functional_notes_idx+6]
            for i, line in enumerate(lines_to_display):
                fontweight = 'bold' if i == 0 else 'normal'
                ax.text(0.65, y_offset - i*0.08, line, ha='left', va='center', fontsize=10, fontweight=fontweight, color='#2A2A2A', usetex=False)

def plot_fusion(ax, exons_left, exons_right, color1, color2, breakpoint1, breakpoint2, x_offset, g1_width, g2_width, direction1, direction2, cds_intervals1, cds_intervals2):
    y0 = 0.0
    y1 = 0.3
    h  = y1 - y0
    yc = y0 + h / 2.0
    def _layout_and_draw(df, width, x_start, color, cds_intervals, breakpoint, direction):
        if df.empty:
            return
        total_len_genome = 0
        for _, r in df.iterrows():
            total_len_genome += int(r["end"]) - int(r["start"])
        coords = []
        curr_x = 0
        for i, r in df.iterrows():
            elen = int(r["end"]) - int(r["start"])
            coords.append((curr_x, elen, int(r["start"]), int(r["end"]), r))
            curr_x += elen
            if i != len(df) - 1:
                curr_x += int(0.02 * total_len_genome)
        if curr_x <= 0:
            curr_x = 1.0
        scale = width / float(curr_x)
        coords_scaled = []
        for local_start, elen, g_start, g_end, row_info in coords:
            xs = x_start + local_start * scale
            w  = elen * scale
            coords_scaled.append((xs, w, g_start, g_end, row_info))
        for xs, w, g_start, g_end, row_info in coords_scaled:
            draw_exon_with_cds_utr(
                ax,
                xs, w,
                g_start, g_end,
                cds_intervals,
                color,
                y0, y1,
                0.5,   # utr_center_frac: keep UTR stripe thinner
                0.5,   # utr_alpha
                0.5,   # cds_alpha
                1      # zorder base
            )
        def _map_pos_to_x(pos):
            for xs, w, gs, ge, _ri in coords_scaled:
                if int(gs) <= int(pos) <= int(ge):
                    glen = max(1, int(ge) - int(gs))
                    return xs + (int(pos) - int(gs)) * (w / glen)
            return None
        coding_parts = []
        if cds_intervals:
            for _, ex in df.iterrows():
                es, ee = int(ex['start']), int(ex['end'])
                sstrand = ex.get('strand', '+')
                for (cs, ce) in cds_intervals:
                    cs, ce = int(cs), int(ce)
                    s = max(es, cs); e = min(ee, ce)
                    if s < e:
                        coding_parts.append({'start': s, 'end': e, 'strand': sstrand})
        coding_df = pd.DataFrame(coding_parts, columns=['start','end','strand']) if coding_parts else pd.DataFrame(columns=['start','end','strand'])
        strand_side = df['strand'].iloc[0] if 'strand' in df.columns and not df.empty else '+'
        retained_coding = determine_retained_from_breakpoint(coding_df, breakpoint, direction, strand_side)
        spans = []
        for _, seg in retained_coding.iterrows():
            xs = _map_pos_to_x(int(seg['start']))
            xe = _map_pos_to_x(int(seg['end']))
            if xs is None or xe is None or xs == xe:
                continue
            spans.append((min(xs, xe), max(xs, xe)))
        use_dotted = False
        if not spans:
            retained_exons = determine_retained_from_breakpoint(df, breakpoint, direction, strand_side)
            for _, seg in retained_exons.iterrows():
                xs = _map_pos_to_x(int(seg['start']))
                xe = _map_pos_to_x(int(seg['end']))
                if xs is None or xe is None or xs == xe:
                    continue
                spans.append((min(xs, xe), max(xs, xe)))
            use_dotted = True
        if spans:
            x_left = min(s[0] for s in spans)
            x_right = max(s[1] for s in spans)
            forward_in_drawn_axis = True
            x0, x1 = (x_left, x_right) if forward_in_drawn_axis else (x_right, x_left)
            length = abs(x1 - x0)
            if length > 1e-6:
                BODY_LW = 3
                HEAD_LEN_PT = 0.1
                HEAD_W_PT = 0.05
                ms = max(12, min(48, 240 * length))
                y_underline = y0 - 0.045
                arrow = mpatches.FancyArrowPatch(
                    (x0, y_underline), (x1, y_underline),
                    arrowstyle=f'-|>,head_length={HEAD_LEN_PT},head_width={HEAD_W_PT}',
                    mutation_scale=ms,
                    linewidth=BODY_LW,
                    color=color,
                    zorder=10,
                    shrinkA=0, shrinkB=0,
                    capstyle='butt', joinstyle='miter',
                    linestyle='--' if use_dotted else '-'
                )
                arrow.set_clip_on(False)
                ax.add_patch(arrow)
        for idx, (xs, w, g_start, g_end, row_info) in enumerate(coords_scaled):
            if "exon_number" in row_info:
                exon_raw = row_info["exon_number"]
            elif "attributes" in row_info:
                m = re.search(r'exon_number\\s+"?(\\d+)"?', str(row_info["attributes"]))
                exon_raw = m.group(1) if m else str(idx + 1)
            else:
                exon_raw = str(idx + 1)
            try:
                exon_int = int(exon_raw)
            except Exception:
                exon_int = idx + 1
            x_center = xs + w / 2.0
            dy = 0.1 * h
            y_text = yc + dy if (exon_int % 2 == 1) else yc - dy
            t = ax.text(
                x_center, y_text, str(exon_raw),
                ha="center", va="center",
                fontsize=11, fontweight="bold",
                color="#FFFFFF", zorder=10
            )
            t.set_path_effects([
                pe.withStroke(linewidth=0.25, foreground="#000000")
            ])
    _layout_and_draw(exons_left,
                     g1_width,
                     x_offset,
                     color1,
                     cds_intervals1,
                     breakpoint=breakpoint1,
                     direction=direction1)
    _layout_and_draw(exons_right,
                     g2_width,
                     x_offset + g1_width,
                     color2,
                     cds_intervals2,
                     breakpoint=breakpoint2,
                     direction=direction2)
    tot = g1_width + g2_width
    ax.set_xlim(x_offset - 0.02 * tot,
                x_offset + tot + 0.02 * tot)
    ax.set_ylim(-0.08, 1.0)
    ax.axis("off")

def _mirror_intervals_within_exon(exon_start, exon_end, intervals):
    out = []
    es, ee = int(exon_start), int(exon_end)
    for s, e in (intervals or []):
        s, e = int(s), int(e)
        cs = max(es, s)
        ce = min(ee, e)
        if cs >= ce:
            out.append((s, e))
            continue
        ms = es + (ee - ce)
        me = es + (ee - cs)
        out.append((min(ms, me), max(ms, me)))
    out = sorted(out, key=lambda t: t[0])
    return out

def flip_edge_cds_within_edge_exons(cds1, cds2, exons_left_df, exons_right_df, strand1, strand2):
    out1 = list(cds1 or [])
    out2 = list(cds2 or [])
    out1 = [(min(int(a), int(b)), max(int(a), int(b))) for a,b in out1]
    out2 = [(min(int(a), int(b)), max(int(a), int(b))) for a,b in out2]
    def _edge_exon_for_side(cds_intervals, exons_df, left_side):
        if not cds_intervals or exons_df is None or exons_df.empty:
            return None, None
        if 'start' not in exons_df.columns or 'end' not in exons_df.columns:
            return None, None
        ex = exons_df.copy()
        ex['start'] = ex['start'].astype(int)
        ex['end'] = ex['end'].astype(int)
        cds_min = min(int(s) for s,_ in cds_intervals)
        cds_max = max(int(e) for _,e in cds_intervals)
        ex = ex.sort_values('start', ascending=True).reset_index(drop=True)
        it = ex.iterrows() if left_side else ex.iloc[::-1].iterrows()
        for _, r in it:
            rs = int(r['start']); re = int(r['end'])
            if rs <= cds_max and re >= cds_min:
                return rs, re
        r = ex.iloc[0 if left_side else -1]
        return int(r['start']), int(r['end'])
    if strand1 == '-' and not exons_left_df.empty:
        if len(cds1) >= len(exons_left_df):
            es = exons_left_df.iloc[0]['start']; ee = exons_left_df.iloc[0]['end']
            out1 = _mirror_intervals_within_exon(es, ee, out1)
        else:
            es, ee = _edge_exon_for_side(out1, exons_left_df, True)
            if es is not None:
                out1 = _mirror_intervals_within_exon(es, ee, out1)
    if strand2 == '-' and not exons_right_df.empty:
        if len(cds2) >= len(exons_right_df):
            es = exons_right_df.iloc[-1]['start']; ee = exons_right_df.iloc[-1]['end']
            out2 = _mirror_intervals_within_exon(es, ee, out2)
        else:
            es, ee = _edge_exon_for_side(out2, exons_right_df, False)
            if es is not None:
                out2 = _mirror_intervals_within_exon(es, ee, out2)
    return out1, out2

def plot_fusion_gene_structure_to_axis(fusion_row, args, ax_e):
    ax_e.axis('off')
    try:
        annotation = pd.read_csv(args.annotation, sep='\t', comment='#', header=None,
                                names=['contig','src','type','start','end','score','strand','frame','attributes'], dtype=str)
    except Exception as e:
        ax_e.text(0.5, 0.5, f"Error reading annotation: {e}", ha='center', va='center', fontsize=14)
        return
    annotation = fix_missing_type_column(annotation)
    annotation['contig']     = annotation['contig'].apply(remove_chr)
    annotation['geneName']   = annotation['attributes'].apply(lambda x: parse_gtf_attribute('gene_name', x))
    annotation['transcript'] = annotation['attributes'].apply(lambda x: parse_gtf_attribute('transcript_id', x))
    annotation['transcript'] = annotation['transcript'].str.replace(r'\.\d+$', '', regex=True)
    for idx, row in fusion_row.iterrows():
        _, exons1, _, cds1, e1 = get_unified_transcript_data(
            annotation, row['contig1'], row['gene_id1'], row['direction1'],
            row['Breakpoint1'], row['transcript_id1'], args.transcriptSelection, args.alignments)
        _, exons2, _, cds2, e2 = get_unified_transcript_data(
            annotation, row['contig2'], row['gene_id2'], row['direction2'],
            row['Breakpoint2'], row['transcript_id2'], args.transcriptSelection, args.alignments)
        strand1_val = exons1['strand'].iloc[0] if not exons1.empty and 'strand' in exons1.columns else '+'
        strand2_val = exons2['strand'].iloc[0] if not exons2.empty and 'strand' in exons2.columns else '+'
        exons_left  = determine_retained_from_breakpoint(exons1, row['Breakpoint1'], row['direction1'], strand1_val)
        exons_right = determine_retained_from_breakpoint(exons2, row['Breakpoint2'], row['direction2'], strand2_val)
        exons_left  = _annotate_retained_with_exon_numbers(exons_left, exons1)
        exons_right = _annotate_retained_with_exon_numbers(exons_right, exons2)
        if 'exon_number' in exons_left.columns and not exons_left.empty:
            exons_left['_exon_order'] = exons_left['exon_number'].apply(
                lambda x: int(str(x)) if str(x).isdigit() else 9999
            )
            exons_left = exons_left.sort_values('_exon_order').drop(columns=['_exon_order']).reset_index(drop=True)
        if 'exon_number' in exons_right.columns and not exons_right.empty:
            exons_right['_exon_order'] = exons_right['exon_number'].apply(
                lambda x: int(str(x)) if str(x).isdigit() else 9999
            )
            exons_right = exons_right.sort_values('_exon_order').drop(columns=['_exon_order']).reset_index(drop=True)
        cds_intervals1, cds_intervals2 = flip_edge_cds_within_edge_exons(cds1, cds2, exons_left, exons_right, strand1_val, strand2_val)
        if exons1.empty or exons2.empty:
            continue
        exon_len1 = sum(int(r['end']) - int(r['start']) for _, r in exons_left.iterrows())
        exon_len2 = sum(int(r['end']) - int(r['start']) for _, r in exons_right.iterrows())
        fixed_intron_width = 0.01
        g1_introns = (len(exons_left) - 1) * fixed_intron_width
        g2_introns = (len(exons_right) - 1) * fixed_intron_width
        g1_span = exon_len1 + g1_introns
        g2_span = exon_len2 + g2_introns
        total_span = g1_span + g2_span
        gap_frac = 0.025
        gap_span = gap_frac * total_span
        total_panel_span = g1_span + gap_span + g2_span
        g1_width = g1_span / total_panel_span
        g2_width = g2_span / total_panel_span
        plot_fusion(ax_e, exons_left, exons_right, args.color1, args.color2, row['Breakpoint1'], row['Breakpoint2'], 0, g1_width, g2_width, row['direction1'], row['direction2'], cds_intervals1, cds_intervals2)

def parse_protein_domains(protein_domains_file):
    cols = ['contig','src','type','start','end','score','strand','frame','attributes']
    try:
        df = pd.read_csv(protein_domains_file, sep='\t', names=cols, header=None, dtype=str)
    except Exception as e:
        print(f"[ERROR] While reading protein domains file {protein_domains_file}: {e}", file=sys.stderr)
        return pd.DataFrame()
    df['start'] = df['start'].astype(int) - 1
    df['end'] = df['end'].astype(int)
    def parse_attr(attr, key):
        m = re.search(rf'{key}=([^;]+)', attr)
        return m.group(1) if m else ''
    df['proteinDomainName'] = df['attributes'].apply(lambda x: urllib.parse.unquote(parse_attr(x, 'Name')))
    df['proteinDomainID']   = df['attributes'].apply(lambda x: parse_attr(x, 'protein_domain_id'))
    df['color']     = df['attributes'].apply(lambda x: parse_attr(x, 'color'))
    df['gene_id']   = df['attributes'].apply(lambda x: parse_attr(x, 'gene_id'))
    df['gene_name'] = df['attributes'].apply(lambda x: parse_attr(x, 'gene_name'))
    return df

def get_domains_for_regions(domains_df, gene_name, retained_exons, lost_exons):
    if domains_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    gene_domains = domains_df[domains_df['gene_name'] == gene_name].copy()
    if gene_domains.empty:
        return pd.DataFrame(), pd.DataFrame()
    retained_domains = []
    lost_domains = []
    for _, domain in gene_domains.iterrows():
        domain_start, domain_end = int(domain['start']), int(domain['end'])
        retained_overlap = any(
            (domain_start <= int(exon['end'])) & (domain_end >= int(exon['start']))
            for _, exon in retained_exons.iterrows()
        )
        lost_overlap = any(
            (domain_start <= int(exon['end'])) & (domain_end >= int(exon['start']))
            for _, exon in lost_exons.iterrows()
        )
        if retained_overlap:
            retained_domains.append(domain)
        if lost_overlap:
            lost_domains.append(domain)
    return pd.DataFrame(retained_domains), pd.DataFrame(lost_domains)

def determine_functional_status(retained_domains1, retained_domains2, lost_domains1, lost_domains2, fusion_data, exons1, exons2):
    critical_domains = {
        'enzymatic': [
            'kinase', 'phosphatase', 'protease', 'transferase', 'hydrolase',
            'ligase', 'lyase',  'isomerase',  'oxidoreductase', 'catalytic',
            'active_site', 'active site', 'enzyme', 'peptidase', 'nuclease',
            'polymerase', 'helicase', 'atpase', 'gtpase'
        ],
        'binding': [
            'dna_binding', 'dna-binding', 'rna_binding', 'rna-binding',
            'protein_binding', 'ligand_binding', 'substrate_binding',
            'cofactor_binding', 'metal_binding', 'zinc_finger', 'helix_turn_helix',
            'leucine_zipper', 'basic_helix_loop_helix'
        ],
        'structural': [
            'transmembrane', 'signal_peptide', 'coiled_coil', 'immunoglobulin',
            'fibronectin', 'collagen', 'elastin', 'keratin'
        ],
        'regulatory': [
            'sh2', 'sh3', 'ph', 'pdz', 'wd40', 'ankyrin', 'tpr', 'armadillo',
            'leucine_rich_repeat', 'death_domain', 'sam_domain'
        ]
    }
    def get_domain_names(domains_df):
        if domains_df.empty:
            return []
        return [str(name).lower().replace(' ', '_') for name in domains_df['proteinDomainName'].fillna('')]
    retained_names1 = get_domain_names(retained_domains1)
    retained_names2 = get_domain_names(retained_domains2)
    lost_names1 = get_domain_names(lost_domains1)
    lost_names2 = get_domain_names(lost_domains2)
    all_retained = retained_names1 + retained_names2
    all_lost = lost_names1 + lost_names2
    functionality_score = 100
    penalties = []
    explanations = []
    reading_frame = fusion_data.get('reading_frame', 'unclear')
    if reading_frame == 'out-of-frame':
        functionality_score -= 80
        penalties.append("Out-of-frame fusion (major penalty)")
        explanations.append("Out-of-frame disruption likely prevents proper translation")
    elif reading_frame == 'stop-codon':
        functionality_score -= 90
        penalties.append("Premature stop codon (critical penalty)")
        explanations.append("Stop codon before fusion junction terminates translation")
    elif reading_frame == 'antisense':
        functionality_score -= 95
        penalties.append("Antisense fusion (critical penalty)")
        explanations.append("Antisense orientation prevents meaningful protein production")
    elif reading_frame == 'in-frame':
        explanations.append("In-frame fusion likely maintains proper translation")
    else:
        functionality_score -= 30
        penalties.append("Unknown reading frame (moderate penalty)")
        explanations.append("Reading frame status unclear")
    critical_lost = []
    for category, domains in critical_domains.items():
        for domain in domains:
            matching_lost = [name for name in all_lost if domain in name]
            if matching_lost:
                critical_lost.extend([(domain, category, matching_lost)])
    if critical_lost:
        enzymatic_lost  = sum(1 for _, cat, _ in critical_lost if cat == 'enzymatic' )
        binding_lost    = sum(1 for _, cat, _ in critical_lost if cat == 'binding'   )
        structural_lost = sum(1 for _, cat, _ in critical_lost if cat == 'structural')
        regulatory_lost = sum(1 for _, cat, _ in critical_lost if cat == 'regulatory')
        if enzymatic_lost > 0:
            penalty = min(50, enzymatic_lost * 25)
            functionality_score -= penalty
            penalties.append(f"Lost {enzymatic_lost} enzymatic domain(s)")
            explanations.append("Loss of catalytic domains severely impacts protein function")
        if binding_lost > 0:
            penalty = min(40, binding_lost * 20)
            functionality_score -= penalty
            penalties.append(f"Lost {binding_lost} binding domain(s)")
            explanations.append("Loss of binding domains affects protein interactions")
        if structural_lost > 0:
            penalty = min(30, structural_lost * 15)
            functionality_score -= penalty
            penalties.append(f"Lost {structural_lost} structural domain(s)")
            explanations.append("Loss of structural domains may affect protein stability")
        if regulatory_lost > 0:
            penalty = min(25, regulatory_lost * 10)
            functionality_score -= penalty
            penalties.append(f"Lost {regulatory_lost} regulatory domain(s)")
            explanations.append("Loss of regulatory domains affects protein regulation")
    critical_retained = []
    for category, domains in critical_domains.items():
        for domain in domains:
            matching_retained = [name for name in all_retained if domain in name]
            if matching_retained:
                critical_retained.append((domain, category))
    if critical_retained:
        enzymatic_retained = sum(1 for _, cat in critical_retained if cat == 'enzymatic')
        binding_retained = sum(1 for _, cat in critical_retained if cat == 'binding')
        if enzymatic_retained > 0:
            explanations.append(f"Retained {enzymatic_retained} enzymatic domain(s)")
        if binding_retained > 0:
            explanations.append(f"Retained {binding_retained} binding domain(s)")
    total_exons1 = len(exons1) if not exons1.empty else 0
    total_exons2 = len(exons2) if not exons2.empty else 0
    if total_exons1 == 0 and total_exons2 == 0:
        functionality_score -= 100
        penalties.append("No coding exons retained")
        explanations.append("Complete loss of coding sequence")
    elif total_exons1 == 0 or total_exons2 == 0:
        functionality_score -= 20
        penalties.append("One gene completely lost")
        explanations.append("Fusion involves only one functional gene")
    fusion_type = fusion_data.get('type', 'Unknown')
    if fusion_type in ['Deletion', 'DEL']:
        functionality_score -= 10
        explanations.append("Deletion may cause loss of regulatory elements")
    elif fusion_type in ['Inversion', 'INV']:
        functionality_score -= 15
        explanations.append("Inversion may disrupt gene structure")
    functionality_score = max(0, min(100, functionality_score))
    if functionality_score >= 70:
        status = "Likely Functional"
        status_color = 'green'
    elif functionality_score >= 40:
        status = "Potentially Functional"
        status_color = 'orange'
    elif functionality_score >= 15:
        status = "Likely Non-Functional"
        status_color = 'red'
    else:
        status = "Non-Functional"
        status_color = 'darkred'
    confidence_factors = []
    if reading_frame in ['in-frame', 'out-of-frame', 'stop-codon', 'antisense']:
        confidence_factors.append("Reading frame determined")
    if len(all_retained) > 0 or len(all_lost) > 0:
        confidence_factors.append("Domain information available")
    if total_exons1 > 0 and total_exons2 > 0:
        confidence_factors.append("Exon structure analyzed")
    confidence = len(confidence_factors) / 3.0
    explanation_text = f"Score: {functionality_score:.0f}/100"
    if penalties:
        explanation_text += f" | Penalties: {'; '.join(penalties[:3])}"
    if len(penalties) > 3:
        explanation_text += f" (+{len(penalties)-3} more)"
    return status, confidence, explanation_text, status_color, explanations[:5]

def get_functional_prediction_summary(status, confidence, explanation, explanations, gene1=None, gene2=None, tx1=None, tx2=None, is_non_coding=False):
    confidence_text = f"Confidence: {'High' if confidence > 0.8 else 'Medium' if confidence > 0.5 else 'Low'}"
    summary_lines = [
        f"Predicted Status: {status}",
        confidence_text,
        explanation
    ]
    if explanations or is_non_coding or (gene1 and tx1) or (gene2 and tx2):
        summary_lines.append("Functional Notes:")
        if is_non_coding:
            summary_lines.append("• Non-coding RNA fusion")
        if gene1 and tx1:
            summary_lines.append(rf"• Gene: $\mathit{{{gene1}}}$ Transcript: {tx1}")
        if gene2 and tx2:
            summary_lines.append(rf"• Gene: $\mathit{{{gene2}}}$ Transcript: {tx2}")
        for exp in explanations[:3]:
            summary_lines.append(f"• {exp}")
    return summary_lines

def parse_fusion_frame(fusion_text):
    if pd.isna(fusion_text) or fusion_text == '.' or fusion_text == '-' or fusion_text == 'NA':
        return 'unclear'
    fusion_str = str(fusion_text).lower()
    if 'in-frame' in fusion_str:
        return 'in-frame'
    elif 'isense' in fusion_str:
        return 'antisense'
    elif 'out-of-frame' in fusion_str:
        return 'out-of-frame'
    elif 'stop' in fusion_str:
        return 'stop-codon'
    else:
        return fusion_str

def count_split_reads(bam_path: str, contig: str, pos: int, window: int = 10) -> np.ndarray:
    bam = pysam.AlignmentFile(bam_path, "rb")
    half = int(window) // 2
    start = max(int(pos) - half, 0)
    end = int(pos) + half
    length = end - start
    if length <= 0:
        bam.close()
        return np.zeros(0, dtype=int)
    counts = np.zeros(length, dtype=int)
    for read in bam.fetch(contig, start, end):
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue
        if read.mapping_quality < 10:
            continue
        has_SA = read.has_tag("SA")
        has_significant_softclip = False
        if read.cigartuples:
            for op, length_cig in read.cigartuples:
                if op == 4 and length_cig >= 10:
                    has_significant_softclip = True
                    break
        if has_SA or has_significant_softclip:
            for i in range(read.reference_start, read.reference_end):
                if start <= i < end:
                    counts[i - start] += 1
    bam.close()
    return counts

def _build_linear_exon_layout(exons: pd.DataFrame, x_offset: float, width: float):
    if exons.empty:
        return [], (lambda _g: None), 0.0
    ex = exons.sort_values('start').copy()
    total_exon_len = int((ex['end'].astype(int) - ex['start'].astype(int)).sum())
    gap_bp = max(1, int(0.02 * total_exon_len))
    coords = []
    cur = 0
    idxs = list(ex.index)
    for i, r in enumerate(ex.itertuples(index=False)):
        s, e = int(r.start), int(r.end)
        w = e - s
        coords.append((cur, w, s, e))
        cur += w
        if i != len(idxs) - 1:
            cur += gap_bp
    scale = width / max(cur, 1)
    coords_scaled = [(x_offset + s * scale, w * scale, gs, ge) for (s, w, gs, ge) in coords]
    def g2x(g):
        for (xs, w, gs, ge) in coords_scaled:
            if gs <= g <= ge:
                return xs + (g - gs) / max(ge - gs, 1) * w
        for i in range(len(coords_scaled) - 1):
            xs0, w0, _, ge0 = coords_scaled[i]
            xs1, _, gs1, _ = coords_scaled[i + 1]
            if ge0 < g < gs1:
                return (xs0 + w0 + xs1) / 2.0
        return None
    total_width = (coords_scaled[-1][0] + coords_scaled[-1][1]) - coords_scaled[0][0]
    return coords_scaled, g2x, total_width

def collect_fusion_bridge_strength(bam_path, contig1, bp1, contig2, bp2, window, min_mapq):
    if not bam_path:
        return 0
    w1 = (max(0, bp1 - window), bp1 + window)
    w2 = (max(0, bp2 - window), bp2 + window)
    seen_1 = set()
    seen_2 = set()
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for r in bam.fetch(contig1, w1[0], w1[1]):
            if r.is_unmapped or r.mapping_quality < min_mapq or not r.has_tag('SA'):
                continue
            hit = False
            for aln in r.get_tag('SA').split(';'):
                if not aln.strip():
                    continue
                rname, pstr, strand, cigar, mapq, nm = aln.split(',')[:6]
                try:
                    p = int(pstr); mq = int(mapq)
                except Exception:
                    continue
                if rname == contig2 and mq >= min_mapq and (w2[0] <= p <= w2[1]):
                    hit = True
            if hit:
                seen_1.add(r.query_name)
        for r in bam.fetch(contig2, w2[0], w2[1]):
            if r.is_unmapped or r.mapping_quality < min_mapq or not r.has_tag('SA'):
                continue
            hit = False
            for aln in r.get_tag('SA').split(';'):
                if not aln.strip():
                    continue
                rname, pstr, strand, cigar, mapq, nm = aln.split(',')[:6]
                try:
                    p = int(pstr); mq = int(mapq)
                except Exception:
                    continue
                if rname == contig1 and mq >= min_mapq and (w1[0] <= p <= w1[1]):
                    hit = True
            if hit:
                seen_2.add(r.query_name)
    strength = len(seen_1 | seen_2)
    return strength

def plot_circos_to_axis(fusions_norm, chromosomes_df, args, ax):
    fusion_colors = {'TRA': '#474747', 'Translocation': '#474747', 'DUP': '#45BA4B', 'Duplication': '#45BA4B', 'DEL': '#EE2B2B', 'Deletion': '#EE2B2B', 'INV': '#2B45EE', 'Inversion': '#2B45EE'}
    if chromosomes_df.empty:
        ax.text(0.5, 0.5, "No chromosomes data available", ha='center', va='center', fontsize=14)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        return
    sectors = {}
    for _, row in chromosomes_df.iterrows():
        chrom_name = str(row['chr']).replace('chr', '')
        chrom_length = int(row['end'])
        sectors[chrom_name] = chrom_length
    for idx, fusion in fusions_norm.iterrows():
        try:
            circos = Circos(sectors, space=2, start=0, end=360)
            ColorCycler.set_cmap("hsv")
            chr_names = [s.name for s in circos.sectors]
            colors = ColorCycler.get_color_list(len(chr_names))
            chr_name2color = {name: color for name, color in zip(chr_names, colors)}
            for sector in circos.sectors:
                colour = chr_name2color[sector.name]
                main_track = sector.add_track((75, 100))
                main_track.axis(fc=colour, ec="#000000", lw=0.5)
                sector.text(f"{sector.name}", r=105, size=8, color=chr_name2color[sector.name])
            if args.cytobands and os.path.isfile(args.cytobands):
                try:
                    cytobands = pd.read_csv(args.cytobands, sep='\t', dtype=str)
                    cytobands['contig'] = cytobands['contig'].apply(remove_chr)
                    for sector in circos.sectors:
                        chrom_bands = cytobands[cytobands['contig'] == sector.name]
                        if not chrom_bands.empty:
                            cyto_track = sector.add_track((75, 100))
                            for _, band in chrom_bands.iterrows():
                                start = int(band['start'])
                                end = int(band['end'])
                                giemsa = band.get('giemsa', 'gneg')
                                colors_dict = {'gneg': '#FFFFFF', 'gpos25': '#EEEEEE', 'gpos50': '#BBBBBB', 'gpos75': '#777777', 'gpos100': '#000000', 'gvar': '#FFFFFF', 'stalk': '#C01E27', 'acen': '#D82322'}
                                color = colors_dict.get(giemsa, '#FFFFFF')
                                cyto_track.rect(start, end, r_lim=(75, 100), fc=color, ec='#000000', lw=0.1, alpha=0.75)
                except Exception:
                    pass
            chrom1 = str(fusion['contig1'])
            pos1 = int(fusion['Breakpoint1'])
            chrom2 = str(fusion['contig2'])
            pos2 = int(fusion['Breakpoint2'])
            fusion_type = fusion.get('type', 'TRA')
            gene1 = str(fusion.get('gene_id1', 'Gene1'))
            gene2 = str(fusion.get('gene_id2', 'Gene2'))
            if chrom1 in [s.name for s in circos.sectors] and chrom2 in [s.name for s in circos.sectors]:
                link_color = fusion_colors.get(fusion_type, '#666666')
                region_size = 500000
                region1 = (chrom1, max(0, pos1 - region_size), pos1 + region_size)
                region2 = (chrom2, max(0, pos2 - region_size), pos2 + region_size)
                circos.link(region1, region2, color=link_color, alpha=1, direction=1 if chrom1 != chrom2 else 0, lw=2)
                sector1 = next(s for s in circos.sectors if s.name == chrom1)
                sector2 = next(s for s in circos.sectors if s.name == chrom2)
                bp_track1 = sector1.add_track((30, 70))
                bp_track2 = sector2.add_track((30, 70))
                bp_track1.rect(max(0, pos1 - 1000000), pos1 + 1000000, fc=link_color, alpha=1, r_lim=(30, 70))
                bp_track2.rect(max(0, pos2 - 1000000), pos2 + 1000000, fc=link_color, alpha=1, r_lim=(30, 70))
                if chrom1 == chrom2:
                    gene_text = f"{gene1}-{gene2}"
                    mid_pos = (pos1 + pos2) // 2
                    sector1.text(gene_text, mid_pos, r=120, size=8, color='#000000', style='italic', orientation='horizontal')
                else:
                    sector1.text(gene1, pos1, r=120, size=8, color='#000000', style='italic', orientation='horizontal')
                    sector2.text(gene2, pos2, r=120, size=8, color='#000000', style='italic', orientation='horizontal')
            temp_fig = circos.plotfig(figsize=(2.75, 2.75))
            temp_fig.savefig('temp_circos.png', dpi=300, bbox_inches='tight', transparent=True)
            img = plt.imread('temp_circos.png')
            ax.imshow(img, extent=[0, 1, 0, 1])
            os.remove('temp_circos.png')
            plt.close(temp_fig)
            ax.axis('off')
            top_left_legend  = [mpatches.Patch(color='#45BA4B', label='DUP')]
            top_right_legend = [mpatches.Patch(color='#EE2B2B', label='DEL')]
            btm_left_legend  = [mpatches.Patch(color='#474747', label='TRA')]
            btm_right_legend = [mpatches.Patch(color='#2B45EE', label='INV')]
            legend1 = ax.legend(handles=top_left_legend,  loc='upper left',  bbox_to_anchor=(-0.1, 0.80), frameon=False, fontsize=11)
            ax.add_artist(legend1)
            legend2 = ax.legend(handles=top_right_legend, loc='upper right', bbox_to_anchor=(1.10, 0.80), frameon=False, fontsize=11)
            ax.add_artist(legend2)
            legend3 = ax.legend(handles=btm_left_legend,  loc='lower left',  bbox_to_anchor=(-0.1, 0.20), frameon=False, fontsize=11)
            ax.add_artist(legend3)
            legend4 = ax.legend(handles=btm_right_legend, loc='lower right', bbox_to_anchor=(1.10, 0.20), frameon=False, fontsize=11)
            ax.add_artist(legend4)
        except Exception as e:
            ax.text(0.5, 0.5, f"Error: {str(e)}", ha='center', va='center', fontsize=10)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')

def plot_combined_panels_vector(fusions_norm, annotation, protein_domains, cytobands, chromosomes_df, args, output_pdf):
    with PdfPages(output_pdf) as pdf:
        for idx, row in fusions_norm.iterrows():
            print(f"\nGenerating vector combined layout for fusion {idx + 1}/{len(fusions_norm)}...")
            fig = plt.figure(figsize=(args.pdfWidth, args.pdfHeight))
            try:
                ax_a = fig.add_axes([0.00, 1.00, 1.00, 0.15])  # Second from top x2: Chromosome ideograms
                ax_b = fig.add_axes([0.00, 0.60, 1.00, 0.30])  # Middle: Gene structure
                ax_c = fig.add_axes([0.00, 0.05, 1.00, 0.35])  # Bottom: Protein structure
                ax_d = fig.add_axes([0.25, 0.90, 0.50, 0.40])  # TopCenter: Circos
                ax_e = fig.add_axes([0.25, 0.40, 0.50, 0.25])  # 3rd from Top: Fusion Transcript
                ax_a.set_zorder(1200)
                ax_b.set_zorder(1200)
                ax_c.set_zorder(1200)
                ax_d.set_zorder(1200)
                ax_e.set_zorder(1200)
                ax_a.patch.set_facecolor('none')
                ax_b.patch.set_facecolor('none')
                ax_c.patch.set_facecolor('none')
                ax_d.patch.set_facecolor('none')
                ax_e.patch.set_facecolor('none')
                single_fusion = fusions_norm.iloc[[idx]]
                plot_chromosome_ideograms_to_axis(single_fusion, cytobands, ax_a)
                plot_coverage_gene_structure_to_axis(single_fusion, args, ax_b)
                plot_protein_structure_to_axis(single_fusion, annotation, protein_domains, args, ax_c)
                plot_circos_to_axis(single_fusion, chromosomes_df, args, ax_d)
                plot_fusion_gene_structure_to_axis(single_fusion, args, ax_e)
                strand_left  = row['strand1'] if 'strand1' in row else '+'
                strand_right = row['strand2'] if 'strand2' in row else '+'
                fusion_type = row['type']
                collect_and_draw_zoom_overlays(fig, ax_a, ax_b, ax_e, ax_c, strand_left, strand_right, args.color1, args.color2, fusion_type)
            except Exception as e:
                print(f"[WARNING] Could not render combined layout for fusion {idx + 1}: {e}", file=sys.stderr)
                fig.clf()
                ax_err = fig.add_axes([0, 0, 1, 1])
                ax_err.axis('off')
                gene1 = row['contig1']
                gene2 = row['contig2']
                ax_err.text(0.5, 0.5, f"Error rendering {gene1}-{gene2} event: ", ha='center', va='center', fontsize=14)
            pdf.savefig(fig, dpi=300, bbox_inches='tight', transparent=True)
            plt.close(fig)

def _data_bbox_to_figfrac(fig, ax, bbox):
    (x0p, y0p) = ax.transData.transform((bbox['x0'], bbox['y0']))
    (x1p, y1p) = ax.transData.transform((bbox['x1'], bbox['y1']))
    inv = fig.transFigure.inverted()
    (xf0, yf0) = inv.transform((x0p, y0p))
    (xf1, yf1) = inv.transform((x1p, y1p))
    return {'x0': xf0, 'x1': xf1, 'y0': yf0, 'y1': yf1}

def _collect_bbox_panelA_left(ax):
    target_lines = []
    for line in ax.lines:
        col = line.get_color()
        if isinstance(col, str) and col == '#F400A1':
            xd = np.array(line.get_xdata())
            yd = np.array(line.get_ydata())
            if xd.size == 0:
                continue
            target_lines.append((float(np.mean(xd)), yd.min(), yd.max()))
    if not target_lines:
        xl, xr = ax.get_xlim()
        yb, yt = ax.get_ylim()
        x = xl
        return {'x0': x-0.02, 'x1': x+0.02, 'y0': yb, 'y1': yt}
    x, y0, y1 = sorted(target_lines, key=lambda t: t[0])[0]
    return {'x0': x, 'x1': x, 'y0': y0, 'y1': y1}

def _collect_bbox_panelA_right(ax):
    target_lines = []
    for line in ax.lines:
        col = line.get_color()
        if isinstance(col, str) and col == '#F400A1':
            xd = np.array(line.get_xdata())
            yd = np.array(line.get_ydata())
            if xd.size == 0:
                continue
            target_lines.append((float(np.mean(xd)), yd.min(), yd.max()))
    if not target_lines:
        xl, xr = ax.get_xlim()
        yb, yt = ax.get_ylim()
        x = xr
        y = yb
        return {'x0': x-0.02, 'x1': x+0.02, 'y0': y, 'y1': y}
    x, y0, y1 = sorted(target_lines, key=lambda t: t[0])[-1]
    return {'x0': x, 'x1': x, 'y0': y0, 'y1': y1}

def _collect_bbox_PanelB_arrows(ax, base_color):
    target_rgba = mcolors.to_rgba(base_color)
    for p in ax.patches:
        if isinstance(p, mpatches.FancyArrowPatch):
            ec = p.get_edgecolor()
            if (
                isinstance(ec, tuple)
                and len(ec) >= 3
                and np.allclose(ec[:3], target_rgba[:3], atol=0.05)
            ):
                verts = p.get_path().vertices
                (x0, y0) = verts[0]
                (x1, y1) = verts[2]
                return {
                    'x0': x0,
                    'x1': x1,
                    'y0': y0 - 0.015,
                    'y1': y1 - 0.015
                }
    return None

def _collect_bbox_PanelB_arrows_fig(fig, ax, base_color):
    bbox_ax = _collect_bbox_PanelB_arrows(ax, base_color)
    if bbox_ax is None:
        xl, xr = ax.get_xlim()
        yb, yt = ax.get_ylim()
        corners_data = np.array([
            [xl, yb],
            [xl + 0.5*(xr-xl), yt]
        ])
        corners_disp = ax.transData.transform(corners_data)
        corners_fig = fig.transFigure.inverted().transform(corners_disp)
        x0_fig, y0_fig = corners_fig[0]
        x1_fig, y1_fig = corners_fig[1]
        return {
            'x0': x0_fig,
            'x1': x1_fig,
            'y0': y0_fig,
            'y1': y1_fig
        }
    corners_data = np.array([
        [bbox_ax['x0'], bbox_ax['y0']],
        [bbox_ax['x1'], bbox_ax['y1']]
    ])
    corners_disp = ax.transData.transform(corners_data)
    corners_fig = fig.transFigure.inverted().transform(corners_disp)
    x0_fig, y0_fig = corners_fig[0]
    x1_fig, y1_fig = corners_fig[1]
    return {
        'x0': x0_fig,
        'x1': x1_fig,
        'y0': y0_fig,
        'y1': y1_fig
    }

def _collect_bbox_PanelE_arrows(ax, base_color):
    target_rgba = mcolors.to_rgba(base_color)
    for p in ax.patches:
        if isinstance(p, mpatches.FancyArrowPatch):
            ec = p.get_edgecolor()
            if (
                isinstance(ec, tuple)
                and len(ec) >= 3
                and np.allclose(ec[:3], target_rgba[:3], atol=0.05)
            ):
                verts = p.get_path().vertices
                (x0, y0) = verts[0]
                (x1, y1) = verts[2]
                return {
                    'x0': x0,
                    'x1': x1,
                    'y0': y0 - 0.05,
                    'y1': y1 - 0.05
                }
    return None

def _collect_bbox_PanelE_arrows_fig(fig, ax, base_color):
    bbox_ax = _collect_bbox_PanelE_arrows(ax, base_color)
    if bbox_ax is None:
        xl, xr = ax.get_xlim()
        yb, yt = ax.get_ylim()
        corners_data = np.array([
            [xl, yb],
            [xl + 0.5*(xr-xl), yt]
        ])
        corners_disp = ax.transData.transform(corners_data)
        corners_fig = fig.transFigure.inverted().transform(corners_disp)
        x0_fig, y0_fig = corners_fig[0]
        x1_fig, y1_fig = corners_fig[1]
        return {
            'x0': x0_fig,
            'x1': x1_fig,
            'y0': y0_fig,
            'y1': y1_fig
        }
    corners_data = np.array([
        [bbox_ax['x0'], bbox_ax['y0']],
        [bbox_ax['x1'], bbox_ax['y1']]
    ])
    corners_disp = ax.transData.transform(corners_data)
    corners_fig = fig.transFigure.inverted().transform(corners_disp)
    x0_fig, y0_fig = corners_fig[0]
    x1_fig, y1_fig = corners_fig[1]
    return {
        'x0': x0_fig,
        'x1': x1_fig,
        'y0': y0_fig,
        'y1': y1_fig
    }

def _collect_bbox_panelC_block_fig(fig, ax, base_color):
    bbox_ax = _collect_bbox_panelC_block(ax, base_color)
    corners_data = np.array([
        [bbox_ax['x0'], bbox_ax['y0']],
        [bbox_ax['x1'], bbox_ax['y1']]
    ])
    corners_disp = ax.transData.transform(corners_data)
    corners_fig = fig.transFigure.inverted().transform(corners_disp)
    x0_fig, y0_fig = corners_fig[0]
    x1_fig, y1_fig = corners_fig[1]
    return {
        'x0': x0_fig,
        'x1': x1_fig,
        'y0': y0_fig,
        'y1': y1_fig
    }

def _collect_bbox_panelB_block(ax, base_color):
    target_rgb = mcolors.to_rgba(base_color)[:3]
    xs = []
    ys = []
    for p in ax.patches:
        if isinstance(p, plt.Rectangle):
            fc = p.get_facecolor()
            if (isinstance(fc, tuple) and len(fc) >= 3 and
                np.allclose(fc[:3], target_rgb, atol=0.05)):
                x, y = p.get_xy()
                w = p.get_width()
                h = p.get_height()
                if y + h < 0:
                    continue
                xs += [x, x + w]
                ys += [y, y + h]
    if not xs:
        xl, xr = ax.get_xlim()
        yb, yt = ax.get_ylim()
        return {
            'x0': xl,
            'x1': xl + 0.5*(xr-xl),
            'y0': max(yb, 0.0),
            'y1': yt
        }
    return {
        'x0': min(xs),
        'x1': max(xs),
        'y0': min(ys),
        'y1': max(ys)
    }

def _collect_bbox_panelC_block(ax, base_color):
    target_rgb = mcolors.to_rgba(base_color)[:3]
    xs = []
    ys = []
    for p in ax.patches:
        if isinstance(p, plt.Rectangle):
            fc = p.get_facecolor()
            if (isinstance(fc, tuple) and len(fc) >= 3 and
                np.allclose(fc[:3], target_rgb, atol=0.075)):
                x, y = p.get_xy()
                w = p.get_width()
                h = p.get_height()
                xs += [x, x + w]
                ys += [y, y + h]
    if not xs:
        xl, xr = ax.get_xlim()
        yb, yt = ax.get_ylim()
        return {
            'x0': xl,
            'x1': xl + 0.5*(xr-xl),
            'y0': yb,
            'y1': yt
        }
    return {
        'x0': min(xs),
        'x1': max(xs),
        'y0': min(ys),
        'y1': max(ys)
    }

def _collect_bbox_panelE_block(ax, base_color):
    target_rgb = mcolors.to_rgba(base_color)[:3]
    xs, ys = [], []
    for p in ax.patches:
        if isinstance(p, mpatches.Rectangle):
            fc = p.get_facecolor()
            if (isinstance(fc, tuple) and len(fc) >= 3 and
                np.allclose(fc[:3], target_rgb, atol=0.05)):
                x, y = p.get_xy()
                w, h = p.get_width(), p.get_height()
                xs += [x, x + w]
                ys += [y, y + h]
    if not xs:
        xl, xr = ax.get_xlim()
        yb, yt = ax.get_ylim()
        return {'x0': xl, 'x1': xl + 0.5*(xr-xl), 'y0': yb, 'y1': yt}
    return {'x0': min(xs), 'x1': max(xs), 'y0': min(ys), 'y1': max(ys)}

def _collect_bbox_panelE_blocks_fig(fig, ax, base_color):
    bbox_ax = _collect_bbox_panelE_block(ax, base_color)
    corners_data = np.array([[bbox_ax['x0'], bbox_ax['y0']], [bbox_ax['x1'], bbox_ax['y1']]])
    corners_disp = ax.transData.transform(corners_data)
    corners_fig = fig.transFigure.inverted().transform(corners_disp)
    x0_fig, y0_fig = corners_fig[0]
    x1_fig, y1_fig = corners_fig[1]
    return {'x0': min(x0_fig, x1_fig), 'x1': max(x0_fig, x1_fig),
            'y0': min(y0_fig, y1_fig), 'y1': max(y0_fig, y1_fig)}

def _draw_highlight_box(ax, bbox, fill_color, alpha=0.001, zorder=999):
    x0, y0 = bbox['x0'], bbox['y0']
    w = bbox['x1'] - bbox['x0']
    h = bbox['y1'] - bbox['y0']
    rect = plt.Rectangle(
        (x0, y0),
        w,
        h,
        facecolor=mcolors.to_rgba(fill_color, alpha),
        edgecolor='none',
        linewidth=0,
        zorder=zorder
    )
    return rect

def _draw_ribbon_connector_twist(fig, bbox_start, bbox_end, base_color, twist=False, alpha=0.5, zorder=1100, n_segments=100, n_gradient_steps=100):
    y_end          =   bbox_end['y1']
    end_left       =   bbox_end['x0']
    end_right      =   bbox_end['x1']
    y_start        = bbox_start['y0']
    start_left     = bbox_start['x0']
    start_right    = bbox_start['x1']
    start_inverted = start_left > start_right
    if start_inverted:
        start_left, start_right = start_right, start_left
    start_x = (start_left + start_right) / 2
    end_x = (end_left + end_right) / 2
    width = (start_right - start_left) / 20
    twist_amount = 1.01 if (twist or start_inverted) else 0.0
    t = np.linspace(0, 1, n_segments)
    x_center = start_x + (end_x - start_x) * t
    y_center = y_start + (y_end - y_start) * t
    width_profile = width * (1 - 0.3 * t)
    twist_angle = np.pi * t * twist_amount
    x_left = x_center - width_profile * np.cos(twist_angle)
    y_left = y_center - width_profile * np.sin(twist_angle)
    x_right = x_center + width_profile * np.cos(twist_angle)
    y_right = y_center + width_profile * np.sin(twist_angle)
    generated_width_start = abs(x_right[0] - x_left[0])
    generated_width_end = abs(x_right[-1] - x_left[-1])
    desired_width_start = abs(bbox_start['x1'] - bbox_start['x0'])
    desired_width_end = abs(bbox_end['x1'] - bbox_end['x0'])
    scale_start = desired_width_start / generated_width_start if generated_width_start > 0 else 1
    scale_end = desired_width_end / generated_width_end if generated_width_end > 0 else 1
    for i in range(len(t)):
        scale_factor = scale_start + (scale_end - scale_start) * t[i]
        current_center_x = x_center[i]
        x_left[i] = current_center_x - (current_center_x - x_left[i]) * scale_factor
        x_right[i] = current_center_x + (x_right[i] - current_center_x) * scale_factor
    base_color = to_rgba(base_color) if isinstance(base_color, str) else base_color
    artists = []
    for i in range(len(t) - 1):
        segment_x = [x_left[i], x_left[i+1], x_right[i+1], x_right[i]]
        segment_y = [y_left[i], y_left[i+1], y_right[i+1], y_right[i]]
        mid_x = (segment_x[0] + segment_x[3]) / 2
        mid_y = (segment_y[0] + segment_y[3]) / 2
        den = max(n_gradient_steps - 1, 1)
        for j in range(n_gradient_steps):
            scale = 1 - j * (1.0 / n_gradient_steps)
            frac = j / den
            center_weight = 1.0 - abs(2.0 * frac - 1.0)
            center_weight = max(center_weight, 0.0)
            w = center_weight ** 0.7
            whiteness = 0.3 + 0.5 * w
            gradient_x = [mid_x + scale * (sx - mid_x) for sx in segment_x]
            gradient_y = [mid_y + scale * (sy - mid_y) for sy in segment_y]
            segment_alpha = alpha
            r = base_color[0] * (1.0 - whiteness) + whiteness
            g = base_color[1] * (1.0 - whiteness) + whiteness
            b = base_color[2] * (1.0 - whiteness) + whiteness
            segment_color = (r, g, b, segment_alpha)
            poly = mpatches.Polygon(
                list(zip(gradient_x, gradient_y)),
                closed=True,
                facecolor=segment_color,
                edgecolor='none',
                linewidth=0,
                transform=fig.transFigure,
                zorder=zorder
            )
            fig.add_artist(poly)
            artists.append(poly)
    return artists

def _draw_A2B_connector(fig, bbox_start, bbox_end, base_color, alpha_g=0.25, zorder=1100):
    y_end          =   bbox_end['y1']
    end_left       =   bbox_end['x0']
    end_right      =   bbox_end['x1']
    y_start        = bbox_start['y0'] - 0.03
    start_left     = bbox_start['x0']
    start_right    = bbox_start['x1']
    poly_color = mpatches.Polygon(
        [(start_left,  y_start),
         (start_right, y_start),
         (end_right, y_end),
         (end_left,  y_end)],
        closed=True,
        alpha=alpha_g,
        facecolor=base_color,
        edgecolor='none',
        linewidth=0,
        transform=fig.transFigure,
        zorder=zorder + 1
    )
    return [poly_color]

def _draw_B2E_connector(fig, bbox_start, bbox_end, base_color, twist=False, alpha=0.25, zorder=1100, steps=100):
    y_end          =   bbox_end['y1']
    end_left       =   bbox_end['x0']
    end_right      =   bbox_end['x1']
    y_start        = bbox_start['y0']
    start_left     = bbox_start['x0']
    start_right    = bbox_start['x1']
    start_inverted = start_left > start_right
    if twist or start_inverted:
        patches = _draw_ribbon_connector_twist(fig, bbox_start, bbox_end, base_color, twist=True, zorder=zorder, n_segments=100, n_gradient_steps=steps)
        return patches
    else:
        poly_x = [
            start_left,
            start_right,
            end_right,
            end_left
        ]
        poly_y = [
            y_start,
            y_start,
            y_end,
            y_end
        ]
        tint = mcolors.to_rgba(base_color, alpha)
        poly = mpatches.Polygon(
            list(zip(poly_x, poly_y)),
            closed=True,
            facecolor=tint,
            edgecolor='none',
            linewidth=0,
            transform=fig.transFigure,
            zorder=zorder
        )
    return poly

def _draw_E2C_connector(fig, bbox_start, bbox_end, base_color, twist=False, alpha=0.25, zorder=1100):
    start_left  = bbox_start['x0']
    start_right = bbox_start['x1']
    end_left    = bbox_end['x0']
    end_right   = bbox_end['x1']
    if twist:
        tmp_left  = end_left
        tmp_right = end_right
        end_left  = tmp_right
        end_right = tmp_left
    y_start_edge  = bbox_start['y0']
    y_end_edge    = bbox_end['y1']
    poly_x = [
        start_left,
        start_right,
        end_right,
        end_left
    ]
    poly_y = [
        y_start_edge,
        y_start_edge,
        y_end_edge,
        y_end_edge
    ]
    tint = mcolors.to_rgba(base_color, alpha)
    poly = mpatches.Polygon(
        list(zip(poly_x, poly_y)),
        closed=True,
        facecolor=tint,
        edgecolor='none',
        linewidth=0,
        transform=fig.transFigure,
        zorder=zorder
    )
    return poly

def collect_and_draw_zoom_overlays(fig,
                                   ax_a,
                                   ax_b,
                                   ax_e,
                                   ax_c,
                                   strand_left,
                                   strand_right,
                                   color_left,
                                   color_right,
                                   fusion_type='TRA'
                                   ):
    bboxA_left_ax   = _collect_bbox_panelA_left( ax_a)
    bboxA_right_ax  = _collect_bbox_panelA_right(ax_a)
    bboxA_left_fig  = _data_bbox_to_figfrac(fig, ax_a, bboxA_left_ax )
    bboxA_right_fig = _data_bbox_to_figfrac(fig, ax_a, bboxA_right_ax)
    bboxB_left_ax   = _collect_bbox_panelB_block(ax_b, color_left )
    bboxB_right_ax  = _collect_bbox_panelB_block(ax_b, color_right)
    bboxB_left_fig  = _collect_bbox_PanelB_arrows_fig(fig, ax_b, color_left )
    bboxB_right_fig = _collect_bbox_PanelB_arrows_fig(fig, ax_b, color_right)
    bboxB_left_blk  = _data_bbox_to_figfrac(fig, ax_b, bboxB_left_ax )
    bboxB_right_blk = _data_bbox_to_figfrac(fig, ax_b, bboxB_right_ax)
    bboxC_left_fig  = _collect_bbox_panelC_block_fig( fig, ax_c, color_left )
    bboxC_right_fig = _collect_bbox_panelC_block_fig( fig, ax_c, color_right)
    bboxE_left_ax   = _collect_bbox_panelE_blocks_fig(fig, ax_e, color_left )
    bboxE_right_ax  = _collect_bbox_panelE_blocks_fig(fig, ax_e, color_right)
    bboxE_left_fig  = _collect_bbox_PanelE_arrows_fig(fig, ax_e, color_left )
    bboxE_right_fig = _collect_bbox_PanelE_arrows_fig(fig, ax_e, color_right)
    rectB_left  =  _draw_highlight_box(ax_b, _collect_bbox_panelB_block(ax_b, color_left),  color_left,  alpha=0.1, zorder=1); ax_b.add_patch(rectB_left)
    rectB_right =  _draw_highlight_box(ax_b, _collect_bbox_panelB_block(ax_b, color_right), color_right, alpha=0.1, zorder=1); ax_b.add_patch(rectB_right)
    rectC_left  =  _draw_highlight_box(ax_c, _collect_bbox_panelC_block(ax_c, color_left),  color_left,  alpha=0.1, zorder=1); ax_c.add_patch(rectC_left)
    rectC_right =  _draw_highlight_box(ax_c, _collect_bbox_panelC_block(ax_c, color_right), color_right, alpha=0.1, zorder=1); ax_c.add_patch(rectC_right)
    rectE_left  =  _draw_highlight_box(ax_e, _collect_bbox_panelE_block(ax_e, color_left),  color_left,  alpha=0.1, zorder=1); ax_e.add_patch(rectE_left)
    rectE_right =  _draw_highlight_box(ax_e, _collect_bbox_panelE_block(ax_e, color_right), color_right, alpha=0.1, zorder=1); ax_e.add_patch(rectE_right)
    for artist in  _draw_A2B_connector(fig, bboxA_left_fig,  bboxB_left_blk,  base_color=color_left ): fig.add_artist(artist)
    for artist in  _draw_A2B_connector(fig, bboxA_right_fig, bboxB_right_blk, base_color=color_right): fig.add_artist(artist)
    fig.add_artist(_draw_E2C_connector(fig, bboxE_left_fig,  bboxC_left_fig,  base_color=color_left ))
    fig.add_artist(_draw_E2C_connector(fig, bboxE_right_fig, bboxC_right_fig, base_color=color_right))
    conector_B2E_left  = _draw_B2E_connector(fig, bboxB_left_fig,  bboxE_left_ax,  base_color=color_left )
    conector_B2E_right = _draw_B2E_connector(fig, bboxB_right_fig, bboxE_right_ax, base_color=color_right)
    if isinstance(conector_B2E_left, mpatches.Polygon):
        fig.add_artist(conector_B2E_left)
    else:
        for artist in conector_B2E_left: fig.add_artist(artist)
    if isinstance(conector_B2E_right, mpatches.Polygon):
        fig.add_artist(conector_B2E_right)
    else:
        for artist in conector_B2E_right: fig.add_artist(artist)

def create_error_pdf(output_path, error_message):
    with PdfPages(output_path) as pdf:
        fig = plt.figure(figsize=(11.692, 8.267))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')
        lines = error_message.split('\n')
        y_start = 0.5 + (len(lines) - 1) * 0.05
        for i, line in enumerate(lines):
            ax.text(0.5, y_start - i * 0.1, line, ha='center', va='center',
                    fontsize=16, fontweight='bold', color='#CC0000')
        pdf.savefig(fig, dpi=300, bbox_inches='tight')
        plt.close(fig)

def main():
    args = parse_args()
    print("\nArguments parsed and validated successfully.")
    print("\nLoading Fusions file...")
    def _load_fusions():
        if args.fusionsFormat == 'vcf':
            return parse_vcf_to_df(args.fusions)
        if args.fusionsFormat == 'tsv':
            return parse_tsv_to_df(args.fusions)
        if args.fusionsFormat == 'txt':
            return parse_txt_to_df(args.fusions)
        if args.fusionsFormat == 'cff':
            return parse_cff_to_df(args.fusions)
        fext = os.path.splitext(args.fusions)[1].lower()
        if fext == '.vcf':
            return parse_vcf_to_df(args.fusions)
        elif fext == '.tsv':
            return parse_tsv_to_df(args.fusions)
        elif fext == '.txt':
            return parse_txt_to_df(args.fusions)
        elif fext == '.cff':
            return parse_cff_to_df(args.fusions)
        with open(args.fusions, 'r', encoding='utf-8', errors='ignore') as fh:
            first_line = fh.readline().strip()
            if first_line.startswith('##'):
                return parse_vcf_to_df(args.fusions)
            if 'gene5_chr' in first_line and 'gene3_chr' in first_line:
                return parse_cff_to_df(args.fusions)
            fh.seek(0)
            for line in fh:
                if line.startswith('##') or not line.strip():
                    continue
                if line.startswith('#CHROM') or (line.count('\t') >= 8 and 'INFO' in line):
                    return parse_vcf_to_df(args.fusions)
                break
        return parse_tsv_to_df(args.fusions)
    try:
        fusions = _load_fusions()
    except Exception as e:
        print(f"[ERROR] While reading Fusions file {args.fusions}: {e}", file=sys.stderr)
        sys.exit(1)
    if fusions.empty or len(fusions) == 0:
        print(f"[ERROR] Fusion file is empty (no fusions found).", file=sys.stderr)
        create_error_pdf(args.output, "ERROR: Fusion file empty.\n\nNo fusions found in the input file.")
        print(f"\nError PDF generated: {args.output}")
        return
    print("\nLoading annotation GTF (exons/CDS only)...")
    gtf_cols = ['contig','src','type','start','end','score','strand','frame','attributes']
    try:
        annotation = pd.read_csv(
            args.annotation, sep='\t', comment='#', header=None, names=gtf_cols, dtype=str
        )
    except Exception as e:
        print(f"[ERROR] While reading annotation file {args.annotation}: {e}", file=sys.stderr)
        sys.exit(1)
    annotation = fix_missing_type_column(annotation)
    annotation = annotation[annotation['type'].isin(['exon','CDS'])].copy()
    annotation['contig'] = annotation['contig'].apply(remove_chr)
    annotation['geneID'] = annotation['attributes'].apply(lambda x: parse_gtf_attribute('gene_id', x))
    annotation['geneName']   = annotation['attributes'].apply(lambda x: parse_gtf_attribute('gene_name', x))
    annotation['transcript'] = annotation['attributes'].apply(lambda x: parse_gtf_attribute('transcript_id', x))
    cytobands = None
    if args.cytobands:
        try:
            cytobands = pd.read_csv(args.cytobands, sep='\t', dtype=str)
        except Exception as e:
            print(f"[ERROR] While reading cytobands file {args.cytobands}: {e}", file=sys.stderr)
            cytobands = None
    protein_domains = pd.DataFrame()
    if args.proteinDomains:
        protein_domains = parse_protein_domains(args.proteinDomains)
    chromosomes_df = pd.DataFrame()
    if args.chromosomes:
        try:
            chromosomes_df = pd.read_csv(args.chromosomes)
        except Exception as e:
            print(f"[ERROR] While reading chromosomes file {args.chromosomes}: {e}", file=sys.stderr)
            sys.exit(1)
    print("\nNormalizing and parsing fusions...")
    fusions_norm = normalize_fusions(fusions)
    print("\nGenerating combined layout with all 5 panels...")
    plot_combined_panels_vector(fusions_norm, annotation, protein_domains, cytobands, chromosomes_df, args, args.output)
    print(f"\nVector-based combined PDF generated: {args.output}")
    print("\nLayout from top to bottom:\n\nPanel D - Circos Plot (center)\nPanel A - Chromosome idiograms (right and left)\nPanel B - Gene structure and exomic coverage x 2\nPanel E - Fusion transcript\nPanel C - Fusion protein structure")
    print("\nAll done!")

if __name__ == "__main__":
    main()
