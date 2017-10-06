#!/usr/bin/env python3
"""
   kin2ped.py - create a plink ped file from a plink genome file
   Copyright (C) 2016 Giulio Genovese

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Written by Giulio Genovese <giulio.genovese@gmail.com>
"""

import argparse, pandas as pd

parser = argparse.ArgumentParser(description = 'kin2ped.py: create a plink ped file from a kinship file (Apr 17th 2017)', add_help = False, usage = 'kin2ped.py [options]')
parser.add_argument('--kin', metavar = '<kinship>', type = str, required = True, help = 'Specify kinship file')
parser.add_argument('--zip', action = 'store_true', default = False, help = 'whether input kinship file is compressed [FALSE]')
parser.add_argument('--fam', metavar = '[filename]', type = str, help = 'Specify fam file')
parser.add_argument('--out', metavar = '[filename]', type = str, required = 'True', help = 'Specify output filename')
parser.add_argument('--pdf', metavar = '[filename]', type = str, help = 'Specify output filename for pdf')
parser.add_argument('--min-dup', metavar = 'FLOAT', type = float, default = 0.35, help = 'maximum kinship for siblings [0.35]')
parser.add_argument('--max-par', metavar = 'FLOAT', type = float, default = 0.01, help = 'maximum IBS0 for parent child duos [0.01]')
parser.add_argument('--min-kin', metavar = 'FLOAT', type = float, default = 0.09, help = 'minimum kinship for half siblings [0.09]')

try:
  parser.error = parser.exit
  args = parser.parse_args()
except SystemExit:
  parser.print_help()
  exit(2)

if args.zip:
  df = pd.read_csv(args.kin, delim_whitespace = True, compression = 'gzip', low_memory = False)
else:
  df = pd.read_csv(args.kin, delim_whitespace = True, low_memory = False)

if df.empty:
  exit()

fam = pd.read_csv(args.fam, delim_whitespace = True, header = None, names = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHE'], low_memory = False)
fid = {x[0]: x[1] for x in zip(fam['IID'], fam['FID'])}
sex = {x[0]: x[1] for x in zip(fam['IID'], fam['SEX'])}
phe = {x[0]: x[1] for x in zip(fam['IID'], fam['PHE'])}

dflite = df[df['Kinship'] > args.min_kin] # subset the dataframe
if dflite.empty:
  exit()

# generate a dictionary structure for each sample using the lite version of the database
parents = {x:[] for x in sex}
rels = {x:[] for x in sex}

# for each sample check whether 
for idx in dflite.index:
  if (dflite['IBS0'][idx] < args.max_par) & (dflite['Kinship'][idx] < args.min_dup):
    parents[dflite['ID1'][idx]] += [dflite['ID2'][idx]]
    parents[dflite['ID2'][idx]] += [dflite['ID1'][idx]]
  rels[dflite['ID1'][idx]] += [dflite['ID2'][idx]]
  rels[dflite['ID2'][idx]] += [dflite['ID1'][idx]]

trios = []
for proband in parents:
  # test all pairs of parents
  for pat in parents[proband]:
    for mat in parents[proband]:
      if pat != mat and not mat in rels[pat] and sex[pat] == 1 and sex[mat] == 2:
        trios += [(fid[proband], proband, pat, mat, sex[proband], phe[proband])]

if trios:
  pd.DataFrame(trios).to_csv(args.out, header = False, index = False, sep = '\t')

if args.pdf:
  from matplotlib.backends.backend_pdf import PdfPages
  with PdfPages(args.pdf) as pdf:
    ax = df.plot(kind='scatter', x='IBS0', y='Kinship', xlim=(0,.006), ylim=(0,1))
    ax.set_xlabel('IBS0')
    ax.set_ylabel('Kinship')
    ax.axhline(args.min_dup, color='r', linestyle='--', lw=2)
    ax.axhline(args.min_kin, color='r', linestyle='--', lw=2)
    ax.axvline(args.max_par, color='r', linestyle='--', lw=2)
    pdf.savefig(ax.get_figure())
    
