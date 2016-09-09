#!/usr/bin/python
#
# Analyze the .bed/.vcf files for positives and no-calls.
# Arguments: names of the .bed files. .vcf files must also be present.
#

import argparse
import os
import sys

def analyzeVcf(file):
  """Returns a dict of position -> mutation mappings"""
  with open(os.path.splitext(file)[0] + '.vcf') as vcffile:
    result = {}
    for line in vcffile:
      fields = line.split()
      if (fields[0] == 'chrY' and fields[6] == 'PASS' and fields[3] != '.'
          and fields[4] != '.'):
        result[int(fields[1])] = fields[1] + '.' + fields[3] + '.' + fields[4]
  return result


def analyzeBed(file):
  """Returns an array of path segments."""
  with open(os.path.splitext(file)[0] + '.bed') as bedfile:
    result = []
    for line in bedfile:
      fields = line.split()
      if (fields[0] == 'chrY'):
        result.append((int(fields[1]), int(fields[2])))
  return result


def makeCall(pos, index_container, bed_calls):
  call = ';nc'
  for bed_index in xrange(index_container[0], len(bed_calls)):
    pos_pair = bed_calls[bed_index]
    index_container[0] = bed_index
    if pos_pair[0] > pos:
      # We have not reached the right segment.
      continue
    if pos_pair[1] >= pos:
      # Position is within this segment.
      if pos_pair[0] == pos_pair[1] and pos_pair[0] == pos:
        call = ';cblu'
      elif pos_pair[0] == pos:
        call = ';cbl'
      elif pos_pair[1] == pos:
        call = ';cbu'
      else:
        call = ''
      return call
    # Else the position is after this segment. Continue.
  return call # No segment found.

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('files', nargs='*')
  args = parser.parse_args()

  d = []
  s = []

  with open('variant-list.txt') as line_headings:
    for line in line_headings:
      d.append(line.rstrip())
      x = line.split(',')
      s.append(int(x[0]))  # s holds the genome position for each line

  for file in args.files:
    vcf_calls = analyzeVcf(file)
    bed_calls = analyzeBed(file)
    bed_index = [0]
    for lineno in xrange(len(d)):
      d[lineno] += ','
      if s[lineno] in vcf_calls:
        d[lineno] += vcf_calls[s[lineno]]
      d[lineno] += makeCall(s[lineno], bed_index, bed_calls)

  for line in d:
    print line

if __name__ == '__main__':
  sys.exit(main())
