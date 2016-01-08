#
# Classes and objects for handling SNPs in the BigY spreadsheets
#

import collections
import csv
import re

class KitCollection(object):
  def __init__(self, persondict):
    self.persons = persondict

  def count(self):
    return len(self.persons)

  def add_person(self, id, mutations):
    self.persons[id] = mutations

  def split(self, snp):
    selected = {}
    not_selected = {}
    uncertain = {}
    for person in self.persons.keys():
      if snp in self.persons[person]:
        if self.persons[person][snp] == snp:
          selected[person] = self.persons[person]
        else:
          uncertain[person] = self.persons[person]
      else:
        not_selected[person] = self.persons[person]
    return (KitCollection(selected),
            KitCollection(not_selected),
            KitCollection(uncertain))


  def filter(self, snp):
    return self.split(snp)[0]

  def snps(self):
    found = set()
    for person in self.persons.values():
      for snp in person.keys():
        found.add(snp)
    return found

  def consistent_snps(self):
    """SNPs that have the same value over the whole set of persons."""
    candidates = self.snps()
    for person in self.persons.values():
      candidates &= set([snp for snp in person.keys() if person[snp] == snp])
    return candidates

  def uncertain_snps(self):
    """SNPs that have some people with a "no call" or "cb" marking."""
    candidates = set()
    for person in self.persons.values():
      for snp in person:
        if person[snp] != snp:
          #if not snp in candidates:
          #  print 'Uncertainty added:', snp, person[snp]
          candidates.add(snp)
    return candidates

  def inconsistent_snps(self):
    """SNPs that vary over this set."""
    candidates = self.snps() - self.consistent_snps()
    verified_inconsistent = set()
    for person in self.persons.values():
      # If we know an SNP is inconsistent, don't look at it again.
      candidates.difference_update(verified_inconsistent)
      for snp in candidates:
        if not snp in person:
          # This means a negative call on the SNP, not just a no-call.
          verified_inconsistent.add(snp)
    return verified_inconsistent
    # Old - this is still wrong for the case of conflicts despite NCs.
    #return self.snps() - self.consistent_snps() - self.uncertain_snps()

  def subclade_candidates(self):
    """SNPs that may be first-level subclades.

    Returns a dict of firstlevel: set(equivalents).
    """
    consistent_for_this = {}
    for snp in self.inconsistent_snps():
      next_level = self.filter(snp)
      if next_level.count() > 0:
        consistent_for_this[snp] = next_level.consistent_snps()
    equivalents = collections.defaultdict(set)

    first_level_candidates = set(consistent_for_this.keys())
    second_look = set()
    for snp1 in consistent_for_this:
      for snp2 in consistent_for_this:
        if snp1 == snp2:
          continue
        if consistent_for_this[snp1] == consistent_for_this[snp2]:
          equivalents[snp1].add(snp2)
        elif consistent_for_this[snp1] > consistent_for_this[snp2]:
          # SNP1 is a subclade of SNP2: All consistent in snp2 are also
          # consistent in snp1, and then some.
          first_level_candidates.discard(snp1)
        elif consistent_for_this[snp1] < consistent_for_this[snp2]:
          # print snp1, snp2, 'subset'
          first_level_candidates.discard(snp2)
        else:
          # Inconsistency, we can't decide anything from this.
          #print snp1, snp2, 'inconsistent'
          second_look.add(snp1)
          second_look.add(snp2)
    # Take a second look at inconsistent pairs that remain at top level.
    # Here we look at people, not SNPs.
    second_look = second_look & first_level_candidates
    if second_look:
      for snp1 in second_look:
        for snp2 in second_look:
          if (not consistent_for_this[snp1] < consistent_for_this[snp2]
              and not consistent_for_this[snp1] > consistent_for_this[snp2]
              and not consistent_for_this[snp1] == consistent_for_this[snp2]):
            (snp1_yes, snp1_no, snp1_maybe) = self.split(snp1)
            (snp2_yes, snp2_no, snp2_maybe) = self.split(snp2)
            # If all the certain ones in snp1 are true or uncertain in snp2,
            # and there are some snp2 (real or uncertain) not in snp1,
            # we think snp1 is a subset of snp2, and drop it from consideration.
            # Note that potential equivalents are missed by this method.
            if set(snp1_yes.persons.keys()) < set(snp2_yes.persons.keys()
                                                  + snp2_maybe.persons.keys()):
              first_level_candidates.discard(snp1)

    # Remove those in equivalent sets
    do_not_present = set()
    for snp in sorted(first_level_candidates):
      if not snp in do_not_present:
        do_not_present |= equivalents[snp]
    result = {}
    for snp in (first_level_candidates - do_not_present):
      result[snp] = equivalents[snp]
    return result

def load_spreadsheet(filename):
  """Load the SNP spreadsheet.

  The input spreadsheet has 1 column per participant and 1 row per mutation.
  The first few rows are header, with the very first line having kitnames.
  """
  with open(filename) as csvfile:
    reader = csv.reader(csvfile)

    kitnames = reader.next()
    new_kitnames = []
    for kitname in kitnames:
      matches = re.match(r'^\S+ (\S+)\.zip$', kitname)
      if matches:
        kitname = matches.group(1)
      new_kitnames.append(kitname)
    kitnames = new_kitnames
          
    kits = collections.defaultdict(dict)
    # Skip header
    for row in reader:
      if row[0] == 'SNP Number':
        break

    # Include only the part down to the first comment block
    for row in reader:
      if len(row) < 2:
        break
      position = int(row[0])
      if row[1]:
        name = row[1].split(' ')[0]
      else:
        name = position
      for pos in xrange(4,len(row)):
        if row[pos]:
          #print 'set', pos, kitnames[pos], name
          kits[kitnames[pos]][name] = row[pos]

  return KitCollection(kits)
