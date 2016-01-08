#!/usr/bin/python
#
# Unit tests for the Big Y SNP manipulation module.
#

import unittest

import bigy_snp

def test_tree():
  """Returns a sample tree for testing.

  The clade tree is mut1, mut2 / mut3 and mut31 as sub of mut3
  We have persons with all sets except the ancestral."""
  personset = {}
  personset['mut1'] = {'mut1': 'mut1'}
  personset['mut2'] = {'mut1': 'mut1', 'mut2': 'mut2'}
  personset['mut3'] = {'mut1': 'mut1', 'mut3': 'mut3'}
  personset['mut31'] = {'mut1': 'mut1', 'mut3': 'mut3', 'mut31': 'mut31'}
  return personset
  

class TestBigY(unittest.TestCase):

  def test_nothing(self):
    pass

  def test_empty_set(self):
    collection = bigy_snp.KitCollection({})
    self.assertEquals(0, collection.count())
    self.assertEquals(set(), collection.snps())
    self.assertEquals(set(), collection.consistent_snps())
    self.assertEquals(set(), collection.uncertain_snps())
    self.assertEquals(set(), collection.inconsistent_snps())

  def test_basicset(self):
    collection = bigy_snp.KitCollection(test_tree())
    self.assertEquals(4, collection.count())
    # The basic set has 4 members and 4 mutations.
    # mut1 is common to all. The rest are inconsistent.
    self.assertEquals(4, len(collection.snps()))
    self.assertEquals(1, len(collection.consistent_snps()))
    self.assertEquals(3, len(collection.inconsistent_snps()))
    self.assertEquals(0, len(collection.uncertain_snps()))

  def test_uncertainty(self):
    collection = bigy_snp.KitCollection(test_tree())
    collection.add_person('uncertain', {'mut1': 'mut1', 'mut2': 'nocall'})
    self.assertEquals(5, collection.count())
    self.assertEquals(1, len(collection.uncertain_snps()))
    # The uncertain one shouldn't change any of the other counts.
    self.assertEquals(1, len(collection.consistent_snps()))
    self.assertEquals(3, len(collection.inconsistent_snps()))
    
  def test_split(self):
    collection = bigy_snp.KitCollection(test_tree())
    filtered = collection.filter('mut2')
    self.assertEquals(1, filtered.count())
    filtered = collection.filter('mut3')
    self.assertEquals(2, filtered.count())
    # A no-call person should not be included in the filtered set.
    collection.add_person('uncertain', {'mut1': 'mut1', 'mut2': 'nocall'})
    filtered = collection.filter('mut2')
    self.assertEquals(1, filtered.count())

  def test_subclade_candidates(self):
    collection = bigy_snp.KitCollection(test_tree())
    subclades = collection.subclade_candidates()
    # The result should be that we find mut2 and mut3,
    # since mut1 doesn't split the set.
    self.assertEquals(set(subclades.keys()), set(('mut2', 'mut3')))
    # Add a person without mut1. This should lead to mut1 and "mut0" being
    # found.
    collection.add_person('toplevel', {'mut0': 'mut0'})
    subclades = collection.subclade_candidates()
    self.assertEquals(set(subclades.keys()), set(('mut0', 'mut1')))
    # An empty collection has no subclades.
    collection = bigy_snp.KitCollection({})
    self.assertEquals({}, collection.subclade_candidates())
    # One person can't have subclades either.
    collection.add_person('random', {'mut0': 'mut0'})
    self.assertEquals({}, collection.subclade_candidates())

  def test_subclade_with_nc(self):
    collection = bigy_snp.KitCollection(test_tree())
    # A test person with a nocall on mut3.
    # The desired result is that because mut31 is still a possible
    # subclade of mut3, mut31 does NOT bubble up.
    collection.add_person('nocall-mut3', {'mut1': 'mut1',
                                          'mut3': 'nc', 'mut31': 'mut31'})
    subclades = collection.subclade_candidates()
    self.assertEquals(set(subclades.keys()), set(('mut2', 'mut3')))

    

if __name__ == '__main__':
  unittest.main()
