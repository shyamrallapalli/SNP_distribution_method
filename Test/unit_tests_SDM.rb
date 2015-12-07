#!/usr/bin/env ruby
require 'test/unit'
require_relative '../lib/SDM'

class TestSDM < Test::Unit::TestCase
  def setup
    @dic_ratio_inv = {11.374979 => ['a'], 11.842904 => ['b'], 10.120768 => ['c'], 5.022447 => ['d', 'e', 'f', 'g', 'h', 'i', 'j']}
    @dic_hm_inv = {11.374979 => ['a'], 11.842904 => ['b'], 10.120768 => ['c'], 5.022447 => ['d', 'e', 'f', 'g', 'h', 'i', 'j']}
    @dic_pos_hm = {'a' => [1], 'b' => [2], 'c' => [3, 4], 'd' => [5, 6], 'e' => [7], 'f' => [10], }
  end

  def test_split
    inkeys = @dic_hm_inv.keys
    right, left, inkeys = SDM.split(@dic_hm_inv, right = [], left = [], inkeys, 0)
    assert_equal( ['e', 'f', 'g', 'd'], right)
    assert_equal( ['h', 'i', 'j'], left)
    assert_equal( [11.374979, 11.842904, 10.120768], inkeys)
  end

  def test_sort
    perm_back, mut_back = SDM.sort(@dic_hm_inv, cross = 'back', average_contig = 2500)
    assert_equal(['e', 'd', 'a', 'b', 'c', 'f'], perm_back)
    assert_equal(['e', 'd', 'a', 'b', 'c', 'f'], mut_back)
  end

  def test_arrange
    perm_hm, perm_ratio, mut, hyp = SDM.arrange(@dic_hm_inv, @dic_ratio_inv, @dic_pos_hm, cross = 'back', average_contig = 2500)
    assert_equal(['e', 'd', 'a', 'b', 'c', 'f'], perm_hm)
    assert_equal(['e', 'd', 'a', 'b', 'c', 'f'], perm_ratio)
    assert_equal(['e', 'd', 'a', 'b', 'c', 'f'], mut)
    assert_equal([1, 2, 3, 4, 5, 6, 7, 10], hyp)
  end

end
