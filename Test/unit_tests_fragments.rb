#!/usr/bin/env ruby
require 'test/unit'
require_relative '../lib/fragments'

class TestFragments < Test::Unit::TestCase
  def setup
    @dic_ratios = {11.374979=>['a'], 11.842904=>['b'], 10.120768=>['c'], 5.022447=>%w(d e f g h i j)}
  end

  def test_split
    left = []
    right = []
    inkeys = @dic_ratios.keys
    left, right, inkeys = Fragments.split(@dic_ratios, left, right, inkeys, 0)
    assert_equal(['d', %w(e f g)], left)
    assert_equal([%w(h i j)], right)
    assert_equal([11.374979, 11.842904, 10.120768], inkeys)
  end

  def test_arrange
    perm_ratio, mut = Fragments.arrange(@dic_ratios, cross = 'back', average_contig = 2500)
    assert_equal(%w(d e f g a b c j i h), perm_ratio)
    assert_equal(%w(d e f g a b c j i h), mut)
  end

end
