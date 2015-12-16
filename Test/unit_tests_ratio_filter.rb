#encoding: utf-8
require_relative '../lib/ratio_filter'
require 'test/unit'
require 'yaml'

class TestRatioFilter < Test::Unit::TestCase

  def setup
    @frag_hash = YAML.load_file('data/frags_hash.yml')
  end

  def test_get_ratios
    ratios = RatioFilter.get_ratios(@frag_hash)
    assert_kind_of(Array, ratios)
    assert_equal([0.3333333333333333,
                  0.3333333333333333,
                  0.42857142857142855,
                  3.0,
                  3.0,
                  1.0,
                  0.38461538461538464],
                 ratios)
  end

  def test_ratio_hash
    ratios_hash = RatioFilter.ratio_hash(@frag_hash)
    assert_kind_of(Hash, ratios_hash)
    assert_equal({0.3333333333333333=>%w(frag68 frag49),
                  0.38461538461538464=>['frag1'],
                  0.42857142857142855=>['frag46'],
                  1.0=>['frag64'],
                  3.0=>%w(frag12 frag39)},
                 ratios_hash)
  end

	def test_selected_ratios
		threshold = 20
		ratios_hash = RatioFilter.selected_ratios(@frag_hash, threshold)
		assert_kind_of(Hash, ratios_hash)
		assert_equal({3.0=>%w(frag12 frag39), 1.0=>['frag64']},
                 ratios_hash)
	end

end
