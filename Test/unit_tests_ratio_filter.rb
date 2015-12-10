require 'test/unit'
require_relative '../lib/ratio_filter'

class TestRatioFilter < Test::Unit::TestCase
	def test_selected_ratios
		snp_hm = [0, 14, 20, 2]
		snp_ht = [5, 4, 2, 5]
		threshold = 0
		adjust = 1
		ids = ["frag1", "frag2", "frag3", "frag4"]
		dic_ratios, ratios = RatioFilter.selected_ratios(snp_hm, snp_ht, ids, threshold, adjust)
		assert_kind_of(Hash, dic_ratios)
		assert_kind_of(Array, ratios)
		assert_equal(dic_ratios, {"frag1" => 1/6.to_f, "frag2" => 3.to_f, "frag3" => 7.to_f, "frag4" => 1/2.to_f})
		assert_equal(ratios, [1/6.to_f, 3.to_f, 7.to_f, 1/2.to_f])
	end

end
