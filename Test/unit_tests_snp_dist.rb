#encoding: utf-8
require_relative '../lib/snp_dist'
require 'test/unit'

class TestSNPdist < Test::Unit::TestCase

	def setup
		@ratios = [0.1, 0.1, 0.2, 0.5, 1.0, 0.2, 0.1]
		@contig_size = 10
	end

	def test_general_positions
		pos = SNPdist.general_positions(@contig_size, @ratios)
		assert_equal([10, 20, 30, 40, 50, 60, 70, 70], pos)
	end

	def test_densities_pos
		@pos = [10, 20, 30, 40, 50, 60, 70]
		densities = SNPdist.densities_pos(@ratios, @pos)
		assert_equal([10, 20, 30, 30, 40, 40, 40, 40, 40, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 60, 60, 70], densities)
	end

end