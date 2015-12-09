#encoding: utf-8
require_relative '../lib/mutation'
require 'test/unit'

class TestLocateMutation < Test::Unit::TestCase

	def setup
		@example_snps = [105,109,87,96,111,95,100,88,110,92]
    @ratios = [0.1, 0.1, 0.2, 0.5, 1.0, 0.2, 0.1]
    @contig_size = 10
	end

  def test_general_positions
    pos = Mutation.general_positions(@contig_size, @ratios)
    assert_equal([10, 20, 30, 40, 50, 60, 70, 70], pos)
  end

  def test_densities_pos
    @pos = [10, 20, 30, 40, 50, 60, 70]
    densities = Mutation.densities_pos(@ratios, @pos)
    assert_equal([10, 20, 30, 30, 40, 40, 40, 40, 40, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 60, 60, 70], densities)
  end

	def test_closest
		assert_equal(95, Mutation.closest_snp(@example_snps, @example_snps, 512))
	end
end