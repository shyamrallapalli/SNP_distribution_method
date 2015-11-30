#encoding: utf-8
require_relative '../lib/plot'
require_relative '../lib/write_it'
require 'test/unit'

class TestPlot < Test::Unit::TestCase

	def setup
		@ratios = WriteIt.file_to_floats_array("test/ratios_example.txt")
		@contig_size = 10
	end

	def test_ylim
		ylim = Plot.get_ylim([1,2,3,4,5], 2000)
		assert_equal(0.25, ylim)
	end

end
