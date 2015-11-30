#encoding: utf-8
require_relative '../lib/file_rw'
require 'test/unit'

class TestFileRW < Test::Unit::TestCase

	def setup
		@file = "test/ratio_values.txt"
	end

	def test_to_array
		contents = FileRW.to_array(@file)
		assert_kind_of(Array, contents)
		assert_kind_of(String, contents[0])
	end

	def test_to_array_int
		contents = FileRW.to_array_int(@file)
		assert_kind_of(Array, contents)
		assert_kind_of(Integer, contents[0])
	end

	def test_to_array_float
		contents = FileRW.to_array_float(@file)
		assert_kind_of(Array, contents)
		assert_kind_of(Float, contents[0])
	end

end

