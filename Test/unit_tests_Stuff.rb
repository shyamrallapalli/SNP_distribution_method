#!/usr/bin/env ruby
require 'test/unit'
# require_relative '../lib/reform_ratio'
require_relative '../lib/stuff'

class TestSuff < Test::Unit::TestCase
	def setup
		@vcf_file = "test/test.vcf"
		@fasta_file = "test/test.fasta"
		@f_array = FastaHandle.fasta_array("test/test.fasta")
	end

	def test_safe_invert
		hash = {"frag1" => 1, "frag2" => 1, "frag3" => 1, "frag4" => 2, "frag5" => 1}
		hash_inv = Stuff.safe_invert(hash)
		assert_equal({1 => ["frag1",  "frag2", "frag3", "frag5"], 2 => ["frag4"]}, hash_inv)
	end

	def test_create_hash_number
		hm = ["frag1", "frag1"]
		ht = ["frag2", "frag3"]
		dic_hm = Stuff.create_hash_number(hm)
		dic_ht = Stuff.create_hash_number(ht)
		assert_kind_of(Hash, dic_hm)
		assert_kind_of(Hash, dic_ht)
		assert_equal(dic_hm, {"frag1"=>"2"})
		assert_equal(dic_ht, {"frag2"=>"1", "frag3"=>"1"})
	end

end

