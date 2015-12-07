#!/usr/bin/env ruby
require 'test/unit'
require_relative '../lib/fasta_handle'

class TestSuff < Test::Unit::TestCase
	def setup
		@vcf_file = "test/test.vcf"
		@fasta_file = "test/test.fasta"
		@f_array = FastaHandle.array("test/test.fasta")
	end

	def test_fasta_id_n_lengths
		ids_n_lengths = FastaHandle.fasta_id_n_lengths(@f_array)
		ids = ids_n_lengths[0]
		lengths = ids_n_lengths[1]
		assert_equal(['frag1', 'frag2', 'frag3'], ids)
		assert_equal([11,8,8], lengths)
	end

  def test_fasta_array
    fasta_array = FastaHandle.array('test/test.fasta')
    assert_equal('frag1', fasta_array[0].entry_id)
    assert_equal('AAAAAAAA', fasta_array[1].seq)
    assert_equal(8, fasta_array[2].length)
  end

	def test_genome_length
		total_length = FastaHandle.genome_length("test/test.fasta")
		assert_equal(total_length, 27)
	end

  def test_create_perm_fasta
    perm = []
    @fasta_array = FastaHandle.array("test/test2.fasta")
    ids, lengths = FastaHandle.fasta_id_n_lengths(@fasta_array)
    perm = ["frag1", "frag3", "frag2"]
    fasta_perm = FastaHandle.create_perm_fasta(perm, @fasta_array, ids)
    assert_equal(fasta_perm, [">frag1 Length = 8", "CCAAATAC\n", ">frag3 Length = 7", "ACGACAC\n", ">frag2 Length = 8", "GCAATCGG\n"])
  end

end
