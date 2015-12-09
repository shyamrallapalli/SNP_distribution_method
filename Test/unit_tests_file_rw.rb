#encoding: utf-8
require_relative '../lib/file_rw'
require 'test/unit'

class TestFileRW < Test::Unit::TestCase

  def setup
    @file = 'test/ratio_values.txt'
    @vcf_file = 'test/test.vcf'
    @fasta_file = 'test/test.fasta'
    infasta = "test/test2.fasta"
    @sequences, @total_length = FileRW.fasta_parse(infasta)
  end

  def test_fasta_parse
    assert_equal(@sequences[:seq], {"frag1" => ">frag1 Length = 8\nCCAAATAC\n", "frag2" => ">frag2 Length = 8\nGCAATCGG\n", "frag3" => ">frag3 Length = 7\nACGACAC\n"})
    assert_equal(@sequences[:len], {"frag1" => 8, "frag2" => 8, "frag3" => 7})
    assert_equal(@total_length, 23)
  end

  def test_write_order
    ids = ['frag1', 'frag2', 'frag3']
    filename = "test/test3.fasta"
    perm = ["frag1", "frag3", "frag2"]
    fasta_perm = FileRW.write_order(perm, @sequences[:seq], filename)
    seqs, total_len = FileRW.fasta_parse("test/test3.fasta")
    assert_equal(seqs[:seq], {"frag1" => ">frag1 Length = 8\nCCAAATAC\n", "frag3" => ">frag3 Length = 7\nACGACAC\n", "frag2" => ">frag2 Length = 8\nGCAATCGG\n"})
    assert_equal(seqs[:len], {"frag1" => 8, "frag3" => 7, "frag2" => 8})
  end

end
