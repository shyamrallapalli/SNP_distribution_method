#!/usr/bin/env ruby
require 'test/unit'
require_relative '../lib/fasta_handle'

class TestSuff < Test::Unit::TestCase
  def setup
    @fasta_file = "test/test2.fasta"
    @sequences, @total_length = FastaHandle.file_parse(@fasta_file)
  end

  def test_file_parse
    assert_equal(@sequences[:seq], {"frag1" => ">frag1 Length = 8", "CCAAATAC\n", "frag2" => ">frag2 Length = 8", "GCAATCGG\n", "frag3" => ">frag3 Length = 7", "ACGACAC\n"})
    assert_equal(@sequences[:len], {"frag1" => 8, "frag2" => 8, "frag3" => 7})
    assert_equal(@total_length, 23)
  end

  def test_write_order
    ids = ['frag1', 'frag2', 'frag3']
    filename = "test/test3.fasta"
    perm = ["frag1", "frag3", "frag2"]
    fasta_perm = FastaHandle.write_order(perm, @sequences[;seq], filename)
    seqs, total_len = FastaHandle.file_parse("test/test3.fasta")
    assert_equal(seqs[:seq], {"frag1" => ">frag1 Length = 8", "CCAAATAC\n", "frag3" => ">frag3 Length = 7", "ACGACAC\n", "frag2" => ">frag2 Length = 8", "GCAATCGG\n"})
    assert_equal(seqs[:len], {"frag1" => 8, "frag3" => 7, "frag2" => 8})
  end

end
