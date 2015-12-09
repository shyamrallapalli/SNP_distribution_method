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

  # def test_safe_invert
  #   hash = {'frag1' => 1, 'frag2' => 1, 'frag3' => 1, 'frag4' => 2, 'frag5' => 1}
  #   hash_inv = FileRW.safe_invert(hash)
  #   assert_equal({1 => ['frag1',  'frag2', 'frag3', 'frag5'], 2 => ['frag4']}, hash_inv)
  # end

  def test_create_hash_number
    array = ['frag1', 'frag1', 'frag2', 'frag3']
    hash = FileRW.create_hash_number(array)
    assert_kind_of(Hash, hash)
    assert_equal(hash, {'frag1'=>'2', 'frag2'=>'1', 'frag3'=>'1'})
  end

  def test_fasta_parse
    assert_equal(@sequences[:seq], {"frag1" => ">frag1 Length = 8", "CCAAATAC\n", "frag2" => ">frag2 Length = 8", "GCAATCGG\n", "frag3" => ">frag3 Length = 7", "ACGACAC\n"})
    assert_equal(@sequences[:len], {"frag1" => 8, "frag2" => 8, "frag3" => 7})
    assert_equal(@total_length, 23)
  end

  def test_write_order
    ids = ['frag1', 'frag2', 'frag3']
    filename = "test/test3.fasta"
    perm = ["frag1", "frag3", "frag2"]
    fasta_perm = FileRW.write_order(perm, @sequences[:seq], filename)
    seqs, total_len = FileRW.fasta_parse("test/test3.fasta")
    assert_equal(seqs[:seq], {"frag1" => ">frag1 Length = 8", "CCAAATAC\n", "frag3" => ">frag3 Length = 7", "ACGACAC\n", "frag2" => ">frag2 Length = 8", "GCAATCGG\n"})
    assert_equal(seqs[:len], {"frag1" => 8, "frag3" => 7, "frag2" => 8})
  end

end
