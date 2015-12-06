#encoding: utf-8
require_relative '../lib/file_rw'
require 'test/unit'

class TestFileRW < Test::Unit::TestCase

  def setup
    @file = 'test/ratio_values.txt'
    @vcf_file = 'test/test.vcf'
    @fasta_file = 'test/test.fasta'
    @f_array = FastaHandle.fasta_array('test/test.fasta')
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

  def test_safe_invert
    hash = {'frag1' => 1, 'frag2' => 1, 'frag3' => 1, 'frag4' => 2, 'frag5' => 1}
    hash_inv = FileRW.safe_invert(hash)
    assert_equal({1 => ['frag1',  'frag2', 'frag3', 'frag5'], 2 => ['frag4']}, hash_inv)
  end

  def test_create_hash_number
    hm = ['frag1', 'frag1']
    ht = ['frag2', 'frag3']
    dic_hm = FileRW.create_hash_number(hm)
    dic_ht = FileRW.create_hash_number(ht)
    assert_kind_of(Hash, dic_hm)
    assert_kind_of(Hash, dic_ht)
    assert_equal(dic_hm, {'frag1'=>'2'})
    assert_equal(dic_ht, {'frag2'=>'1', 'frag3'=>'1'})
  end

end
