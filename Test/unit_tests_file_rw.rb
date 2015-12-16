#encoding: utf-8
require_relative '../lib/file_rw'
require 'test/unit'

class TestFileRW < Test::Unit::TestCase

  def setup
    @file = 'data/fragments.txt'
    @fragments = %w(frag1 frag2 frag3)
    infasta = 'data/test2.fasta'
    @sequences = FileRW.fasta_parse(infasta)
  end

  def test_to_array
    fragments = FileRW.to_array(@file)
    assert_equal(fragments, @fragments)
  end

  def test_write
    filename = 'data/array_file.txt'
    FileRW.write(filename, @fragments)
    frags = FileRW.to_array(filename)
    assert_equal(frags, @fragments)
    %x[rm data/array_file.txt]
  end

  def test_write_txt
    filename = 'data/array_file'
    FileRW.write_txt(filename, @fragments)
    frags = FileRW.to_array("#{filename}.txt")
    assert_equal(frags, @fragments)
    %x[rm "#{filename}.txt"]
  end

  def test_fasta_parse
    assert_equal(@sequences[:seq], {"frag1"=>">frag1 Length = 8\nCCAAATAC\n", "frag2"=>">frag2 Length = 8\nGCAATCGG\n", "frag3"=>">frag3 Length = 7\nACGACAC\n"})
    assert_equal(@sequences[:len], {"frag1"=>8, "frag2"=>8, "frag3"=>7})
  end

  def test_write_order
    filename = 'data/test3.fasta'
    perm = %w(frag1 frag3 frag2)
    FileRW.write_order(perm, @sequences[:seq], filename)
    seqs = FileRW.fasta_parse('data/test3.fasta')
    assert_equal(seqs[:seq], {"frag1"=>">frag1 Length = 8\nCCAAATAC\n", "frag3"=>">frag3 Length = 7\nACGACAC\n", "frag2"=>">frag2 Length = 8\nGCAATCGG\n"})
    assert_equal(seqs[:len], {"frag1"=>8, "frag3"=>7, "frag2"=>8})
    %x[rm data/test3.fasta]
  end

end



