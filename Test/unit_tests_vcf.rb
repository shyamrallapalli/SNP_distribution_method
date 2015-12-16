#encoding: utf-8
require_relative '../lib/vcf'
require_relative '../lib/file_rw'
require 'test/unit'
require 'bio-samtools'
require 'yaml'

class TestVcf < Test::Unit::TestCase

	def setup
		@vcf_ngs = 'data/ngs.vcf'
    @vcf = []
    File.open(@vcf_ngs, 'r').each do |line|
      next if line =~ /^#/
      @vcf << Bio::DB::Vcf.new(line)
    end
    @vars_pos = Vcf.get_vars(@vcf_ngs)
    @fasta_file = 'data/test.fasta'
	end

	def test_get_allele_freq
    vcf_obj = @vcf[0]
		allele_freq = Vcf.get_allele_freq(vcf_obj)
    assert_equal(0.4375, allele_freq)
    vcf_obj = @vcf[1]
    allele_freq = Vcf.get_allele_freq(vcf_obj)
    assert_equal(0.44, allele_freq)
  end

  def test_get_vars
    assert_kind_of(Hash, @vars_pos)
    assert_equal({:het=>{'frag1'=>[5], 'frag3'=>[3]}}, @vars_pos)
  end

  def test_filtering
    vcf_bg = 'data/ngs_bg.vcf'
    vars_pos = Vcf.filtering(@vcf_ngs, vcf_bg)
    assert_kind_of(Hash, vars_pos)
    assert_equal({:het=>{'frag1'=>[5]}}, vars_pos)
  end

  def test_varpos_aggregate
    seqs = FileRW.fasta_parse(@fasta_file)
    details = Vcf.varpos_aggregate(@vars_pos, seqs[:len], %w(frag3 frag1), 0.5)
    assert_kind_of(Hash, details)
    assert_equal({'frag3'=>{:hm=>0.5, :ht=>1.5, :hm_pos=>[], :ht_pos=>[3], :ratio=>0.3333333333333333, :len=>8},
                  'frag1'=>{:hm=>0.5, :ht=>1.5, :hm_pos=>[], :ht_pos=>[13], :ratio=>0.3333333333333333, :len=>11}},
                 details)
  end

  def test_varpositions
    outcome = YAML.load_file('data/outcome.yml')
    hm, ht = Vcf.varpositions(outcome)
    assert_kind_of(Array, hm)
    assert_kind_of(Array, ht)
    assert_equal([11057806, 11271450, 11431639, 11530279, 11729301, 11842973, 11908935, 12140542,
                  12273646, 12361325, 12463499, 12683017, 12848989, 12869226, 12951951], hm)
  end

end