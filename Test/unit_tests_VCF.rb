#encoding: utf-8
require_relative '../lib/vcf'
require_relative '../lib/write_it'
require 'test/unit'

class TestVCF < Test::Unit::TestCase
	def setup
		@vcf_ngs = "test/ngs.vcf"
		@chromosome = 1
		@vcfs_info = {"ADP"=>"17", "WT"=>"0", "HET"=>"1", "HOM"=>"0", "NC"=>"0"}, {"ADP"=>"25", "WT"=>"0", "HET"=>"1", "HOM"=>"0", "NC"=>"0"}
		@vcfs_pos = [5, 123]
		@snps = {5 => "HET", 123 => "HET"}
		@vcf = ["1\t5\t.\tC\tA\t.\tPASS\tADP=17;WT=0;HET=1;HOM=0;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t0/1:24:17:17:9:7:41.18%:3.3988E-3:65:52:9:0:1:6\n", 
		"1\t123\t.\tG\tA\t.\tPASS\tADP=25;WT=0;HET=1;HOM=0;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t0/1:39:25:25:14:11:44%:1.1933E-4:69:66:6:8:6:5\n"]
    @vcf_file = "test/test.vcf"
    @fasta_file = "test/test.fasta"
    @f_array = FastaHandle.fasta_array("test/test.fasta")
	end

  def test_snps_in_vcf
    snp_data, hm, ht = Stuff.snps_in_vcf(@vcf_file)
    assert_equal(["frag1", "frag1"], hm)
    assert_equal(["frag2", "frag3"], ht)
    assert_equal([["frag1", "frag1", "frag2", "frag3"], [7, 8, 2, 2], {"frag1" =>2, "frag2" =>1, "frag3" =>1}, [{"AF"=>"1.0"},{"AF"=>"1.0"},{"AF"=>"0.5"},{"AF"=>"0.5"}]], snp_data)
  end

  def test_dic_id_pos
    hm = ["frag1", "frag1", "frag1", "frag4"]
    ht = ["frag1", "frag1", "frag2", "frag4"]
    pos1 = [12, 13, 14, 45]
    pos2 = [1, 4, 25, 40]
    dic_pos_hm = Stuff.dic_id_pos(hm, pos1)
    dic_pos_ht = Stuff.dic_id_pos(ht, pos2)
    assert_kind_of(Hash, dic_pos_hm)
    assert_kind_of(Hash, dic_pos_ht)
    assert_equal(dic_pos_hm, {"frag1"=>[12, 13, 14], "frag4"=>[45]})
    assert_equal(dic_pos_ht, {"frag1"=> [1, 4], "frag2"=>[25], "frag4"=>[40]})
  end

  def test_define_snps
    ids = ["frag1", "frag3", "frag2"]
    dic1 = {"frag1"=>"2"}
    shu_dic1, snps_1  = Stuff.define_snps(ids, dic1)
    assert_kind_of(Hash, shu_dic1)
    assert_kind_of(Array, snps_1)
    assert_kind_of(Float, snps_1[0])
    assert_equal(shu_dic1, {"frag1"=>2.0, "frag3"=>0, "frag2"=>0})
    assert_equal(snps_1, [2.0, 0.0, 0.0])
  end

  def test_positions_by_fragment
    dic = {"frag1" => 1.0, "frag2"=>1.0, "frag3"=>2.0, "frag4"=>0.0}
    snp_list = [15, 18, 20, 25]
    assert_kind_of(Hash, dic)
    assert_kind_of(Array, snp_list)
    dic = Stuff.positions_by_fragment(dic, snp_list)
    assert_equal(dic, {"frag1" =>[15], "frag2"=>[18], "frag3"=>[20, 25]})
  end

	def test_open_vcf
		vcf, chrom, pos, info = Vcf.open_vcf(@vcf_ngs, @chromosome)
		assert_equal( ["1\t5\t.\tC\tA\t.\tPASS\tADP=17;WT=0;HET=1;HOM=0;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t0/1:24:17:17:9:7:41.18%:3.3988E-3:65:52:9:0:1:6\n", 
		"1\t123\t.\tG\tA\t.\tPASS\tADP=25;WT=0;HET=1;HOM=0;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t0/1:39:25:25:14:11:44%:1.1933E-4:69:66:6:8:6:5\n"], vcf)
		assert_equal(["1", "1"], chrom)
		assert_equal([5, 123], pos)
		assert_equal( [{"ADP"=>"17", "WT"=>"0", "HET"=>"1", "HOM"=>"0", "NC"=>"0"}, {"ADP"=>"25", "WT"=>"0", "HET"=>"1", "HOM"=>"0", "NC"=>"0"}], info)
	end

	def test_type_per_pos
		snps, hm, ht = Vcf.type_per_pos(@vcfs_info, @vcfs_pos)
		assert_equal({5 => "HET", 123 => "HET"}, snps)
		assert_equal([], hm)
		assert_equal([5, 123], ht)
	end

	def test_filtering
		snps_p = {5 => "HET", 365 => "HOM"}
		short_vcf = Vcf.filtering(@vcfs_pos, snps_p, @snps, @vcf)
		assert_equal(["1\t123\t.\tG\tA\t.\tPASS\tADP=25;WT=0;HET=1;HOM=0;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t0/1:39:25:25:14:11:44%:1.1933E-4:69:66:6:8:6:5\n"], short_vcf)
	end
end