#encoding: utf-8
require_relative '../lib/mutation'
require_relative '../lib/file_rw'
require 'test/unit'
require 'yaml'

class TestMutation < Test::Unit::TestCase

	def setup
    @mean_contig_len = 12000 # 12 kb
    @ratios = [0.1, 0.1, 0.2, 0.5, 1.0, 0.2, 0.1]
    @putate_dens = [12000, 24000, 36000, 36000, 48000, 48000, 48000, 48000, 48000,\
    60000, 60000, 60000, 60000, 60000, 60000, 60000, 60000, 60000, 60000,\
    72000, 72000, 84000]
	end

  def test_putative_density
    density = Mutation.putative_density(@mean_contig_len, @ratios)
    assert_equal(density, @putate_dens)
  end

  def test_get_candidates
    candidates = {"frag39"=>[107273], "frag41"=>[87355], "frag37"=>[99032], "frag64"=>[118176], "frag46"=>[48121], "frag47"=>[18333], "frag26"=>[68402], "frag69"=>[1006], "frag34"=>[31962], "frag9"=>[94160], "frag21"=>[29408], "frag10"=>[112102], "frag12"=>[139168]}
    frags = FileRW.to_array('data/frags.txt')
    var_pos_hm = YAML.load_file('data/var_pos_hm.yml')
    candi_frags = Mutation.get_candidates(frags, var_pos_hm)
    assert_equal(candi_frags, candidates)
  end

  def test_adjusted_positions
    candi_frags = %w(frag12 frag39 frag34 frag9 frag26 frag21 frag10 frag41 frag47 frag69 frag37 frag64)
    org_pos = [1775592, 5755261, 4919240, 1247713, 3775054, 3057628, 1405085, 6035523, 6983113, 10163762, 5482374, 9507519]
    out_pos = [4239973, 4400162, 4498802, 4697824, 4811496, 4877458, 5109065, 5242169, 5329848, 5432022, 5651540, 5817512]
    outcome = YAML.load_file('data/outcome.yml')
    original = YAML.load_file('data/original.yml')
    original_pos, outcome_pos = Mutation.adjusted_positions(candi_frags, original, outcome)
    assert_equal(original_pos, org_pos)
    assert_equal(outcome_pos, out_pos)
  end

end
