#encoding: utf-8
require_relative 'file_rw'
require_relative 'plot'
require 'simple_stats'

class Mutation

  # function to get separate array of cumulative variant positions and ratios
  # input: a hash of frag ids with all details and variant positions
  # hash input is resulted from varpos_aggregate method from Vcf class
  # output: array of homozygous (hm) and heterozygous (ht) var positions and a false density to make hm/ht ratios
  def self.putative_density(outcome)
    hm_list = []
    ht_list = []
    positions = []
    cumulative_len = 0
    outcome.each_key do | frag |
      ht_list << outcome[frag][:ht_pos]
      hmpos = outcome[frag][:hm_pos]
      hm_list << hmpos
      cumulative_len += outcome[frag][:len].to_i
      ratio = outcome[frag][:ratio]
      next if ratio < 0.1
      multiplier = (ratio * 10).to_i
      if hmpos.empty?
        positions << [cumulative_len] * multiplier
      else
        hmpos.each do | pos |
          positions << [pos] * multiplier
        end
      end
    end
    [hm_list.flatten!, ht_list.flatten!, positions.flatten!]
  end

  # selected frags with likelihood of carrying mutation
  # hash of fragments with variant positions
  # returns hash of candidate fragments with positions as values
  def self.get_candidates(frags, var_pos_hm)
    candidate_frags = {}
    frags.each do |frag|
      if var_pos_hm.has_key?(frag)
        candidate_frags[frag] = var_pos_hm[frag]
      end
    end
    candidate_frags
  end

  def self.density_plot(outcome, dir)
    # create arrays with the  SNP positions in the new ordered file.
    # generate experimental densities from ratios of the outcome order
    snps_hm, snps_ht, exp_order_density = putative_density(outcome)
    FileRW.write_txt("#{dir}/perm_hm", snps_hm)
    FileRW.write_txt("#{dir}/perm_ht", snps_ht)
    FileRW.write_txt("#{dir}/perm_ratio", exp_order_density)
    n = outcome.length

    Plot.densities(snps_hm, snps_ht, exp_order_density, dir, n)
    Plot.qqplot(snps_hm, dir, 'QQplot for hm density', 'Theoretical normal distribution', 'Hypothetical SNP density', 'hm_snps')
    Plot.qqplot(exp_order_density, dir, 'QQplot for the ratios', 'Theoretical normal distribution', 'Hypothetical ratios', 'ratios')

  end

  def self.adjusted_positions(candidate_frags, original, outcome)
    original_pos, outcome_pos = [], []
    candidate_frags.each { |frag|
      original_pos << original[frag][:hm_pos] if original.key?(frag)
      outcome_pos << outcome[frag][:hm_pos] if outcome.key?(frag)
    }
    original_pos.flatten!
    outcome_pos.flatten!
    # calculate mean of orignal and outcome position to adj position on plots
    original_mean = original_pos.mean
    outcome_mean = outcome_pos.mean
    adj_mean = original_mean - outcome_mean
    outcome_pos.map! {|x| x + adj_mean.to_i }
    [original_pos, outcome_pos]
  end

  def self.compare_density(outcome, mut_frags, mean_contig_len, genome_len, dir, original)
    # generate experimental densities from ratios of the outcome order
    exp_order_density = putative_density(outcome)
    real_order_density = putative_density(original)

    original_pos, outcome_pos = adjusted_positions(mut_frags, original, outcome)
    region = mean_contig_len * outcome.keys.length
    ylim = Plot.get_ylim(exp_order_density, region)
    Plot.comparison(real_order_density, exp_order_density, genome_len, dir, ylim, original_pos, outcome_pos)
  end

end
