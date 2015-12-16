#encoding: utf-8
require_relative 'file_rw'
require_relative 'plot'
require 'simple_stats'

class Mutation

  # Input 0: The average contig length
  # Input 1: Array of homozygous (hm) or heterozygous (hm) or ratio (hm/ht) of SNPs per each contig.
  # Output: A list of positions which represents the density of hm or ht or hm/ht ratios
  def self.putative_density(mean_contig_len, ratios)
    positions = []
    ratios.length.times do | i |
      multiple = i + 1
      positions << [mean_contig_len * multiple] * (ratios[i] * 10).to_i
    end
    positions.flatten!
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

  def self.density_plot(outcome, mean_contig_len, dir)
    # create arrays with the  SNP positions in the new ordered file.
    snps_hm, snps_ht = Vcf.varpositions(outcome)
    FileRW.write_txt("#{dir}/perm_hm", snps_hm)
    FileRW.write_txt("#{dir}/perm_ht", snps_ht)

    # generate experimental densities from ratios of the outcome order
    region = mean_contig_len * outcome.keys.length
    outcome_ratios = RatioFilter.get_ratios(outcome)
    exp_order_density = putative_density(mean_contig_len, outcome_ratios)

    # Find the peak in the approximated (hypothetical SNP) distribution
    # candidate_snp = closest_snp(exp_order_density, snps_hm, n)

    Plot.densities(snps_hm, snps_ht, exp_order_density, region, dir)
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
    outcome_ratios = RatioFilter.get_ratios(outcome)
    exp_order_density = putative_density(mean_contig_len, outcome_ratios)
    original_ratios = RatioFilter.get_ratios(original)
    real_order_density = putative_density(mean_contig_len, original_ratios)

    original_pos, outcome_pos = adjusted_positions(mut_frags, original, outcome)
    region = mean_contig_len * outcome.keys.length
    ylim = Plot.get_ylim(exp_order_density, region)
    Plot.comparison(real_order_density, exp_order_density, genome_len, dir, ylim, original_pos, outcome_pos)
  end

end
