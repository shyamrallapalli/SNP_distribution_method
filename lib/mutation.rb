#encoding: utf-8
require_relative 'file_rw'
require_relative 'plot'
require 'rinruby'
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

  # Input 0: List ratios from sdm
  # Input 1: List of homozygous SNP positions
  # Input 2: The number of equally spaced points at which the density is to be estimated. Specify n as a power of two.
  # Output: get the highest kernel density value and returns closest homozygous SNP to the peak
  def self.closest_snp(ratios, hm, n)
    myr = RinRuby.new(:echo=>false)
    myr.n = n
    myr.ratios = ratios
    myr.eval 'kernel_density <- density(ratios, n=n)'
    # this only finds the first index with the max density if there is > 1
    myr.eval 'index <- match(max(kernel_density$y),kernel_density$y)'
    myr.eval 'peak <- kernel_density$x[index]'
    peak = myr.pull 'peak'
    myr.quit
    peak.to_i
    hm.min { |a, b| (peak - a).abs <=> (peak - b).abs }
  end

  # selected frags with likelyhood of carrying mutation
  # hash of fragments with variant positions
  # returns hash of 
  def self.get_candidates(frags, var_pos_hm)
    candidate_frags = {}
    frags.each do |frag|
      if var_pos_hm.has_key?(frag)
        candidate_frags[frag] = var_pos_hm[frag]
      end
    end
    candidate_frags
  end

  def self.adjusted_positions(candidate_frags, original, outcome)
    original_pos, outcome_pos = [], []
    candidate_frags.each_key { |frag|
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

  def self.density_plots(outcome, mut_frags, frag_pos_hm, mean_contig_len, genome_len, dir, original, logdir)
    n = 1048576*4
    # create arrays with the  SNP positions in the new ordered file.
    snps_hm, snps_ht = Vcf.varpositions(outcome)
    FileRW.write_txt("#{dir}/perm_hm", snps_hm)
    FileRW.write_txt("#{dir}/perm_ht", snps_ht)

    # generate experimental densities from ratios of the outcome order
    region = mean_contig_len * outcome.keys.length
    outcome_ratios = Ratio_filtering.get_ratios(outcome)
    FileRW.write_txt("#{logdir}/5_2_expected_ratios", outcome_ratios)
    exp_order_density = putative_density(mean_contig_len, outcome_ratios)

    # Find the peak in the approximated (hypothetical SNP) distribution
    candidate_snp = closest_snp(exp_order_density, snps_hm, n)
    candidate_frags = get_candidates(mut_frags, frag_pos_hm)
    original_pos, outcome_pos = Mutation.adjusted_positions(candidate_frags, original, outcome)

    original_ratios = Ratio_filtering.get_ratios(original)
    FileRW.write_txt("#{logdir}/t_08_ratios", original_ratios)
    real_order_density = putative_density(mean_contig_len, original_ratios)
    ylim = Plot.get_ylim(exp_order_density, region)
    Plot.densities(snps_hm, snps_ht, exp_order_density, region, dir)
    Plot.comparison(real_order_density, exp_order_density, genome_len, dir, ylim, original_pos, outcome_pos)
    Plot.qqplot(snps_hm, dir, 'QQplot for hm density', 'Theoretical normal distribution', 'Hypothetical SNP density', 'hm_snps')
    Plot.qqplot(exp_order_density, dir, 'QQplot for the ratios', 'Theoretical normal distribution', 'Hypothetical ratios', 'ratios')

    # write the candidate mutation and related ratio information to files
    FileRW.write_txt("#{dir}/hyp_ratios", exp_order_density)
    FileRW.write_txt("#{dir}/ratios", real_order_density)
    File.open("#{dir}/mutation.txt", 'w+') do |f|
      f.puts "The length of the group of contigs that form the peak of the distribution is #{region.to_i} bp"
      f.puts "The mutation is likely to be found on the following contigs #{candidate_frags}"
      f.puts "Likeliest candidate responsible for mutation from the mutations is at #{candidate_snp}"
    end
  end

end
