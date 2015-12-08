#encoding: utf-8
require_relative 'file_rw'
require_relative 'plot'
require_relative 'snp_dist'
require 'rinruby'
require 'simple_stats'

class Mutation

  # Input 0: List of SNP positions
  # Input 1: The number of equally spaced points at which the density is to be estimated. Specify n as a power of two.
  # Output: The highest kernel density value for this SNP distribution
  def self.find_peak(snps, n)
    myr = RinRuby.new(:echo=>false)
    myr.n = n
    myr.snps = snps
    myr.eval 'kernel_density <- density(snps, n=n)'
    myr.eval 'index <- match(max(kernel_density$y),kernel_density$y)' # this only finds the first index with the max density if there is > 1
    myr.eval 'peak <- kernel_density$x[index]'
    peak = myr.pull 'peak'
    myr.quit
    peak.to_i
  end

  # Input 0: Value under distribution peak (genome position as a float)
  # Input 1: List of homozygous SNP positions
  # Output: The closest homozygous SNP to the peak
  def self.closest_snp(peak, hm)
    hm.min { |a, b| (peak - a).abs <=> (peak - b).abs }
  end

	def self.test_genomes_define(hm, ht, perm_hm, perm_ht, genome_length, ratios, expected_ratios)
		hm.pop
		# div = 100
		n = 2000
		hyp = SNPdist.hyp_snps(ratios, genome_length)
		peak =  find_peak(hyp, n) # Find the peak in the approximated (hypothetical SNP) distribution
		causal = closest_snp(peak, hm)
		perm_hyp = SNPdist.hyp_snps(expected_ratios, genome_length)

		perm_peak = find_peak(perm_hyp, n)
		candidate = closest_snp(perm_peak, perm_hm)
		normalised = (candidate - causal).abs
		percent = (normalised*100)/genome_length.to_f
		return causal, candidate, percent
	end

	def self.candidate(mut, frag_pos_hm)
		candidate_mutations ={}
		mut.each do |frag|
			if frag_pos_hm.has_key?(frag)
				candidate_mutations.store(frag, frag_pos_hm[frag])
			end
		end
		return candidate_mutations
	end

  def self.adjusted_positions(candidate_mutations, original, outcome)
    original_pos, outcome_pos = [], []
    candidate_mutations.each_key { |frag|
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

  def self.density_plots(contig_size, ratios, expected_ratios, snps_hm, snps_ht, region, genome_len, file, mut, frag_pos_hm, original, outcome)
		n = 1048576*4
		average_positions = SNPdist.general_positions(contig_size, ratios)
		hyp_ratios = SNPdist.densities_pos(expected_ratios, average_positions)
		real_ratios = SNPdist.densities_pos(ratios, average_positions)

  	# Find the peak in the approximated (hypothetical SNP) distribution
    peak =  find_peak(hyp_ratios, n)
		ylim = Plot.get_ylim(hyp_ratios, region)
		candidate_peak = closest_snp(peak, snps_hm)
		candidate_mutations = Mutation.candidate(mut, frag_pos_hm)
    original_pos, outcome_pos = Mutation.adjusted_positions(candidate_mutations, original, outcome)
		Plot.densities(snps_hm, snps_ht, hyp_ratios, region, file)
		Plot.comparison(real_ratios, hyp_ratios, genome_len, file, ylim, original_pos, outcome_pos)
		Plot.qqplot(snps_hm, file, "QQplot for hm density", "Theoretical normal distribution", "Hypothetical SNP density", "hm_snps")
		Plot.qqplot(hyp_ratios, file, "QQplot for the ratios", "Theoretical normal distribution", "Hypothetical ratios", "ratios")

    # write the candidate mutation and related ratio information to files
    FileRW.write_txt("#{file}/hyp_ratios", hyp_ratios)
    FileRW.write_txt("#{file}/ratios", real_ratios)
    File.open("#{file}/mutation.txt", "w+") do |f|
      f.puts "The length of the group of contigs that form the peak of the distribution is #{region.to_i} bp"
      f.puts "The mutation is likely to be found on the following contigs #{candidate_mutations}"
      f.puts "Likeliest candidate responsible for mutation from the mutations is at #{candidate_peak}"
    end
  end

end
