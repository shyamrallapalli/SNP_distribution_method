#encoding: utf-8

class SNPdist
	require 'rinruby'
	### Hypothetical SNP positions ###
	# Input 0: The contig size
	# Input 1: Array of ratios for homozygous to heterozygous SNPs at each contig.
	# Output: A list of positions for the given ratios array size.
	def self.general_positions(average_contig_size, ratios)
		positions = []
    pos = 0
		ratios.length.times do |i|
      pos += average_contig_size
			positions << pos
		end
		positions
	end

  # Input 0: The homozygous, heterozygous or ratio densities in the ordered genome.
	# Input 1: Previous output. A list of "hypothetical"  positions
	# Output: A list of "hypothetical" SNP positions which represents the distribution of hm, ht SNPs or  homozygous/heterozygous SNP density ratio
	def self.densities_pos(raw_densities, positions)
		y = 0
		densities_in_pos = []
		raw_densities.each do |density|
			(density*10).to_i.times do |i|
				densities_in_pos << positions[y]
			end
			y += 1
		end
		return densities_in_pos
	end

	def self.hyp_snps(ratios, genome_length)
		breaks = []
		(1..ratios.length).to_a.each do |i|
			breaks << (genome_length/ratios.length.to_f)*i
		end
		hyp, x = [], 0
		ratios.each do |ratio|
			(ratio*10).to_i.times do
				hyp << rand(genome_length/ratios.length.to_f) + breaks[x] # random value from within the range that the freq has been taken
			end
			x+=1
		end
		return hyp # These don't need to be unique or integers like the real SNPs, since they are just representing a distribution
	end

end
