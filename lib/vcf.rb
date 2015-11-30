# encoding: utf-8
require 'bio'
require 'bio-samtools'
require_relative 'stuff'

class Vcf

  ##Input: vcf file
	##Ouput: lists of hm and ht SNPS
	def self.snps_in_vcf(vcf_file, ht_cutoff=0.5, hm_cutoff=1.0)
		vcfs_chrom, vcfs_pos, vcfs_info, hm, ht = [], [], [], [], []
		frag_pos = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
		File.open(vcf_file, "r").each do |line| # get array of vcf lines, you can call a method on one line
			next if line =~ /^#/
			v = Bio::DB::Vcf.new(line)
			vcfs_chrom << v.chrom
			vcfs_pos << v.pos
			vcfs_info << v.info # so this will be an array of hashes of strings
			allele_freq = v.info["AF"].to_f
      if allele_freq == ht_cutoff
		    if frag_pos[:het].has_key?(v.chrom)
		      frag_pos[:het][v.chrom] << v.pos
		    else
		      frag_pos[:het][v.chrom] = []
		      frag_pos[:het][v.chrom] << v.pos
		    end
		    ht << v.chrom
      elsif allele_freq == hm_cutoff
		    if  frag_pos[:hom].has_key?(v.chrom)
		      frag_pos[:hom][v.chrom] << v.pos
		    else
		      frag_pos[:hom][v.chrom] = []
		      frag_pos[:hom][v.chrom] << v.pos
		    end
		    hm << v.chrom
      end
		end
    # putting the number of snps for each frag into hash
    # frag_id is the key, the number of snps for that frag is the value
    num_snps_frag_hash = Stuff.create_hash_number(vcfs_chrom)
    snp_data = vcfs_chrom, vcfs_pos, num_snps_frag_hash, vcfs_info
		[snp_data, hm, ht, frag_pos]
	end

  def self.dic_id_pos(h_ids, snp_list)
		dic_pos = {}
	  (0..snp_list.length - 1).each do |x|
	  	if dic_pos.has_key?(h_ids[x])
	  		dic_pos[h_ids[x]] << snp_list[x]
	  	else
	  		dic_pos[h_ids[x]] = []
	  		dic_pos[h_ids[x]] << snp_list[x]
	  	end
	 	end
	  dic_pos
  end

  ##Input 1: Array of fragment ids.
	##Input 2: Hash of hm SNPs
	##Input 3: Hash of ht SNPs
	##Assign the number of SNPs to each fragment.
	##If a fragment does not have SNPs, the value assigned will be 0.
	##Output1: New hash with the number of SNPs assigned to the unordered fragments
	##Output2: Array with the number of SNPs
  	def self.define_snps(ids, dic_snps)
		dic_snps_num = {}
		snps = []
		ids.each { |frag|
	      if dic_snps.has_key?(frag)
	        dic_snps_num.store(frag, dic_snps[frag].to_f)
	      else
	        dic_snps_num.store(frag, 0)
	      end
	  	}
		dic_snps_num.each { |id, snp| snps << snp }
		return dic_snps_num, snps
	end

=begin
	def self.define_global_pos(ids, dic_pos, lengths)
		global_array, lens = [], []
	  or_pos, temporal_pos, dic_global_pos, idlen_s = {}, {}, {}, {}

	  # pp ids, dic_pos, lengths
		# ids.each { |frag|
		#     if dic_pos.has_key?(frag)
		#       or_pos.store(frag, dic_pos[frag])
		#     end
	  #  	}
	  #  	keys = or_pos.keys
	  	# temporal_pos.store(keys[0], or_pos[keys[0]])
	  	# or_pos.delete_if { |frag, pos|  temporal_pos.has_key?(frag)}
	  	# id_length.each do |frag, length|
	  	# 	idlen_s.store(frag, length)
	  	# 	lens << length
	  	# end
	  	x = 1
	  	lengths.length.times do |i|
	  		lengths[x] = lengths[x-1].to_i + lengths[x].to_i
	  		x += 1
	  	end

	  	x = 0
    	dic_pos.each do |frag, array|
	  		array.each do |pos|
	  			pos2 = (pos.to_i+lengths[x].to_i)
	  			global_array << pos2
	  		end
	  		x += 1
	  		dic_global_pos.store(frag, global_array)
  			global_array = []
  		end
	  	#dic_global_pos = temporal_pos.merge(dic_global_pos)
	  	all_global_positions = dic_global_pos.values
	  	all_global_positions.flatten!
		return dic_global_pos, all_global_positions
	end
=end

	def self.positions_by_fragment(dic, snp_list)
		dic.each do |id, number|
			dic.store(id, snp_list[0..(number.to_i-1)])
			(number.to_i).times do
		    	snp_list.delete_at(0)
			end
		end
		dic.delete_if { |id, number|  number.empty?}
		return dic
	end

  def self.open_vcf(vcf_file, chromosome)
    vcfs_chrom = []
    vcfspos = []
    vcfsinfo = []
    new_vcf = []
    File.open(vcf_file, 'r').each do |line|
      next if line =~ /^#/
      v = Bio::DB::Vcf.new(line)
      vcfs_chrom << v.chrom
      vcfspos << v.pos
      vcfsinfo << v.info
      # a = line.split("\t")
      new_vcf << line if v.chrom == "#{chromosome}"
    end
    [new_vcf, vcfs_chrom, vcfspos, vcfsinfo]
  end

  def self.type_per_pos(vcfs_info, vcfs_pos)
    snps = {}
    x = 0
    vcfs_info.each do |hash|
      hash.each do |type, number|
        if number == '1'
          snps.store(vcfs_pos[x], type)
          x += 1
        end
      end
    end
    hm = []
    ht = []
    snps.each do |pos, type|
      if type == 'HET'
        ht << pos
      elsif type == 'HOM'
        hm << pos
      end
    end
    [snps, hm, ht]
  end

  def self.filtering(vcfs_pos_c, snps_p, snps_c, child_chr_vcf)
    short_vcfs_pos_c = vcfs_pos_c
    short_vcfs_pos_c.flatten!
    snps_p.each do |pos, _type|
      if snps_c.key?(pos)
        snps_c.delete(pos)
        short_vcfs_pos_c.delete(pos)
      end
    end
    short_child_chr_vcf = []
    child_chr_vcf.each do |line|
      position = line.split("\t")[1].to_i
      short_child_chr_vcf << line if short_vcfs_pos_c.include?(position)
    end
    short_child_chr_vcf
  end

  # function to get cumulative variant positions from the order of fragments
  # input1: hash of frag ids with positions for homozygous and heterozygous variants
  # input2: hash of fragment lengths
  # input3: array of fragment order
  # input4: ratio adjustment factor
  # output: a hash of frag ids with all details and variant positions
  # are accumulated using length and order of fragments
  def self.varpos_aggregate(frag_info, frag_len, frag_order, ratio_adjust)
    details = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    asmbly_len = 0
    frag_order.each { | frag |
      details[frag][:hm] = ratio_adjust
      details[frag][:ht] = ratio_adjust
      details[frag][:hm_pos] = []
      details[frag][:ht_pos] = []
      if frag_info[:hom].key?(frag)
        hm_pos = frag_info[:hom][frag]
        details[frag][:hm] += hm_pos.length
        details[frag][:hm_pos] = hm_pos.map { |position| position + asmbly_len }
      end
      if frag_info[:het].key?(frag)
        ht_pos = frag_info[:het][frag]
        details[frag][:ht] += ht_pos.length
        details[frag][:ht_pos] = ht_pos.map { |position| position + asmbly_len }
      end
      if details[frag][:hm] == ratio_adjust and details[frag][:ht] == ratio_adjust
        details[frag][:ratio] = 0.0
      else
        details[frag][:ratio] = details[frag][:hm]/details[frag][:ht]
      end
      details[frag][:len] = frag_len[frag].to_i
      asmbly_len += frag_len[frag].to_i
    }
    details
  end

end
