# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Vcf

  DEFAULT = {
      ht_low: 0.25,
      ht_high: 0.75,
      cumulate: true,
  }

  def self.get_allele_freq(vcf_obj)
    allele_freq = 0
    # check if the vcf is from samtools (has DP4 and AF1 fields in INFO)
    if vcf_obj.info.key?('DP4')
      # freq = vcf_obj.info['DP4'].split(',')
      # depth = freq.inject { | sum, n | sum + n.to_f }
      # alt = freq[2].to_f + freq[3].to_f
      # allele_freq = alt / depth
      allele_freq = vcf_obj.non_ref_allele_freq
    # check if the vcf has has AF fields in INFO
    elsif vcf_obj.info.key?('AF')
      allele_freq = vcf_obj.info['AF'].to_f
    # check if the vcf is from VarScan (has RD, AD and FREQ fields in FORMAT)
    elsif vcf_obj.samples['1'].key?('RD')
      alt = vcf_obj.samples['1']['AD'].to_f
      depth = vcf_obj.samples['1']['RD'].to_f + alt
      allele_freq = alt / depth
    # check if the vcf is from GATK (has AD and GT fields in FORMAT)
    elsif vcf_obj.samples['1'].key?('AD')
      info =  vcf_obj.samples['1']['AD']
      if info =~ /','/
        freq = vcf_obj.samples['1']['AD'].split(',')
        allele_freq = freq[1].to_f / ( freq[0].to_f + freq[1].to_f )
      elsif vcf_obj.info.key?('AF')
        allele_freq = vcf_obj.info['AF'].to_f
      elsif vcf_obj.samples['1'].key?('GT')
        gt = vcf_obj.samples['1']['GT']
        if gt == '1/1'
          allele_freq = 1.0
        elsif gt == '0/1'
          allele_freq = 0.5
        end
      end
    else
      warn "not a known vcf format and \
check that it is one sample vcf\n"
      exit
    end
    allele_freq
  end


  ##Input: vcf file
  ##Ouput: lists of hm and ht SNPS and hash of all fragments with variants
  def self.get_vars(vcf_file, opts = {})
    opts = DEFAULT.merge(opts)
    ht_low = opts[:ht_low]
    ht_high = opts[:ht_high]

    # hash of :het and :hom with frag ids and respective variant positions
    var_pos = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    File.open(vcf_file, 'r').each do |line|
      next if line =~ /^#/
      v = Bio::DB::Vcf.new(line)
      if v.variant?
        allele_freq = get_allele_freq(v)
        if allele_freq.between?(ht_low, ht_high)
          var_pos[:het][v.chrom][v.pos] = allele_freq
        elsif allele_freq > ht_high
          var_pos[:hom][v.chrom][v.pos] = allele_freq
        end
      end
    end
    var_pos
  end

  def self.filtering(mutant_vcf, parent_vcf)
    var_pos_mut = get_vars(mutant_vcf)
    var_pos_bg = get_vars(parent_vcf)

    var_pos_mut.each_key do | type |
      var_pos_mut[type].each_key do | frag |
        positions = var_pos_mut[type][frag].keys
        pos_bg_bulk = var_pos_bg[type][frag].keys
        positions.each do |pos|
          if pos_bg_bulk.include?(pos)
            var_pos_mut[type][frag].delete(pos)
            positions.delete(pos)
          end
        end
        var_pos_mut[type].delete(frag) if positions.empty?
      end
    end
    var_pos_mut
  end

  # function to get cumulative variant positions from the order of fragments
  # input1: hash of frag ids with positions for homozygous and heterozygous variants
  # input2: hash of fragment lengths
  # input3: array of fragment order
  # input4: ratio adjustment factor
  # output: a hash of frag ids with all details and variant positions
  # are accumulated using length and order of fragments
  def self.varpos_aggregate(frag_info, frag_len, frag_order, ratio_adjust,  opts = {})
    opts = DEFAULT.merge(opts)
    cumulate = opts[:cumulate]

    details = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    cumulate_len = 0
    frag_order.each { | frag |
      details[frag][:hm] = ratio_adjust
      details[frag][:ht] = ratio_adjust
      details[frag][:hm_pos] = []
      details[frag][:ht_pos] = []
      details[frag][:cum_len] = cumulate_len
      if frag_info[:hom].key?(frag)
        hm_pos = frag_info[:hom][frag].keys
        details[frag][:hm] += hm_pos.length
        details[frag][:hm_pos] = hm_pos.map { |position| position + cumulate_len }
      end
      if frag_info[:het].key?(frag)
        ht_pos = frag_info[:het][frag].keys
        details[frag][:ht] += ht_pos.length
        details[frag][:ht_pos] = ht_pos.map { |position| position + cumulate_len }
      end
      if details[frag][:hm] == ratio_adjust and details[frag][:ht] == ratio_adjust
        details[frag][:ratio] = 0.0
      else
        details[frag][:ratio] = details[frag][:hm]/details[frag][:ht]
      end
      details[frag][:len] = frag_len[frag].to_i
      if cumulate
        cumulate_len += frag_len[frag].to_i
      end
    }
    details
  end

end
