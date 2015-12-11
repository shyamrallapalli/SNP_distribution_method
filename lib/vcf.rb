# encoding: utf-8
require 'bio'
require 'bio-samtools'
require_relative 'file_rw'

class Vcf

  def self.get_allele_freq(vcf_obj)
    allele_freq = 0
    if vcf_obj.info.key?('DP4')
      freq = vcf_obj.info['DP4'].split(',')
      depth = freq.inject { | sum, n | sum + n.to_f }
      alt = freq[2].to_f + freq[3].to_f
      allele_freq = alt / depth
    elsif vcf_obj.info.key?('RD')
      alt = vcf_obj.info['AD'].to_f
      depth = vcf_obj.info['RD'].to_f + alt
      allele_freq = alt / depth
    elsif vcf_obj.info.key?('AD')
      info =  vcf_obj.info['AD']
      if info =~ /','/
        freq = vcf_obj.info['DP4'].split(',')
        allele_freq = freq[1].to_f / ( freq[0].to_f + freq[1].to_f )
      elsif vcf_obj.info.key?('GT')
        gt = vcf_obj.info['GT']
        if gt == '1/1'
          allele_freq = 1.0
        elsif gt == '0/1'
          allele_freq = 0.5
        end
      end
    elsif vcf_obj.info.key?('AF')
      allele_freq = vcf_obj.info['AF'].to_f
    else
      warn "not a known vcf format and\
 check that it is one sample vcf\n"
      exit
    end
    warn "#{allele_freq}\n"
    allele_freq
  end

  ##Input: vcf file
  ##Ouput: lists of hm and ht SNPS and hash of all fragments with variants
  def self.get_vars(vcf_file, ht_cutoff=0.5, hm_cutoff=1.0)
    # hash of :het and :hom with frag ids and respective variant positions
    var_pos = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    File.open(vcf_file, 'r').each do |line|
      next if line =~ /^#/
      v = Bio::DB::Vcf.new(line)
      allele_freq = get_allele_freq(v)
      if allele_freq == ht_cutoff
        if var_pos[:het].has_key?(v.chrom)
          var_pos[:het][v.chrom] << v.pos
        else
          var_pos[:het][v.chrom] = []
          var_pos[:het][v.chrom] << v.pos
        end
      elsif allele_freq == hm_cutoff
        if  var_pos[:hom].has_key?(v.chrom)
          var_pos[:hom][v.chrom] << v.pos
        else
          var_pos[:hom][v.chrom] = []
          var_pos[:hom][v.chrom] << v.pos
        end
      end
    end

    var_pos
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
  def self.varpos_aggregate(frag_info, frag_len, frag_order, ratio_adjust, cumulate='yes')
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
      if cumulate == 'yes'
        asmbly_len += frag_len[frag].to_i
      end
    }
    details
  end

  # function to get seperate array of cumulative variant positions
  # input: a hash of frag ids with all details and variant positions
  # hash input is resutled from varpos_aggregate method
  # output: array of heterozygous and homozygous var positions
  def self.varpositions(fragdetails)
    hm_list = []
    ht_list = []
    fragdetails.keys.each { | frag |
      hm_list << fragdetails[frag][:hm_pos]
      ht_list << fragdetails[frag][:ht_pos]
    }
    [hm_list.flatten!, ht_list.flatten!]
  end

end
