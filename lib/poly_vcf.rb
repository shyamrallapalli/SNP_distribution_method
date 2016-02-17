# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Polyploid

  # getting vars from pileup file
  # each var is checked from pileup information
  # added to a hash to return
  def self.vars_in_pileup(pileupfile)
    # hash of frag ids with respective variant positions and their base hash info
    # only snps have base hash info and indels base hash is read bases
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

    # read mpileup file and process each variant
    File.open(pileupfile, 'r').each do |line|
      pileup = Bio::DB::Pileup.new(line)
      if pileup.is_snp?(:ignore_reference_n => true, :min_depth => 6, :min_non_ref_count => 3) and pileup.consensus != pileup.ref_base
        read_bases = Pileup.get_read_bases(pileup)
        basehash = Pileup.read_base_hash(read_bases)
        vars_hash[pileup.ref_name][pileup.pos] = basehash
        # puts "#{pileup.ref_name}\t#{pileup.pos}\t#{pileup.consensus}\t#{basehash}\n"
      end
    end
    vars_hash
  end

  # form hash of base information, [ATGC] counts for snp
  # a hash of base proportion is calculated
  # base proportion hash below a selected depth is empty
  # base proportion below or equal to a noise factor are discarded
  def self.get_base_freq(hash, depth, noise)
    snp_hash = {}
    coverage = hash.values.inject { | sum, n | sum + n.to_f }
    return snp_hash if coverage < depth
    # calculate proportion of each base in coverage
    hash.each_key do | base |
      next if base == :ref
      freq = hash[base].to_f/coverage
      next if freq <= noise
      snp_hash[base] = freq
    end
    snp_hash
  end

  # calculate var zygosity for non-polyploid variants
  # increased range is used for heterozygosity for RNA-seq data
  def self.var_mode(ratio, ht_low = 0.10, ht_high = 0.90)
    mode = ''
    if ratio.between?(ht_low, ht_high)
      mode = :het
    elsif ratio > ht_high
      mode = :hom
    end
    mode
  end


  def self.push_base_hash(base_hash, store_hash, frag, pos, background='')
    # we are only dealing with single element hashes
    if base_hash.keys.length > 1
      warn "#{frag}\t#{pos}\t#{base_hash}\t#{background}\n"
      return store_hash
    end
    base = base_hash.keys[0]
    if background == ''
      mut_type = var_mode(base_hash[base])
      store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
    else
      if background.key?(base)
        mut_type = var_mode(base_hash[base])
        bg_type = var_mode(background[base])
        # if both have the same base type then return original hash
        return store_hash if mut_type == bg_type
        store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
      else
        mut_type = var_mode(base_hash[base])
        store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
      end
    end
    store_hash
  end

  def self.filter_vars(mut_pileup, vars_hash_bg, depth=6, noise=0.1)
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    # read mpileup file and process each variant
    File.open(mut_pileup, 'r').each do |line|
      pileup = Bio::DB::Pileup.new(line)
      if pileup.is_snp?(:ignore_reference_n => true, :min_depth => 6, :min_non_ref_count => 3) and pileup.consensus != pileup.ref_base
        read_bases = Pileup.get_read_bases(pileup)
        data1 = Pileup.read_base_hash(read_bases)
        frag = pileup.ref_name
        pos = pileup.pos
        if data1.instance_of? Hash
          mut_bases = get_base_freq(data1, depth, noise)
          if vars_hash_bg[frag].key?(pos) and vars_hash_bg[frag][pos].instance_of? Hash
            bg_bases = get_base_freq(vars_hash_bg[frag][pos], depth, noise)
            vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos, bg_bases)
          else
            vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos)
          end
        elsif data1.instance_of? String
          mut_ratio = Pileup.get_nonref_ratio(data1)
          mut_type = var_mode(mut_ratio)
          if vars_hash_bg[frag].key?(pos) and vars_hash_bg[frag][pos].instance_of? String
            bg_ratio = Pileup.get_nonref_ratio(vars_hash_bg[frag][pos])
            bg_type = var_mode(bg_ratio)
            next if mut_type == bg_type
            vars_hash = Vcf.push_to_hash(vars_hash, frag, pos, mut_type)
          else
            vars_hash = Vcf.push_to_hash(vars_hash, frag, pos, mut_type)
          end
        else
          warn "I don't know the type\t#{data1.class}"
        end
      end
    end
    vars_hash
  end

end
