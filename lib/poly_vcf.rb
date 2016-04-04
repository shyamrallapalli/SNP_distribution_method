# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Polyploid

  attr_accessor :polyploidy, :ht_low, :ht_high, :min_depth, :noise

  DEFAULT = {
      ignore_reference_n: true,
      min_depth: 6,
      min_non_ref_count: 3,
      noise: 0.1,
      polyploidy: false,
      ht_low: 0.1,
      ht_high: 0.9,
  }

  # getting vars from pileup file
  # each var is checked from pileup information
  # added to a hash to return
  def self.vars_in_pileup(pileupfile, opts = {})
    opts = DEFAULT.merge(opts)
    ignore_reference_n = opts[:ignore_reference_n]
    min_depth  = opts[:min_depth]
    min_non_ref_count = opts[:min_non_ref_count]

    # hash of frag ids with respective variant positions and their base hash info
    # only snps have base hash info and indels base hash is read bases
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

    # read mpileup file and process each variant
    File.open(pileupfile, 'r').each do |line|
      pileup = Bio::DB::Pileup.new(line)
      if pileup.is_snp?(:ignore_reference_n => ignore_reference_n, :min_depth => min_depth,
                        :min_non_ref_count => min_non_ref_count) and
          pileup.consensus != pileup.ref_base
        read_bases = Pileup.get_read_bases(pileup)
        basehash = Pileup.read_bases_to_hash(read_bases)
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
  def self.get_var_base_prop(hash)
    snp_hash = {}
    coverage = hash.values.inject { | sum, n | sum + n.to_f }
    return snp_hash if coverage < @min_depth
    # calculate proportion of each base in coverage
    hash.each_key do | base |
      next if base == :ref
      freq = hash[base].to_f/coverage
      next if freq <= @noise
      snp_hash[base] = freq
    end
    snp_hash
  end

  # calculate var zygosity for non-polyploid variants
  # increased range is used for heterozygosity for RNA-seq data
  def self.var_mode(ratio)
    mode = ''
    if ratio.between?(@ht_low, @ht_high)
      mode = :het
    elsif ratio > @ht_high
      mode = :hom
    end
    mode
  end

  # get total proportion of bases in hash
  def self.polybase_proportion(vars, hash)
    # if polyploidy set then take combination of proportions
    # if not then take maximum value
    if @polyploidy
      prop = 0.0
      hash.each_key { | key |
        if vars.include?(key)
          prop += hash[key]
        end
      }
      # warn "polybase\t#{hash}\t#{vars}\t#{prop}\n"
    else
      prop = hash.values.max
      # warn "non_polybase\t#{hash}\t#{vars}\t#{prop}\n"
    end
    prop
  end

  # method to compare base hash between and background and mutant
  # returns the var_type if single base left after background subtraction
  # otherwise returns empty for multiple vars
  def self.multi_var_hash(base_hash, background='')
    var_type = ''
    mut_vars = base_hash.keys.sort
    # ignore complex variant locations
    # return var_type if mut_vars.length > 2
    if background != ''
      bg_vars = background.keys.sort
      return var_type if mut_vars == bg_vars
      mut_vars.delete_if { |base| bg_vars.include?(base) }
      if mut_vars.length == 1
        var_type = var_mode(base_hash[mut_vars[0]])
      else
        ratio = polybase_proportion(mut_vars, base_hash)
        var_type = var_mode(ratio)
      end
    else
      ratio = polybase_proportion(mut_vars, base_hash)
      var_type = var_mode(ratio)
    end
    var_type
  end

  def self.push_base_hash(base_hash, store_hash, frag, pos, background='')
    # we are only dealing with single element hashes
    # so discard hashes with more than one element and empty hashes
    # empty hash results from position below selected coverage or bases freq below noise
    return store_hash if base_hash.empty?
    if base_hash.length > 1
      vartype = multi_var_hash(base_hash, background)
      # warn "#{frag}\t#{pos}\t#{base_hash}\t#{background}\n"
      if vartype == ''
        return store_hash
      else
        store_hash = Vcf.push_to_hash(store_hash, frag, pos, vartype)
      end
    else
      base = base_hash.keys[0]
      mut_type = var_mode(base_hash[base])
      if background != '' and background.key?(base)
        bg_type = var_mode(background[base])
        # if both have the same base type then return original hash
        return store_hash if mut_type == bg_type
        store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
      else
        store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
      end
    end
    store_hash
  end

  def self.filter_vars(mut_pileup, vars_hash_bg, opts = {})
  opts = DEFAULT.merge(opts)
  ignore_reference_n = opts[:ignore_reference_n]
  min_non_ref_count = opts[:min_non_ref_count]
  @min_depth  = opts[:min_depth]
  @noise = opts[:noise]
  @polyploidy = opts[:polyploidy]
  @ht_low = opts[:ht_low]
  @ht_high  = opts[:ht_high]

  vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    # read mpileup file and process each variant
    File.open(mut_pileup, 'r').each do |line|
      pileup = Bio::DB::Pileup.new(line)
      if pileup.is_snp?(:ignore_reference_n => ignore_reference_n, :min_depth => @min_depth, :min_non_ref_count => min_non_ref_count) and
          pileup.consensus != pileup.ref_base
        read_bases = Pileup.get_read_bases(pileup)
        data1 = Pileup.read_bases_to_hash(read_bases)
        frag = pileup.ref_name
        pos = pileup.pos
        mut_bases = get_var_base_prop(data1)
        if vars_hash_bg[frag].key?(pos)
          bg_bases = get_var_base_prop(vars_hash_bg[frag][pos])
          vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos, bg_bases)
        else
          vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos)
        end
      end
    end
    vars_hash
  end

end
