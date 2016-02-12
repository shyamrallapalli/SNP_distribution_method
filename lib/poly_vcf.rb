# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Polyploid

  # getting vars from vcf file
  # each var is checked using pileup information from bam file
  # added to a hash to return
  def self.vars_in_file(vcf_file, bamfile, fastafile)
    # hash of frag ids with respective variant positions and their base hash info
    # only snps have base hash info and indels base hash is empty
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

    # open bam object to access pileup information
    bam = Bio::DB::Sam.new(:bam=>bamfile, :fasta=>fastafile)
    bam.open

    # read vcf file and process each variant
    File.open(vcf_file, 'r').each do |line|
      next if line =~ /^#/
      v = Bio::DB::Vcf.new(line)
      # some variant callers like Freebays are including some non variants as vars
      # so check ref and alt differ
      if v.variant? and v.ref != v.alt
        pileups = Pileup.get_pileup(bam,v.chrom,v.pos)
        next if pileups.empty?
        basehash = Pileup.read_base_hash(pileups[0])
        vars_hash[v.chrom][v.pos] = basehash
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
      freq = data[base].to_f/coverage
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


  def self.push_base_hash(base_hash, store_hash, frag, pos)
    if base_hash.keys.length > 1
      warn "#{frag}\t#{pos}\t#{base_hash}\n"
      return store_hash
    end
    mut_type = var_mode(base_hash[base_hash.keys[0]])
    store_hash = Vcf.push_to_hash(store_hash, frag, pos, mut_type)
    store_hash
  end


  def self.filter_vars(vars_hash_mut, vars_hash_bg, depth=6, noise=0.1)
    vars_hash = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    vars_hash_mut.each_key do | frag |
      positions = vars_hash_mut[frag].keys
      positions.each do | pos |
        data1 = vars_hash_mut[frag][pos]
        # background bulk has variant at same position
        if vars_hash_bg[frag].key?(pos)
          data2 = vars_hash_bg[frag][pos]
          if data1.instance_of? Hash
            mut_bases = get_base_freq(data1, depth, noise)
            reps = mut_bases.length
            if data2.instance_of? Hash
              bg_bases = get_base_freq(data2, depth, noise)
              if reps == 1
                base = mut_bases.keys[0]
                if bg_bases.key?(base)
                  mut_type = var_mode(mut_bases[base])
                  bg_type = var_mode(bg_bases[base])
                  if mut_type == bg_type
                    next
                  else
                    vars_hash = Vcf.push_to_hash(vars_hash, frag, pos, mut_type)
                  end
                else
                  mut_type = var_mode(mut_bases[base])
                  vars_hash = Vcf.push_to_hash(vars_hash, frag, pos, mut_type)
                end
              else
                warn "#{frag}\t#{pos}\t#{mut_bases}\t#{bg_bases}\n"
              end
            else
              vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos)
            end
          elsif data1.instance_of? String
            mut_ratio = Pileup.get_nonref_ratio(data1)
            mut_type = var_mode(mut_ratio)
            if data2.instance_of? String
              bg_ratio = Pileup.get_nonref_ratio(data2)
              bg_type = var_mode(bg_ratio)
              if mut_type == bg_type
                next
              else
                vars_hash = Vcf.push_to_hash(vars_hash, frag, pos, mut_type)
              end
            else
              vars_hash = Vcf.push_to_hash(vars_hash, frag, pos, mut_type)
            end
          end
        # only mut bulk has variant at this position
        else
          # var is snp
          if data1.instance_of? Hash
            mut_bases = get_base_freq(data1, depth, noise)
            vars_hash = push_base_hash(mut_bases, vars_hash, frag, pos)
          # var is indel
          elsif data1.instance_of? String
            mut_ratio = Pileup.get_nonref_ratio(data1)
            mut_type = var_mode(mut_ratio)
            vars_hash = Vcf.push_to_hash(vars_hash, frag, pos, mut_type)
          end
        end
      end
    end
    vars_hash_mut
  end

end
