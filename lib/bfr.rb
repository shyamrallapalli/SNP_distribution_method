# encoding: utf-8

class BFR

  attr_accessor :polyploidy, :ht_low, :ht_high, :min_depth
  attr_accessor :noise, :bfr_adj, :parent_hemi_hash

  DEFAULT = {
      ignore_reference_n: true,
      min_depth: 6,
      min_non_ref_count: 3,
      noise: 0.1,
      polyploidy: false,
      ht_low: 0.1,
      ht_high: 0.9,
      bfr_adj: 0.01,
      parent_hemi_hash: '',
  }

  # get bulk frequency ratio (bfr) for marked hemi snps only
  # ignore postions with complex variants
  def self.get_bfr(mut_hash, bg_hash='')
    if bg_hash != ''
      # checking if only two vars in base hash and that includes ref
      # checking if only one var in hemi snp
      # suggests enrichment for one of two alleles
      if mut_hash.length == 2 and mut_hash.key?(:ref) or mut_hash.length == 1
        bfr = calculate_bfr(mut_hash, bg_hash)
      elsif bg_hash.length == 2  and bg_hash.key?(:ref) or bg_hash.length == 1
        bfr = calculate_bfr(bg_hash, mut_hash)
      else # complex
        bfr = ''
      end
    elsif mut_hash.length == 2 and mut_hash.key?(:ref) or mut_hash.length == 1
      bfr = calc_fraction(mut_hash)[0]/ @bfr_adj
    else
      bfr = ''
    end
    bfr
  end

  # calculate bfr using both mutant and background bulk information
  def self.calculate_bfr(two_key_hash, other_hash)
    # fix :ref value if absent due to below noise depth
    unless two_key_hash.key?(:ref)
      two_key_hash[:ref] = 0
    end
    unless other_hash.key?(:ref)
      other_hash[:ref] = 0
    end
    frac_1, base = calc_fraction(two_key_hash)
    if other_hash.key?(base)
      sum = other_hash[base] + other_hash[:ref] + @bfr_adj
      frac_2 = (other_hash[base] + @bfr_adj)/sum
    else
      sum = other_hash[:ref] + @bfr_adj
      frac_2 = @bfr_adj/sum
    end
    if frac_1 > frac_2
      bfr = frac_1/frac_2
    else
      bfr = frac_1/frac_2
    end
    bfr
  end

  def self.calc_fraction(hash)
    unless hash.key?(:ref)
      hash[:ref] = 0
    end
    array = hash.keys
    sum = hash[array[0]] + hash[array[1]] + @bfr_adj
    if array[0] == :ref
      frac = (hash[array[1]] + @bfr_adj)/sum
      base = array[1]
    else
      frac = (hash[array[0]] + @bfr_adj)/sum
      base = array[0]
    end
    [frac, base]
  end

end
