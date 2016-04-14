#encoding: utf-8
require 'csv'

class RatioFilter

  DEFAULT = {
      only_frag_with_vars: true,
      ratio_adjust: 0.5,
      ratio_type: 'ratio',
  }

  def self.filter_hash(inhash, adjust, ratio_type, only_frag_with_vars)
    discard = 0
    inhash.each_key do | frag |
      numhm = inhash[frag][:hm]
      numht = inhash[frag][:ht]
      numbfr = inhash[frag][:bfr]
      if only_frag_with_vars
        # selecting fragments which have a variant
        if ratio_type == :ratio and numht + numhm <= 2 * adjust
          inhash[frag][:discard] = true
          discard += 1
          # inhash.delete(frag)
        # in polyploidy scenario selecting fragments with at least one bfr position
        elsif ratio_type == :bfr_rat and numbfr <= 0
          inhash[frag][:discard] = true
          discard += 1
          # inhash.delete(frag)
        else
          inhash[frag][:discard] = false
        end
      else
        inhash[frag][:discard] = false
      end
    end
    if only_frag_with_vars
      warn "Discarded #{discard} out of #{inhash.length} fragments, which lack any variant\n"
    else
      warn "No filtering was applied to fragments that lack any variant or bfr\n"
    end
    inhash
  end

  def self.selected_ratios(inhash, opts = {})
    opts = DEFAULT.merge(opts)
    only_frag_with_vars = opts[:only_frag_with_vars]
    adjust = opts[:ratio_adjust]
    ratio_type = opts[:ratio_type].to_sym
    ratios_hash = {}
    # select only fragments with variants
    # this is to discard fragments which may not contain
    # any information about causative mutations
    inhash = filter_hash(inhash, adjust, ratio_type, only_frag_with_vars)
    inhash.keys.each do | frag |
      unless inhash[frag][:discard]
        ratio = inhash[frag][ratio_type]
        if ratios_hash.key?(ratio)
          ratios_hash[ratio] << frag
        else
          ratios_hash[ratio] = []
          ratios_hash[ratio] << frag
        end
      end
    end
    ratios_hash
  end

end
