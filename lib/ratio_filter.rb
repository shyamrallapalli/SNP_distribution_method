#encoding: utf-8
require 'csv'

class RatioFilter

  def self.prune_hash(inhash, adjust)
    first = inhash.length
    inhash.each_key do | frag |
      numhm = inhash[frag][:hm]
      numht = inhash[frag][:ht]
      # selecting fragments which have a variant
      if numht + numhm <= 2 * adjust
        inhash.delete(frag)
      end
    end
    last = inhash.length
    warn "Discarded #{first-last} out of #{first} fragments, which lack any variant\n"
    inhash
  end

  def self.selected_ratios(inhash, adjust, select_hmes=true)
    ratios_hash = {}
    # select only frag with higher hmes (homozygous enrichment score)
    # or select only fragments with a variant
    # this is to discard fragments which may not contain
    # any information about causative mutations
    if select_hmes
      inhash = prune_hash(inhash, adjust)
    end
    inhash.keys.each do | frag |
      ratio = inhash[frag][:ratio]
      if ratios_hash.key?(ratio)
        ratios_hash[ratio] << frag
      else
        ratios_hash[ratio] = []
        ratios_hash[ratio] << frag
      end
    end
    ratios_hash
  end

end
