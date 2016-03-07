#encoding: utf-8
require 'csv'

class RatioFilter

  def self.prune_hash(inhash)
    first = inhash.length
    inhash.each_key do | frag |
      numhm = inhash[frag][:hm]
      numht = inhash[frag][:ht]
      if numht >= numhm
        inhash.delete(frag)
      end
    end
    last = inhash.length
    warn "Discarded #{first-last} out of #{first} fragments\n"
    inhash
  end

  def self.selected_ratios(inhash, select_hmes=true)
    ratios_hash = {}
    # select only frag with higher hmes (homozygous enrichment score)
    # this is to discard fragments which may not contain information
    #  about causative mutations
    if select_hmes
      inhash = prune_hash(inhash)
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
