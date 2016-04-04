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

end
