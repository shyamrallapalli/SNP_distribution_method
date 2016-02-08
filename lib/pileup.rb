# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class Pileup

  def self.pick_frag_vars(bamfile,infasta,keyfrags,input_frags)
    bam = Bio::DB::Sam.new(:bam=>bamfile, :fasta=>infasta)
    bam.open
    sortfrags = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    keyfrags.each do | selfrag |
      positions = input_frags[selfrag][:hm_pos]
      positions.each do | mutpos |
        bam.mpileup(:r => "#{selfrag}:#{mutpos}-#{mutpos}", :Q => 15, :q => 20) do |pileup|
          ratio = 0
          if pileup.is_snp?(:ignore_reference_n => true, :min_depth => 6, :min_non_ref_count => 3)
            # if defined? pileup.non_ref_count
            ratio = 1.0 - (pileup.ref_count/pileup.coverage)
          end
          sortfrags[ratio][selfrag][mutpos] = pileup
        end
      end
    end
    sortfrags
  end

end