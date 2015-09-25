# encoding: utf-8

class Vcf
  require 'bio'
  require 'bio-samtools'

  def self.open_vcf(vcf_file, chromosome)
    vcfs_chrom = []
    vcfspos = []
    vcfsinfo = []
    new_vcf = []
    File.open(vcf_file, 'r').each do |line|
      next if line =~ /^#/
      v = Bio::DB::Vcf.new(line)
      vcfs_chrom << v.chrom
      vcfspos << v.pos
      vcfsinfo << v.info
      # a = line.split("\t")
      new_vcf << line if v.chrom == "#{chromosome}"
    end
    [new_vcf, vcfs_chrom, vcfspos, vcfsinfo]
  end

  def self.type_per_pos(vcfs_info, vcfs_pos)
    snps = {}
    x = 0
    vcfs_info.each do |hash|
      hash.each do |type, number|
        if number == '1'
          snps.store(vcfs_pos[x], type)
          x += 1
        end
      end
    end
    hm = []
    ht = []
    snps.each do |pos, type|
      if type == 'HET'
        ht << pos
      elsif type == 'HOM'
        hm << pos
      end
    end
    [snps, hm, ht]
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
end
