#encoding: utf-8
require 'bio'

class FastaHandle

  # Input: Array of Bio::FastaFormat entries
  # Output 0: Array of identifiers
  # Output 1: Array of lengths (integers)
  def self.fasta_id_n_lengths(fasta)
    ids, lengths = [], []
    id_len = {}
    fasta.each do |i|
      ids << i.entry_id
      lengths << i.length
      id_len.store(i.entry_id, i.length)
    end
    return ids, lengths, id_len
  end

  # Input: FASTA file
  # Output: Array of Bio::FastaFormat entries
  def self.array(fasta_file)
    fasta = [] # we have the lengths of each fasta, but the frags are different to those of the vcf/hash(this only has the frags w snps)
    Bio::FastaFormat.open(fasta_file).each do |i| # get array of fasta format frags, ##  WE NEED TO REORDER THE FASTA FRAGS HERE, TO TEST DIFFERENT ARRANGEMENTS
      fasta << i
    end
    return fasta
  end

  # Input: FASTA file
  # Output: Integer of the length of the genome
  def self.genome_length(fasta_file)
    lengths = []
    FastaHandle.array(fasta_file).each do |frag|
      lengths << frag.length
    end
    return lengths.inject(:+)
  end

  # Input: FASTA file
  # Output: hash of sequence ids with lengths and seqeunces and total assmebly length
  def self.fasta_parse(fasta_file)
    sequences = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    assembly_len = 0
    Bio::FastaFormat.open(fasta_file).each do |inseq|
      sequences[:seq][inseq.entry_id] = inseq.entry
      sequences[:len][inseq.entry_id] = inseq.length
      assembly_len += inseq.length
    end
    return sequences, assembly_len
  end

  # Input1: permutation array of frag ids after SDM
  # Input2: Fasta file seq hash with ids as keys
  # Output: permutation of fragments after SDM with the data (lengths, etc) obtained from the original fasta file
  def self.create_perm_fasta(perm, frags)
    fasta_perm = []
    perm.each do |frag|
      fasta_perm << Bio::FastaFormat.new(frags[frag].to_s) if frags.key?(frag)
    end
    fasta_perm
  end

end
