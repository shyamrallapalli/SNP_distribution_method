#encoding: utf-8
require 'bio'

class FastaHandle

  # Input: FASTA file
  # Output: hash of sequence ids with lengths and seqeunces and total assmebly length
  def self.file_parse(fasta_file)
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
  def self.write_order(perm, frags, filename='ordered_frags.fasta')
    fasta_perm = []
    File.open(filename, 'w+') do |f|
      perm.each do |frag|
        element = Bio::FastaFormat.new(frags[frag].to_s)
        seqout = Bio::Sequence::NA.new(element.seq).upcase
        f.puts seqout.to_fasta(element.definition, 80)
      end
    end
  end

end
