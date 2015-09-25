# encoding: utf-8

# Input1: permutation array of frag ids after SDM
# Input2: Fasta file seq hash with ids as keys
# Output: permutation of fragments after SDM with the data (lengths, etc) obtained from the original fasta file
class Output
  def self.create_perm_fasta(perm, frags)
    fasta_perm = []
    perm.each do |frag|
      fasta_perm << Bio::FastaFormat.new(frags[frag].to_s) if frags.key?(frag)
    end
    fasta_perm
  end
end
