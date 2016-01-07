# encoding: utf-8
require 'bio'
require 'bio-samtools'
require 'bio-gngm'

class FileRW


  # Input: file location
  # Output: array with each line of the file as an entry (strings)
  def self.to_array(file)
    IO.foreach(file).collect {|l| l.chomp }
  end

  # Input 0: Filename by which to save an array with filetype extension, one value per line
  # Input 1: Array to save
  def self.write(filename, array)
    File.open("#{filename}", 'w+') do |f|
      array.each { |i| f.puts(i) }
    end
  end

  # Input 0: Filename by which to save an array to a .txt file, one value per line
  # Input 1: Array to save
  def self.write_txt(filename, array)
    write("#{filename}.txt", array)
  end

  # Input: FASTA file
  # Output: hash of sequence ids with lengths and sequences
  def self.fasta_parse(fasta_file)
    sequences = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
    Bio::FastaFormat.open(fasta_file).each do |inseq|
      sequences[:seq][inseq.entry_id] = inseq.entry
      sequences[:len][inseq.entry_id] = inseq.length
    end
    sequences
  end

  # Input1: permutation array of frag ids after SDM
  # Input2: sequences hash of frag ids and a filename to write
  # Output will be written to file and if no filename is given
  # output will be written to "ordered_frags.fasta"
  def self.write_order(array, fasta_file, filename='ordered_frags.fasta')
    fasta_db = Bio::DB::Fasta::FastaFile.new({:fasta=>fasta_file})
    samtools = fasta_db.instance_variable_get(:@samtools)
    File.open(filename, 'w+') do |f|
      array.each do |frag|
        # element = Bio::FastaFormat.new(seqhash[frag].to_s)
        # seqout = Bio::Sequence::NA.new(element.seq).upcase
        # f.puts seqout.to_fasta(element.definition, 80)
        command = "#{samtools} faidx #{fasta_file} \"#{frag}\""
        f.puts `#{command}`
      end
    end
  end

end
