#!/usr/bin/env ruby

#
# find regions that are specific to a cell type
#
# 
#

require 'rubygems'
require 'trollop'

require_relative 'splice'
require_relative 'exon'

opts = Trollop::options do
  version "v0.0.1a"
  opt :gff, "Gff", :required => true, :type => String
  opt :data, "RSEM data as csv", :required => true, :type => String
  opt  :cell, "Cell type", :required => true, :type => String
  opt :output, "Output", :type => String
  opt :verbose, "Be verbose"
end

Trollop::die :gff, "must exist" if !File.exist?(opts[:gff]) if opts[:gff]
Trollop::die :data, "must exist" if !File.exist?(opts[:data]) if opts[:data]
#Trollop::die :output, "mustn't exist" if File.exist?(opts[:output]) if opts[:output]

# # # # # # # # # # ##  # # # # # # # # # # # # # # # #
# Load the mRNA and exon  lines from the gff file
# into a hash[][]  called  genes
#
genes = Hash.new 
print "Loading gff file..." if opts.verbose

File.open("#{opts.gff}", "r").each_line do |line|
  line.chomp!
  cols = line.split(/\t/)
  if cols[2] =~ /mRNA/
    if cols[8] =~ /ID=(\S+);Parent=(\S+);Name/
      splicename = $1
      genename = $2
      if !genes.has_key?(genename)
        genes[genename] = Hash.new
      end
      s = Splice.new(splicename, cols[0], cols[3], cols[4], cols[6])
      genes[genename][splicename] = s 
    end

  elsif cols[2] =~ /exon/
    if cols[8] =~ /Parent=(\S+)/
      splicename = $1
      genename = splicename.split(/\./).first
      if genes.has_key?(genename)
        genes[genename][splicename].exons << Exon.new(cols[3], cols[4], "exon")
      end
    end
  end
end

puts "Done" if opts.verbose

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Load the RSEM data with isoform expression
#

data = Hash.new
print "Loading gene list..." if opts.verbose
File.open("#{opts.data}", "r").each_line do |line|
  line.chomp!
  cols = line.split(/\t/)
  splicename = cols[0]

  genename = splicename.split(/\./).first
  if genes.has_key?(genename)
    if genes[genename].has_key?(splicename)
      genes[genename][splicename].cell = cols[12]    # cell specific?
      genes[genename][splicename].gdc = cols[9].to_f # expression
      genes[genename][splicename].s35 = cols[10].to_f
      if cols[9].to_f > 0 or cols[10].to_f > 0
        data[genename]=1
      end
    end
  end
end
puts "Done" if opts.verbose

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
output=""
uni=[]
data.each_key do |genename|
  # puts "genename = #{genename}" if opts.verbose
  # make list of splices that are DE in cell type opts.cell
  splice_list_cell = []
  splice_list_not  = []
  genes[genename].each_pair do |splicename, splice|
    if splice.cell == opts.cell 
      splice_list_cell << splice
    else
      splice_list_not << splice
    end
  end
  
  # create intersection of all cell specific isoforms
  if splice_list_cell.length > 0
    common_splice = splice_list_cell.first.dup
    splice_list_cell.each do |splice|
      common_splice = common_splice.intersect(splice)
    end

    if splice_list_not.length > 0  
      splice_list_not.each do |splice_B|
        common_splice.minus!(splice_B)
      end
      uni << common_splice
      output << "#{common_splice}\n"
    end
  end
end # data.each_key 

File.open("#{opts.output}", "w") {|file| file.write(output)}

# list = Array.new(1001) {|i| i=0}

# uni.each do |splice|
#   if splice.exons.length >0
#     start = 0
#     stop = splice.stop - splice.start
#     length = stop+1
#     splice.exons.each do |exon|
#       exon_start = exon.start - splice.start
#       exon_stop = exon.stop - splice.start
#       from = (1000*exon_start.to_f / length).round
#       to   = (1000*exon_stop.to_f  / length).round
#       (from..to).each do |i|
#         list[i]+=1
#       end
#     end
#   end
# end

# # 0-1000 graph showing where specific exons are
# list.each_index{|i| puts "#{i}\t#{list[i]}"}

