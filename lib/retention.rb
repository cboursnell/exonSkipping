#!/usr/bin/env ruby

#
# find intron retention and exon skipping
#

# requires stuff

require 'rubygems'
require 'trollop'

require_relative 'splice'
require_relative 'exon'

if __FILE__ == $0
  opts = Trollop::options do
    version "v0.0.1a"
    opt :gff, "Gff file", :required => true, :type => String
    opt :verbose, "Be verbose"
  end

  Trollop::die :gff, "must exist" if !File.exist?(opts[:gff]) if opts[:gff]

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

  retained = Hash.new
  genes.each_pair do |genename, splices|
    # puts "genename: #{genename}"
    list =[]
    splices.each_pair do |splicename, splice|
      list << splice
    end
    list.each_index do |i|
      list.each_index do |j|
        # puts "#{i} #{j}"
        exon_list = list[i].retention(list[j])
        unless exon_list==[]
          exon_list.each do |retained_exon|
            found=false
            retained.each_pair do |g, g_list|
              g_list.each do |e|
                if e.start == retained_exon.start and e.stop == retained_exon.stop
                  found=true
                end
              end
            end
            if !found
              if !retained.has_key?(list[i].name)
                retained[list[i].name]=[]
              end
              retained[list[i].name] << retained_exon
            end
          end
        end
      end
    end
  end

  retained.each_pair do |splicename, exon_list|
    puts "splicename: #{splicename}"
    exon_list.each do |exon|
      puts exon
    end
  end

end

