#!/usr/bin/env ruby

#
# find intron retention and exon skipping
#


require 'rubygems'
require 'trollop'

require_relative 'splice'
require_relative 'exon'

if __FILE__ == $0
  opts = Trollop::options do
    version "v0.0.1a"
    opt :gff, "Gff file", :required => true, :type => String
    opt :data, "Data file", :required => true, :type => String
    opt :output, "Output prefix", :default => "output", :type => String
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
        genes[genename][splicename].cell = cols[12]     # cell specific?
        genes[genename][splicename].gdc = cols[9].to_f  # expression
        genes[genename][splicename].s35 = cols[10].to_f # expression
        if cols[9].to_f > 0 or cols[10].to_f > 0
          data[genename]=1
        end
      end
    end
  end
  puts "Done" if opts.verbose

  skipped = Hash.new
  retained = Hash.new
  five_prime = Hash.new
  three_prime = Hash.new
  genes.each_pair do |genename, splices|
    # puts "genename: #{genename}"
    list =[]
    splices.each_pair do |splicename, splice|
      # puts "  splicename: #{splicename} #{splice.cell}"
      list << splice
    end
    list.each_index do |i|
      list.each_index do |j|
        if i!=j
          # puts "#{i} #{j} #{list[i].cell}"
          exon_list = list[i].skipping(list[j]) # returns list of exons from list[i] that have no overlap in list[j]
          exon_list_ret = list[i].retention(list[j])
          exon_list_five_prime  = list[i].five_prime(list[j])
          exon_list_three_prime = list[i].three_prime(list[j])
          # exon skipping
          unless exon_list==[] 
            exon_list.each do |skipped_exon|
              found=false
              if skipped.has_key?(list[i].name)
                skipped[list[i].name].each do |e|
                  if e.start == skipped_exon.start and e.stop == skipped_exon.stop
                    found=true # found, therefore don't add the same exon to 'skipped' again
                  end
                end
              end

              if !found
                if !skipped.has_key?(list[i].name)
                  skipped[list[i].name]=[]
                end
                skipped_exon.cell = list[i].cell
                skipped[list[i].name] << skipped_exon
              end
            end
          end
          # intron retention
          unless exon_list_ret==[] 
            exon_list_ret.each do |retained_exon|
              found=false
              if retained.has_key?(list[i].name)
                retained[list[i].name].each do |e|
                  if e.start == retained_exon.start and e.stop == retained_exon.stop
                    found=true 
                  end
                end
              end
              if !found
                if !retained.has_key?(list[i].name)
                  retained[list[i].name]=[]
                end
                retained_exon.cell = list[i].cell
                retained[list[i].name] << retained_exon
              end
            end
          end

          # five prime
          unless exon_list_five_prime==[] 
            exon_list_five_prime.each do |five_exon|
              found=false
              if five_prime.has_key?(list[i].name)
                five_prime[list[i].name].each do |e|
                  if e.start == five_exon.start and e.stop == five_exon.stop
                    found=true 
                  end
                end
              end
              if !found
                if !five_prime.has_key?(list[i].name)
                  five_prime[list[i].name]=[]
                end
                five_exon.cell = list[i].cell
                five_prime[list[i].name] << five_exon
              end
            end
          end

          
          # three prime
          unless exon_list_three_prime==[] 
            exon_list_three_prime.each do |three_exon|
              found=false
              if three_prime.has_key?(list[i].name)
                three_prime[list[i].name].each do |e|
                  if e.start == three_exon.start and e.stop == three_exon.stop
                    found=true 
                  end
                end
              end
              if !found
                if !three_prime.has_key?(list[i].name)
                  three_prime[list[i].name]=[]
                end
                three_exon.cell = list[i].cell
                three_prime[list[i].name] << three_exon
              end
            end
          end
        end
      end
    end
  end

  skipped_output=""
  skipped.each_pair do |splicename, exon_list|
    if genes[splicename.split(".").first][splicename].strand == "+"
      skipped_output << "#{splicename}\t#{genes[splicename.split(".").first][splicename].cell}\t"
      exon_list.each do |exon|
        skipped_output << "#{exon.start-genes[splicename.split(".").first][splicename].start}..#{exon.stop-genes[splicename.split(".").first][splicename].start}\t"
      end
      skipped_output << "\n"
    else
      skipped_output << "#{splicename}\t#{genes[splicename.split(".").first][splicename].cell}\t"
      exon_list.each do |exon|
        skipped_output << "#{genes[splicename.split(".").first][splicename].stop-exon.stop}..#{genes[splicename.split(".").first][splicename].stop-exon.start}\t"
      end
      skipped_output << "\n"
    end
  end

  retained_output=""
  retained.each_pair do |splicename, exon_list|
    if genes[splicename.split(".").first][splicename].strand == "+"
      retained_output << "#{splicename}\t#{genes[splicename.split(".").first][splicename].cell}\t"
      exon_list.each do |exon|
        retained_output << "#{exon.start-genes[splicename.split(".").first][splicename].start}..#{exon.stop-genes[splicename.split(".").first][splicename].start}\t"
      end
      retained_output << "\n"
    else
      retained_output << "#{splicename}\t#{genes[splicename.split(".").first][splicename].cell}\t"
      exon_list.each do |exon|
        retained_output << "#{genes[splicename.split(".").first][splicename].stop-exon.stop}..#{genes[splicename.split(".").first][splicename].stop-exon.start}\t"
      end
      retained_output << "\n"
    end
  end

  five_prime_output=""
  five_prime.each_pair do |splicename, exon_list|
    if genes[splicename.split(".").first][splicename].strand == "+"
      five_prime_output << "#{splicename}\t#{genes[splicename.split(".").first][splicename].cell}\t"
      exon_list.each do |exon|
        five_prime_output << "#{exon.start-genes[splicename.split(".").first][splicename].start}..#{exon.stop-genes[splicename.split(".").first][splicename].start}\t"
      end
      five_prime_output << "\n"
    else
      five_prime_output << "#{splicename}\t#{genes[splicename.split(".").first][splicename].cell}\t"
      exon_list.each do |exon|
        five_prime_output << "#{genes[splicename.split(".").first][splicename].stop-exon.stop}..#{genes[splicename.split(".").first][splicename].stop-exon.start}\t"
      end
      five_prime_output << "\n"
    end
  end

  three_prime_output=""
  three_prime.each_pair do |splicename, exon_list|
    if genes[splicename.split(".").first][splicename].strand == "+"
      three_prime_output << "#{splicename}\t#{genes[splicename.split(".").first][splicename].cell}\t"
      exon_list.each do |exon|
        three_prime_output << "#{exon.start-genes[splicename.split(".").first][splicename].start}..#{exon.stop-genes[splicename.split(".").first][splicename].start}\t"
      end
      three_prime_output << "\n"
    else
      three_prime_output << "#{splicename}\t#{genes[splicename.split(".").first][splicename].cell}\t"
      exon_list.each do |exon|
        three_prime_output << "#{genes[splicename.split(".").first][splicename].stop-exon.stop}..#{genes[splicename.split(".").first][splicename].stop-exon.start}\t"
      end
      three_prime_output << "\n"
    end
  end

  File.open("#{opts.output}-skipped.txt", "w") {|file| file.write(skipped_output)}
  File.open("#{opts.output}-retained.txt", "w") {|file| file.write(retained_output)}
  File.open("#{opts.output}-five_prime.txt", "w") {|file| file.write(five_prime_output)}
  File.open("#{opts.output}-three_prime.txt", "w") {|file| file.write(three_prime_output)}

end

prefixes = ["genes-with-exon-skipping","genes-with-intron-retention","genes-with-five-prime-alternate-splicing","genes-with-three-prime-alternate-splicing"]

cmd = "ruby lib/results.rb -i #{opts.output}-skipped.txt > #{prefixes[0]}_#{opts.output}.txt"
`#{cmd}`
cmd = "ruby lib/results.rb -i #{opts.output}-retained.txt > #{prefixes[1]}_#{opts.output}.txt"
`#{cmd}`
cmd = "ruby lib/results.rb -i #{opts.output}-five_prime.txt > #{prefixes[2]}_#{opts.output}.txt"
`#{cmd}`
cmd = "ruby lib/results.rb -i #{opts.output}-three_prime.txt > #{prefixes[3]}_#{opts.output}.txt"
`#{cmd}`

results=Hash.new
(0..3).each do |i|
  cmd1 = "grep 35S #{prefixes[i]}_#{opts.output}.txt | wc -l"
  a= `#{cmd1}`.chomp.to_i
  cmd2 = "grep GDC #{prefixes[i]}_#{opts.output}.txt | wc -l"
  b= `#{cmd2}`.chomp.to_i
  cmd3 = "grep both #{prefixes[i]}_#{opts.output}.txt | wc -l"
  c= `#{cmd3}`.chomp.to_i
  cmd4 = "grep -v \"35S\\|GDC\\|both\" #{prefixes[i]}_#{opts.output}.txt | wc -l"
  d= `#{cmd4}`.chomp.to_i
  results[prefixes[i]]={:s35=>a, :gdc=>b, :both=>c, :neither=>d}
end

results.each_pair do |key, value|
  puts "#{key}"
  value.each do |key2, value2|
    puts "#{key2}\t#{value2}"
  end
end