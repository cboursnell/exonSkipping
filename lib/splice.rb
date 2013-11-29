class Splice
  attr_accessor :name, :chromosome, :start, :stop, :strand, :cell, :exons, :s35, :gdc

  def initialize(name, chromosome, start, stop, strand)
    @name = name
    @chromosome = chromosome
    @start = start.to_i  # derived from the mRNA line in the gff file
    @stop = stop.to_i
    @strand = strand
    @cell = ""
    @exons = []
    @s35 = -1
    @gdc = -1 # -1 indicates these haven't been set yet...
  end

  def five_prime(splice)
    # warning! make sure you're looking at the right end
    # have to consider strand?
    #
    five = []
    if @strand != splice.strand
      abort "comparing 2 splice variants from different strands... i'm confused..."
    end
    if strand == "+"
      self.exons.each do |exon_a|
        overlap = false
        splice.exons.each do |exon_b|
          if exon_a.stop == exon_b.stop and exon_a.start!=exon_b.start 
            overlap = true
          end
        end
        if overlap
          five << exon_a.dup
        end
      end
    else
      self.exons.each do |exon_a|
        overlap = false
        splice.exons.each do |exon_b|
          if exon_a.start == exon_b.start and exon_a.stop!=exon_b.stop
            overlap = true
          end
        end
        if overlap
          five << exon_a.dup
        end
      end
    end
    return five
  end

  def three_prime(splice)
    three = []
    if @strand != splice.strand
      abort "comparing 2 splice variants from different strands... i'm confused..."
    end
    if strand == "+"
      self.exons.each do |exon_a|
        overlap = false
        splice.exons.each do |exon_b|
          if exon_a.start == exon_b.start and exon_a.stop!=exon_b.stop 
            overlap = true
          end
        end
        if overlap
          three << exon_a.dup
        end
      end
    else
      self.exons.each do |exon_a|
        overlap = false
        splice.exons.each do |exon_b|
          if exon_a.stop == exon_b.stop and exon_a.start!=exon_b.start
            overlap = true
          end
        end
        if overlap
          three << exon_a.dup
        end
      end
    end
    return three
  end

  def skipping(splice)
    # Exon Skipping:
    #
    # a)  --[===]---[===]--[===]---
    # b)  --[===]----------[===]---

    # exon skipping:
    # - exon should overlap with nothing on the other splice
    # - does it have to have an exon either side of it?
    skipped = []

    self.exons.each do |exon_a|
      overlap = false
      splice.exons.each do |exon_b|
        o =exon_a.overlap(exon_b)
        if o >= 0 and o <= 4
          overlap = true
        end
      end
      if !overlap
        skipped << exon_a.dup
      end
    end
    
    return skipped
  end

  def retention(splice)
    # Intron Retention
    #
    # a)  --[===========]---
    # b)  --[===]---[===]---

    # intron retention:
    # - should overlap with 1 exon and have the same start
    # - should overlap with another exon and have the same stop
    # - could overlap with another exon and have lower start and higher stop
    exons = []

    self.exons.each do |exon_a|
      overlap_1 = false
      overlap_2 = false
      splice.exons.each do |exon_b|
        o = exon_a.overlap(exon_b)
        if o == 1
          overlap_1=true
        elsif o == 2
          overlap_2=true
        end 
      end

      if overlap_1 and overlap_2
        exons << exon_a.dup
      end
    end
    return exons
  end

  def intersect(splice)
    intersect = self.dup

    intersect.exons.each do |exon_a|
      i=false
      splice.exons.each do |exon_b|
        ab = exon_a.overlap(exon_b)
        if ab==0     # perfect overlap
          # A |-----|     => I |-----|
          # B |-----|     
          i=true
        elsif ab==1
          # A    |-----|  => I    |--|
          # B |-----|     => B |-----|     
          exon_a.stop = exon_b.stop
          i=true
        elsif ab==2
          # A |-----|     => I    |--|
          # B    |-----|  => B    |-----|
          exon_a.start = exon_b.start
          i=true
        elsif ab==3
          # A    |---|    => I    |---|
          # B |---------| => B |---------|
          # nothing to do
          i=true
        elsif ab==4
          # A |---------|  => I    |---|   
          # B    |---|     => B    |---|       
          exon_a.start = exon_b.start
          exon_a.stop = exon_b.stop
          i=true
        elsif ab==5
          # A |---|       => I 
          # B       |---| => B       |---| 
        elsif ab==6
          # A       |---| => I
          # B |---|       => B |---|
        end
      end
      if i==false
        exon_a.start = -1
        exon_a.stop  = -1
      end
    end
    intersect.exons.reject! { |exon| exon.start == -1 and exon.stop == -1 }
    return intersect
  end

  def minus!(splice) # changes self
    # do self minus splice
    @exons.each do |exon_a|
      splice.exons.each do |exon_b|
        ab = exon_a.overlap(exon_b)
        if ab==0     # perfect overlap
          exon_a.start = -1
          exon_a.stop = -1
        elsif ab==1
          # A    |-----|  => A         |-|
          # B |-----|     => B  |-----|
          exon_a.start = exon_b.stop+1
        elsif ab==2
          # A |-----|     => A |-|
          # B    |-----|  => B    |-----|
          exon_a.stop = exon_b.start-1
        elsif ab==3
          # A    |---|    => A    
          # B |---------| => B |---------|
          exon_a.start = -1
          exon_a.stop  = -1
        elsif ab==4
          # A |---------|  => A |-|     |-| 
          # B    |---|     => B    |---|     
          exon_c = exon_a.dup
          exon_a.stop = exon_b.start - 1
          exon_c.start = exon_b.stop + 1
          @exons<<exon_c
        elsif ab==5 or ab==6
          # a and b don't overlap

        end

      end
    end
    # print "before:"
    # p @exons
    @exons.reject! { |exon| exon.start == -1 and exon.stop == -1 }
    # print "after:"
    # p @exons
    # puts "size: #{@exons.size}"
  end

  def length
    sum=0
    @exons.each do |exon|
      sum+=exon.length
    end
    sum
  end

  def to_s
    out = "#{@name}\t#{@chromosome}\t#{@start}..#{@stop} #{@strand}\t#{@cell}\n"
    @exons.each do |exon|
      out += "  #{exon}\n"
    end
    out
  end
end