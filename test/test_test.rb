#!/usr/bin/env	ruby

require 'helper'

class TestExtractCellSpecificExons < Test::Unit::TestCase

  context "splice" do

    setup do
      name = "AT1G01640.1"
      chromosome = 1
      start = 100
      stop = 1000
      strand = "+"
      @a = Splice.new(name, chromosome, start, stop, strand)
      @a.exons << Exon.new(100, 200, "exon")
      @a.exons << Exon.new(300, 500, "exon")
      @a.exons << Exon.new(600, 800, "exon")
      @a.exons << Exon.new(900, 1000, "exon")

      name = "AT1G01640.2"
      chromosome = 1
      start = 100
      stop = 1000
      strand = "+"
      @b = Splice.new(name, chromosome, start, stop, strand)
      @b.exons << Exon.new(100, 200, "exon")
      @b.exons << Exon.new(300, 400, "exon")
      @b.exons << Exon.new(600, 800, "exon")
      @b.exons << Exon.new(900, 1000, "exon")
      
      name = "AT1G01640.3"
      chromosome = 1
      start = 100
      stop = 1000
      strand = "+"
      @c = Splice.new(name, chromosome, start, stop, strand)
      @c.exons << Exon.new(100, 200, "exon")
      @c.exons << Exon.new(300, 400, "exon")
      @c.exons << Exon.new(600, 800, "exon")
      @c.exons << Exon.new(920, 980, "exon")
      
      name = "AT1G01640.4"
      chromosome = 1
      start = 100
      stop = 1000
      strand = "+"
      @d = Splice.new(name, chromosome, start, stop, strand)
      @d.exons << Exon.new(100, 200, "exon")
      @d.exons << Exon.new(350, 400, "exon")
      @d.exons << Exon.new(500, 700, "exon")
      @d.exons << Exon.new(900, 960, "exon")    

      name = "AT2G02567.1"
      chromosome = 2
      start = 100
      stop = 1000
      strand = "+"
      @e = Splice.new(name, chromosome, start, stop, strand)
      @e.exons << Exon.new(100, 200, "exon")
      @e.exons << Exon.new(350, 400, "exon") # this exon is skipped in @f
      @e.exons << Exon.new(600, 800, "exon")
      @e.exons << Exon.new(920, 980, "exon")
      
      name = "AT2G02567.2"
      chromosome = 2
      start = 100
      stop = 1000
      strand = "+"
      @f = Splice.new(name, chromosome, start, stop, strand)
      @f.exons << Exon.new(100, 200, "exon")
      @f.exons << Exon.new(500, 700, "exon")
      @f.exons << Exon.new(900, 960, "exon")

      name = "AT3G06789.1"
      chromosome = 3
      start = 100
      stop = 1000
      strand = "+"
      @g = Splice.new(name, chromosome, start, stop, strand)
      @g.exons << Exon.new(100, 200, "exon")
      @g.exons << Exon.new(350, 800, "exon") # intron here is 
      @g.exons << Exon.new(920, 980, "exon")
      
      name = "AT3G06789.2"
      chromosome = 3
      start = 100
      stop = 1000
      strand = "+"
      @h = Splice.new(name, chromosome, start, stop, strand)
      @h.exons << Exon.new(100, 200, "exon")
      @h.exons << Exon.new(350, 400, "exon") 
      @h.exons << Exon.new(500, 800, "exon")
      @h.exons << Exon.new(900, 960, "exon")

      name = "AT4G02345.1"
      chromosome = 4
      start = 100
      stop = 1000
      strand = "-"
      @i = Splice.new(name, chromosome, start, stop, strand)
      @i.exons << Exon.new(100, 200, "exon")
      @i.exons << Exon.new(350, 400, "exon") 
      @i.exons << Exon.new(500, 800, "exon")
      @i.exons << Exon.new(900, 960, "exon")

      name = "AT4G02345.2"
      chromosome = 4
      start = 100
      stop = 1000
      strand = "-"
      @j = Splice.new(name, chromosome, start, stop, strand)
      @j.exons << Exon.new(100, 200, "exon")
      @j.exons << Exon.new(330, 400, "exon") 
      @j.exons << Exon.new(500, 820, "exon")
      @j.exons << Exon.new(900, 960, "exon")
    end

    should "detect intron retention" do
      exons = @g.retention(@h)
      unless exons==nil
        assert_equal exons.size, 1
        assert_equal exons[0].start, 350
        assert_equal exons[0].stop, 800
      end
    end

    should "detect exon skipping" do
      exons = @e.skipping(@f)
      unless exons==nil
        assert_equal exons.size, 1
        assert_equal exons[0].start, 350
        assert_equal exons[0].stop, 400
      end
    end

    should "test length of exon" do
      assert_equal @a.exons.first.length, 101
    end

    should "create intersection of c and d" do
      i = @c.intersect(@d)
      l = 101+51+101+41
      assert_equal i.length, l
    end

    should "overlap" do
      assert_equal @a.exons[0].overlap(@b.exons[0]) , 0
    end

    should "not overlap" do
      assert_equal @a.exons[0].overlap(@b.exons[1]), 5
    end

    should "partially overlap" do
      assert_equal @a.exons[1].overlap(@b.exons[1]), 1
    end

    should "still have exons left after doing a minus!" do
      @a.minus!(@b)
      assert @a.exons.length >0
    end

    should "not have any exons that overlap after doing a minus b" do
      @a.minus!(@b)
      found=false
      @a.exons do |exon_a|
        @b.exons do |exon_b|
          if exon_a.overlap(exon_b) <= 4
            found=true
          end
        end
      end
      assert_equal found, false
    end

    should "have one exon from 400 to 500 after a minus b" do
      @a.minus!(@b)
      assert_equal @a.exons[0].start, 401
      assert_equal @a.exons[0].stop, 500
    end

    should "have no exons left after doing b minus a" do
      @b.minus!(@a)
      assert_equal @b.exons.size, 0
    end

    should "create new exons" do
      @b.minus!(@c)
      assert_equal @b.exons.size, 2
    end

    should "find three prime splicing variant" do
      a = @a.three_prime(@b)
      assert a.length==1, "expected 1 but found #{a.length}\n#{a}"
    end

    should "find three prime splicing variant on reverse strand" do
      a = @i.three_prime(@j)
      assert a.length==1, "expected 1 but found #{a.length}\n#{a}"
    end

    should "find five prime splicing variant" do
      a = @i.five_prime(@j)
      assert a.length==1, "expected 1 but found #{a.length}\n#{a}"
    end
  end
end