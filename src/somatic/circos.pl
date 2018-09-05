#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
use Cwd;
use Getopt::Long;
sub usage
{
        print STDERR <<USAGE;
==============================================================================
Description    
Options
	-seg <s>:segment file,result of facets
        -tumor <s>: tumor prefix
        -outdir <s>:outdir of result
        -max <s>: max threshold for circos
        -min <s>: min threshold for circos
        -circos <s> :pathway of circos software
        -h|?|help : Show this help
==============================================================================
USAGE
}

my ($help,$max,$min,$outdir,$circos,$tumor,$seg,$wd);
GetOptions(
        "h|?|help"=>\$help,
        "outdir=s"=>\$wd,
        "seg=s"=>\$seg,
        "max=s"=>\$max,
        "min=s"=>\$min,
        "circos=s"=>\$circos,
        "tumor=s"   => \$tumor,
);
if( defined($help)||!defined ($seg) ||!defined ($tumor)){
        &usage;
        exit 0;
}
$circos ||= "/GPFS01/softwares/circos-0.69-5/bin/circos";
$max ||= 0.4;
$min ||= -0.4;
$wd ||= getcwd();
`awk -F "," '{if(\$5 > $max || \$5 < $min){print "hs"\$1"\t"\$10"\t"\$11"\t"\$5}}' $seg|sed '1d' >$wd/$tumor\_circos_input`;
open (OUT,">$wd\/$tumor\_circos_input.modified");
open (IN,"$wd\/$tumor\_circos_input");
my $bin = 100000000000000;
while (<IN>){
	chomp;
	my ($chr,$start,$end,$ratio) = split /\t+|\s+/,$_;
	if ($start =~ "e\+"){
		my @start = split /e\+/,$start;
		$start[1] =~ s/^0// if ($start[1] =~ /^0/);
		$start = $start[0] * (10**$start[1]);
	}
	if ($end =~ "e\+"){
                my @end=split /e\+/,$end;
                $end[1] =~ s/^0// if ($end[1] =~ /^0/);
	#	print "$end[1]\n";
                $end = $end[0] * (10**$end[1]);
		print "$end\n";
        }
	my $space = $end - $start;
	if ($space > $bin){
		my $end_tmp = 0;
		while ($end_tmp <= ($end - $bin)){
			$end_tmp = $start + $bin;
			print OUT "$chr\t$start\t$end_tmp\t$ratio\n";
			$start = $end_tmp + 1;
		}
		print OUT "$chr\t$start\t$end\t$ratio\n";
	}else {
		print OUT "$_\n";
	}
}
close OUT;
close IN;
my $ciros_used = "$wd/$tumor\_circos_input.modified";
open (CONF,">$wd/$tumor\_circos_input.config");
my $ab_config = "$wd/$tumor\_circos_input.config";
my $CON  =<<"END";
<<include etc/colors_fonts_patterns.conf>>
	
<ideogram>

<spacing>
# spacing between ideograms
default = 0.005r
break=15u
</spacing>
# ideogram position, thickness and fill
radius           = 0.65r
thickness        = 30p
fill             = yes
show_label=yes
#label_font=bold
label_radius=dims(ideogram,radius)+0.15r
label_size=30
label_parallel=yes
</ideogram>
	
<<include ticks.conf>>

<image>
<<include etc/image.conf>>
file*   = revise.png
radius* = 1000p
</image>

karyotype = /GPFS01/softwares/circos-0.69-5/etc/karyotype.human.txt 

chromosomes_units           = 1000000
chromosomes                 = -hsX;-hsY
chromosomes_display_default = yes

<plots>
# Data out of bounds should be hidden. Otherwise the
# # # default is to clip the data to range min/max.
range = hide

# scatter plot for values [-3,-0.6]
<plot>
type = scatter
file = $ciros_used
r0   = 0.65r
r1   = 0.77r
min  = -3
max  = $min
glyph = circle
glyph_size = 8
color = green

<axes>
<axis>
color     = lgreen
thickness = 2
spacing   = 0.3r
</axis>
</axes>

<backgrounds>
<background>
color = vlred_a5
</background>
</backgrounds>

<rules>
<rule>
condition  = 1
glyph_size = eval( 8 + 4*abs(var(value)))
flow       = continue
</rule>
<rule>
condition  = var(value) < -2
stroke_color = gray
stroke_thickness = 2
</rule>
</rules>
</plot>

# scatter plot for values [0.8,3]
<plot>
type = scatter
file = $ciros_used
r0   = 0.77r
r1   = 0.88r
min  = $max
max  = 3
glyph = circle
glyph_size = 8
color = red

<axes>
<axis>
color     = lred
thickness = 2
spacing   = 0.3r
</axis>
</axes>

<backgrounds>
<background>
color = vlred_a5
</background>
</backgrounds>
<rules>
<rule>
condition  = 1
glyph_size = eval( 8 + 4*abs(var(value)))
flow       = continue
</rule>
<rule>
condition    = var(value) < 3
stroke_color = black
stroke_thickness = 2
</rule>
</rules>

</plot>

####
type            = tile
layers_overflow = hide

<plot>
file        = $ciros_used
r1          = 0.96r
r0          = 0.91r
orientation = out

layers      = 15
margin      = 0.02u
thickness   = 4
padding     = 8

stroke_thickness = 1
stroke_color     = dblue
color            = blue
</plot>
</plots>

<<include etc/housekeeping.conf>>
END
print CONF "$CON";
close CONF;
`$circos -conf $ab_config -outputdir $wd -outputfile $tumor\_CNV_CIRCOS`;
chdir ("$wd");
