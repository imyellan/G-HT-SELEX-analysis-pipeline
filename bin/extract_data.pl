use strict;

if(@ARGV == 0){
die <<USAGE
usage:
perl extract_data.pl <chromDir> <bedfile> <chromcol> <poscol> <widthcol> <statcol> <outfile>
    
    <chromDir>: the directory containing the chromosome FastAs (as 
                typically downloaded from UCSC as "chromFa"): one FastA per
                chromosome, one sequence per file
    <bedfile>:  the file containing the peaks in tabular format, 
                e.g., bed, gff, narrowPeak
    <chromcol>: the column of <bedfile> containing the chromosome
    <poscol>:   the column of <bedfile> containing the position relative to
                the chromosome start
    <widthcol>: either i) <int> the column of <bedfile> containing the width of 
                          the peak or any region to be extracted 
                          symmetrically around <poscol>
                or    ii) f<int> a fixed width of all regions, centered at the
                          position given in <poscol>, e.g., f100 for a fixed
                          width of 100 bp
    <statcol>:  the column of <bedfile> containing the peak statistic
                or a similar measure of confidence
    <outfile>:  the path to the output file, written as FastA
USAGE
}


my $chrom = $ARGV[0];
my $bed = $ARGV[1];
my $chromcol = $ARGV[2]-1;
my $poscol = $ARGV[3]-1;
my $widthcol = $ARGV[4];
my $statcol = $ARGV[5]-1;
my $outfile = $ARGV[6];

my $sort = 1;

my $width = 0;

if($widthcol =~ /f([0-9]+)/i){
	$width = $1;
}else{
	$widthcol --;
}


sub loadSeq{
	my $prefix = shift;
	print $prefix," ";
	open(FA,$chrom."/".$prefix.".fa");
	my $head = "";
	my @lines = ();
	while(<FA>){
		chomp();
		if(/^>/){
			if($head){
				die "Chromosome FastA contained more than one sequence for ".$prefix."\n";
			}
			$head = $_;
		}else{
			push(@lines,lc($_));
		}
	}
	my $str = join("",@lines);
	print "loaded\n";
	return $str;
}



open(IN,$ARGV[1]);

my @lines = ();

while(<IN>){
	chomp();
	my @parts = split("\t",$_);
	$parts[$chromcol] =~ s/chr0/chr/g;
	my @vals = ($parts[$chromcol],$parts[$poscol],$parts[$statcol]);
	if($width){
		push(@vals,$width);
	}else{
		push(@vals,$parts[$widthcol]);
	}
	push(@lines,\@vals);
}

close(IN);
print "Read input file ".$bed."\n";


if($sort){

	@lines = sort { ${$a}[0] cmp ${$b}[0]  } @lines;

}

open(OUT,">".$outfile);

print "Extracting sequences...\n\n";

my $oldchr = "";
my $sequence = "";
for my $line (@lines){
	my @ar = @{$line};
	my $chr = $ar[0];
	unless($chr eq $oldchr){
		$sequence = loadSeq($chr);
	}
	$oldchr = $chr;
	my $w = $ar[3];
	if($w <= 0){
		print $w," -> next\n";
		next;
	}
	if($w % 2 == 0){
		$w = $w/2;
	}else{
		$w = ($w-1)/2;
	}

	my $start = $ar[1]-$w-1;

	my $head = "> chr: ".$chr."; start: ".$start."; peak: ".($ar[1]-$start)."; signal: ".$ar[2]."\n";
	my $curr = substr($sequence,$start,$ar[3]);
	if($curr =~ /[^ACGTacgt]/){
		print "Sequence for\n\t",substr($head,1),"omitted due to ambiguous nucleotides.\n\n";
	}else{
		print OUT $head,$curr,"\n";
	}
}

close(OUT);
print "\nDone.\n";