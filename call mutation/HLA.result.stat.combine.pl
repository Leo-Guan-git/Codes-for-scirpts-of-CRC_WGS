use warnings;
use strict;

my $Pfile = "$ARGV[0]/winners.hla.txt";
open O,">$ARGV[1]" or die "$!";

my $num = 1;
		open I1,$Pfile or die "$!";
		my ($PA,$PB,$PC);
		foreach my $line(<I1>){
			chomp $line;
                	my @tmp = split (/\t/,$line);
                	my @tmp2 = split(/\_/,$tmp[1]);
                	my @tmp3 = split(/\_/,$tmp[2]);
                	if($tmp2[1] eq "a"){
                        	$PA = "HLA-A$tmp2[2]\:$tmp2[3]\nHLA-A$tmp3[2]\:$tmp3[3]";
                	}
                	if($tmp2[1] eq "b"){
                        	$PB = "HLA-B$tmp2[2]\:$tmp2[3]\nHLA-B$tmp3[2]\:$tmp3[3]";
                	}
                	if($tmp2[1] eq "c"){
                        	$PC = "HLA-C$tmp2[2]\:$tmp2[3]\nHLA-C$tmp3[2]\:$tmp3[3]";
                	}
        	}	
		print O "$PA\n$PB\n$PC\n";

