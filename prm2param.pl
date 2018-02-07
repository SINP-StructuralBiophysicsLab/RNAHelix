#!/usr/bin/perl 
  ##########################################################################
  # *.prm is the output of nuparm                                          #
  # this script takes input from *.prm file and prepares parameter.loc     #
  # and sequence.str file                                                  #
  # parameter.loc is the input of nucgen/RNAhelix                          #
  # sequence.str is required for the minimization of renegerated backbone  #
  #                                                                        #
  ##########################################################################
use strict;
use warnings;

my $usage="
  ##########################################################################
    # *.prm is the output of nuparm                                          #
      # this script takes input from *.prm file and prepares parameter.loc     #
        # and sequence.str file                                                  #
          # parameter.loc is the input of nucgen/RNAhelix                          #
            # sequence.str is required for the minimization of renegerated backbone  #
              #                                                                        #
                ##########################################################################
                
[USAGE]
\n$0 file.prm\n\n";
my @lcp;
my @blp;
my @hed;
my @chna;
my @chnb;
my $tests=shift ||die $usage;
open KK,$tests ||die $!;
open LL,">parameter.loc" ||die $!;
open MM,">Sequence.str" ||die $!;
while(<KK>){
	if($_=~ /^LC/){
		my @alc=split /\s+/,$_;
		my $chk=substr($_,19,47)||'   0.00    0.00    0.00    0.00    0.00    0.00';
		push @lcp,$chk;
	}
	if($_=~ /BL/){
		my @abl=split /\s+/,$_;
		my $chk=substr($_,18,48);
		my $nwchk="$abl[2] $abl[10]$abl[11] ";
		my @alseq=split ':',$abl[2];
		push @chna,$alseq[0];
		push @chnb,$alseq[1];
		push @hed,$nwchk;
		push @blp,$chk;
	}
}
close KK;
my $len=scalar @blp;
print LL"  $len   0   0   $tests\n";
print MM"*FILENAME: Sequence.str
*PURPOSE: This file helps in reading the sequence information dynamically
*AUTHOR: Debasish Mukherjee, Dhananjay Bhattacharyya, SINP (May 16, 2016)
*

read sequ card
*Sequence for Chain A \n* \n$len\n";
for(my $i=0;$i<$len;$i++){
	print  LL"$hed[$i]  $lcp[$i]$blp[$i]\n";
	if($chna[$i] eq "A"){
		print MM"ADE ";
	}
	if($chna[$i] eq "U"){
		print MM"URA ";
	}
	if($chna[$i] eq "T"){
		print MM"THY ";
	}
	if($chna[$i] eq "G"){
		print MM"GUA ";
	}
	if($chna[$i] eq "C"){
		print MM"CYT ";
	}
}
close LL;

print MM"\n
gene RNA setup first 5ter last 3ter\n
read sequ card
*Sequence for Chain B \n* \n$len\n";
for(my $i=$len-1;$i>=0;$i--){
        if($chnb[$i] eq "A"){
                print MM"ADE ";
        }
        if($chnb[$i] eq "U"){
                print MM"URA ";
        }
        if($chnb[$i] eq "T"){
                print MM"THY ";
        }
        if($chnb[$i] eq "G"){
                print MM"GUA ";
        }
        if($chnb[$i] eq "C"){
                print MM"CYT ";
        }
}
print MM"\ngene RNA2 setup first 5ter last 3ter\n
join RNA RNA2 renum

return";
