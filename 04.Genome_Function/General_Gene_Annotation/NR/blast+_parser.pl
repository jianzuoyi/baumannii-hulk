#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: blast+_pahse.pl
#
#        USAGE: ./blast+_pahse.pl  
#
#  DESCRIPTION: This script used to extract information from tblastx(blast+) pairwased output
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Liwujiao (), hnnd059@gmail.com
# ORGANIZATION: Sichuan Key Laboratory of Conservation Biology on Endangered Wildlife
#      VERSION: 1.0
#      CREATED: 11/28/2012 06:35:46 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Long;
use Data::Dumper;

my ($Nohead,$Tophit,$Topmatch,$Eval);
my ($Verbose,$Help);
GetOptions(
  "nohead"    =>\$Nohead,
  "tophit:i"  =>\$Tophit,
  "topmatch:i"=>\$Topmatch,
  "eval:f"    =>\$Eval,
  "verbose"   =>\$Verbose,
  "help"      =>\$Help
);

if(@ARGV==0 or $Help){
  print "Usage: perl $0 blast_out > out\n";exit(0);
}

my $blast_file=shift @ARGV;
###convert blast+ raw result to tabular format
&parse_blast($blast_file,$Tophit,$Topmatch,$Eval,$Nohead);

##############Subroutines##################
sub parse_blast{
  my($file,$tophit,$topmatch,$eval,$nohead)=@_;
  open BLAST,$file or die "Could not open blast file:$!\n";
  print"Query_id\tQuery_length\tQuery_start\tQuery_end\tSubject_id\tSubject_length\tSubject_start\tSubject_end\tFrame\t",
  "Identity\tPositive\tGap\tAlign_length\tScore\tE_value\tQuery_annotation\tSubject_annotation\n" unless(defined $nohead);
  
  my $new_line="Effective search";
  $/=$new_line;
  while(<BLAST>){
    chomp;
	next if(/\* No hits found \*/);
	my @cycle=split (/\n>/,$_);#利用"\n>"将Query序列和每一条subject序列比对结果信息分开
	
	my ($pointer,$query,$query_len,$subject,$subject_len,$query_annotation,$subject_annotation,$frame);
    if($cycle[0]=~/Query=\s+(\S+)\s+(.*?)Length=(\d+)/s){
	  $query=$1;
	  $query_annotation=$2;
	  $query_len=$3;
	  $query_annotation =~ s/\s+/ /g; ##处置空白符
	  $query_annotation="--" if(!$query_annotation || $query_annotation eq " ");	#信息为空时用"--"代替
	}#提取Query id,Query length,Query annotation信息

	shift @cycle;# 去掉Query相关信息
	for(my $i=0;$i<@cycle;$i++){
	  last if ((defined $tophit) && $i>$tophit-1);
      if ($cycle[$i]=~/(\S+)\s+(.*?)\s*Length=(\d+)/s){
	    $subject=$1;
		$subject_annotation=$2;
		$subject_len=$3;
		$subject_annotation=~ s/\s+/ /g; ##处置空白符
		$subject_annotation="--" if(!$subject_annotation || $subject_annotation eq " ");	#信息为空时用"--"代替
		}	#提取Subject id,Subject length,Subject annotation信息

		my @cycle_inner=split (/Score =/,$cycle[$i]);	#分开同一个query和同一个subject的多个比对结果
		shift @cycle_inner;# 去掉第一个不含信息的元素
		for (my $j=0;$j<@cycle_inner;$j++){
		  last if((defined $topmatch) && $j>$topmatch-1);
		  $pointer->[$i][$j]{score}=$1 if( $cycle_inner[$j]=~/([\d\,\.]+)\s*bits\s+\(/);
		  $pointer->[$i][$j]{e_value}=$1 if ($cycle_inner[$j]=~/Expect\(*\d*\)*\s*=\s*(\S+),?/);
		  if ($cycle_inner[$j]=~/Identities\s*=\s*[\d\,\.]+\/([\d\,\.]+)\s*\((\S+)\%\)/s) {
		    $pointer->[$i][$j]{align_len}=$1;
			$pointer->[$i][$j]{identity}=$2/100;
    	  }	#提取Score,E value,Identity,Align_len信息
			last if((defined $eval) && $pointer->[$i][$j]{e_value}>$eval);

			$pointer->[$i][$j]{positive}=($cycle_inner[$j]=~/Positives\s*=\s*\S+\s*\((\S+)\%\)/s)? $1/100 : "--";
			#BLASTN结果文件中无Positive信息,其他几种结果文件都有
			exit if (!$pointer->[$i][$j]{positive});
			$pointer->[$i][$j]{gap}=($cycle_inner[$j]=~/Gaps\s*=\s*\S+\s*\((\S+)\%\)/)? $1/100 : 0;
			#Gap信息不是每一组里都有,单独考虑

			$pointer->[$i][$j]{q_start}=$1 if($cycle_inner[$j]=~/Query\s*([\d\,\.]+)\s*/);
     		$pointer->[$i][$j]{s_start}=$1 if($cycle_inner[$j]=~/Sbjct\s*([\d\,\.]+)\s*/);
		  $cycle_inner[$j]=~s/\n\s*Database:\s*.+?$//s;	#最后一组情况中需删除文件末尾的某些多余信息
		  $cycle_inner[$j]=~s/\s*?Lambda.+$//s;
		  $cycle_inner[$j]=~s/\s+$//s;
		  
		  if ($cycle_inner[$j]=~/Query\s*[\d\,\.]+\s*\D+([\d\,\.]+)\D+?Sbjct\s*[\d\,\.]+\s*\D+([\d\,\.]+)$/s){ # blast+ Query:\s*
			  $pointer->[$i][$j]{q_end}=$1;
			  $pointer->[$i][$j]{s_end}=$2;
			  #$pointer->[$i][$j]{q_start} =~ s/,//g; 
			  #$pointer->[$i][$j]{s_start} =~ s/,//g; 
			  #$pointer->[$i][$j]{q_end} =~ s/,//g; 
			  #$pointer->[$i][$j]{s_end} =~ s/,//g; 
		  }
		  $pointer->[$i][$j]{frame}=( ($cycle_inner[$j]=~/Frame =\s*(\S+)/)? $1 : "--");# add the frame information if exist

		
		  print "$query\t$query_len\t$pointer->[$i][$j]{q_start}\t$pointer->[$i][$j]{q_end}\t",
		  "$subject\t$subject_len\t$pointer->[$i][$j]{s_start}\t$pointer->[$i][$j]{s_end}\t$pointer->[$i][$j]{frame}\t",
		  "$pointer->[$i][$j]{identity}\t$pointer->[$i][$j]{positive}\t$pointer->[$i][$j]{gap}\t$pointer->[$i][$j]{align_len}\t",
		  "$pointer->[$i][$j]{score}\t$pointer->[$i][$j]{e_value}\t$query_annotation\t$subject_annotation\n";
      }
    }
  }
  $/="\n";
  close(BLAST);
}
