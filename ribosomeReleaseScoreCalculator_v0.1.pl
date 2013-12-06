#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a perl script to reads calculate ribosome release score using the Ribo-Seq and RNA-Seq count data generated in pileupPerlStorableCounter
#
#	Input
#		--RNACountStorablePath=		file path[compulsory]; path of a list of perl storables generated from pileupPerlStorableCounter of the RNA-Seq sample;
#		--riboCountStorablePath=	file path[compulsory]; path of a list of perl storables generated from pileupPerlStorableCounter of the ribo-Seq sample;
#		--gffPath=					file path[compulsory]; path of the reference GFF for gene annotation;
#		--outDir=					directory path ['./BAMToReadEndPerlStorable/']; output directory;
#
#	v0.1
#		[Thu 12 Sep 2013 21:09:11 CEST] debut, and inherited from pileupPerlStorableCounterStatsCombiner_v0.1;
#
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-09-12 22:54]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/ribosomeReleaseScoreCalculator/v0.1/ribosomeReleaseScoreCalculator_v0.1.pl --RNACountStorablePath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS025/pileupPerlStorableCounterWithGetorf_UTR3Added/M5End.0.M3End.0/storable/normalizedCountHsh.pls --riboCountStorablePath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_pooled/pileupPerlStorableCounterWithGetorf_UTR3Added/M5End.0.M3End.0/storable/normalizedCountHsh.pls --gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/TbruceiTreu927_TriTrypDB-4.2_withGetorf.VC1.CN2.CP70.ML10.KO_yes.UTR3Added.gff --outDir=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/RRScore/BF/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/ribosomeReleaseScoreCalculator/v0.1/ribosomeReleaseScoreCalculator_v0.1.pl
#	--RNACountStorablePath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS028/pileupPerlStorableCounter_final_withNonAnnoORFuORF_UTR3Added/M5End.0.M3End.0/storable/normalizedCountHsh.pls
#	--riboCountStorablePath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS029_30_pooled/pileupPerlStorableCounter_final_withNonAnnoORFuORF_UTR3Added/M5End.0.M3End.0/storable/normalizedCountHsh.pls
#	--gffPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/TbruceiTreu927_TriTrypDB-4.2_withLongestTranscribedNonAnnoORFAnduORF.UTR3Added.gff
#	--outDir=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/RRScore/PF/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $global_scriptDir = dirname(rel2abs($0));
open DEBUGLOG, ">", "$global_scriptDir/debug.log.txt";
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|262, readParameters|435
#	secondaryDependOnSub: currentTime|214
#
#<section ID="startingTasks" num="0">
########################################################################## 
&printCMDLogOrFinishMessage("CMDLog");#->262
my ($RNACountStorablePath, $riboCountStorablePath, $gffPath, $outDir) = &readParameters();#->435
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
my $minCount = 0;
my $minLen = 30;
my $replaceZero = "yes";
my $paramTag = "MC.$minCount.ML.$minLen.RZ.$replaceZero";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="2">
my @mkDirAry;
my $resultDir = "$outDir/$paramTag"; push @mkDirAry, $resultDir;
my $resultLogDir = "$resultDir/log/"; push @mkDirAry, $resultLogDir;
foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_retrieveTheStorables
#	primaryDependOnSub: getGeneRngLength|231, readGFF_oneRNAPerGene|329
#	secondaryDependOnSub: currentTime|214
#
#<section ID="retrieveTheStorables" num="3">
my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);#->329
my ($geneRngLenHsh_ref) = &getGeneRngLength($geneInfoHsh_ref);#->231
my ($riboCountHsh_ref) = retrieve($riboCountStorablePath);
my ($RNACountHsh_ref) = retrieve($RNACountStorablePath);
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_getTheRatioAndScore
#	primaryDependOnSub: calculateRRScore|157
#	secondaryDependOnSub: reportStatus|468
#
#<section ID="getTheRatioAndScore" num="4">
my ($RRScoreDataHsh_ref) = &calculateRRScore($riboCountHsh_ref, $RNACountHsh_ref, $geneRngLenHsh_ref, $geneInfoHsh_ref, $minCount, $minLen, $replaceZero);#->157
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_printRRScoreXls
#	primaryDependOnSub: printRRScore|294
#	secondaryDependOnSub: >none
#
#<section ID="printRRScoreXls" num="5">
&printRRScore($RRScoreDataHsh_ref, $geneRngLenHsh_ref, $resultLogDir, $geneInfoHsh_ref);#->294
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|262
#	secondaryDependOnSub: currentTime|214
#
#<section ID="finishingTasks" num="6">
&printCMDLogOrFinishMessage("finishMessage");#->262
close DEBUGLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

sub calculateRRScore {
#....................................................................................................................................................#
#	dependOnSub: reportStatus|468
#	appearInSub: >none
#	primaryAppearInSection: 4_getTheRatioAndScore|122
#	secondaryAppearInSection: >none
#	input: $RNACountHsh_ref, $geneInfoHsh_ref, $geneRngLenHsh_ref, $minCount, $minLen, $replaceZero, $riboCountHsh_ref
#	output: $RRScoreDataHsh_ref
#	toCall: my ($RRScoreDataHsh_ref) = &calculateRRScore($riboCountHsh_ref, $RNACountHsh_ref, $geneRngLenHsh_ref, $geneInfoHsh_ref, $minCount, $minLen, $replaceZero);
#	calledInLine: 127
#....................................................................................................................................................#
	my ($riboCountHsh_ref, $RNACountHsh_ref, $geneRngLenHsh_ref, $geneInfoHsh_ref, $minCount, $minLen, $replaceZero) = @_;
	
	my %tmpCountHshRefHsh = ('ribo' => $riboCountHsh_ref, 'RNA' => $RNACountHsh_ref);
	my $RRScoreDataHsh_ref = {};
	my $totalGeneNum = (keys %{$geneInfoHsh_ref});
	my $geneProc = 0;
	
	my %tmpGeneIDHsh = ();
	foreach my $riboOrRNA (keys %tmpCountHshRefHsh) {
		foreach my $rngType (keys %{$tmpCountHshRefHsh{$riboOrRNA}}) {
			foreach my $geneID (keys %{$tmpCountHshRefHsh{$riboOrRNA}->{$rngType}}) {
				$tmpGeneIDHsh{$geneID}++;
			}
		}
	}
	
	foreach my $geneID (keys %tmpGeneIDHsh) {
		$geneProc++;
		&reportStatus("RRScore of $geneProc of $totalGeneNum calculated", 10, "\r");#->468
		next if $geneInfoHsh_ref->{$geneID}{'ctgry'} !~ m/mRNA/;
		$RRScoreDataHsh_ref->{$geneID}{'score'} = -1;
		$RRScoreDataHsh_ref->{$geneID}{'RNA_ratio'} = -1;
		$RRScoreDataHsh_ref->{$geneID}{'ribo_ratio'} = -1;
		my $validData = 'yes';
		
		foreach my $riboOrRNA (keys %tmpCountHshRefHsh) {
			foreach my $UTR3RngOrCDSRng (qw/UTR3Rng CDSRng/) {
				$RRScoreDataHsh_ref->{$geneID}{$riboOrRNA."_".$UTR3RngOrCDSRng} = $tmpCountHshRefHsh{$riboOrRNA}->{$UTR3RngOrCDSRng}{$geneID}{'count'}{'s'};
				$RRScoreDataHsh_ref->{$geneID}{$riboOrRNA."_".$UTR3RngOrCDSRng} = 0.5 if $RRScoreDataHsh_ref->{$geneID}{$riboOrRNA."_".$UTR3RngOrCDSRng} == 0 and $replaceZero eq 'yes';
				$validData = 'no' if $RRScoreDataHsh_ref->{$geneID}{$riboOrRNA."_".$UTR3RngOrCDSRng} < $minCount;
				$validData = 'no' if $geneRngLenHsh_ref->{$geneID}{$UTR3RngOrCDSRng} < $minLen;
			}
		}
	
		$RRScoreDataHsh_ref->{$geneID}{'validData'} = $validData;

		if ($validData eq 'yes') {
			foreach my $riboOrRNA (keys %tmpCountHshRefHsh) {
				$RRScoreDataHsh_ref->{$geneID}{$riboOrRNA.'_ratio'} = sprintf "%.5f", ($RRScoreDataHsh_ref->{$geneID}{$riboOrRNA."_CDSRng"}/$geneRngLenHsh_ref->{$geneID}{'CDSRng'})/($RRScoreDataHsh_ref->{$geneID}{$riboOrRNA."_UTR3Rng"}/$geneRngLenHsh_ref->{$geneID}{'UTR3Rng'});
			}
			$RRScoreDataHsh_ref->{$geneID}{'score'} = sprintf "%.5f", $RRScoreDataHsh_ref->{$geneID}{'ribo_ratio'}/$RRScoreDataHsh_ref->{$geneID}{'RNA_ratio'};
		}
	}

	return ($RRScoreDataHsh_ref);
}
sub currentTime {
#....................................................................................................................................................#
#	dependOnSub: >none
#	appearInSub: printCMDLogOrFinishMessage|262, readGFF_oneRNAPerGene|329, reportStatus|468
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|71, 3_retrieveTheStorables|109, 6_finishingTasks|142
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 281, 284, 289, 348, 483
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub getGeneRngLength {
#....................................................................................................................................................#
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 3_retrieveTheStorables|109
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: $geneRngLenHsh_ref
#	toCall: my ($geneRngLenHsh_ref) = &getGeneRngLength($geneInfoHsh_ref);
#	calledInLine: 115
#....................................................................................................................................................#
	my ($geneInfoHsh_ref) = @_;
	
	my $geneRngLenHsh_ref = {};
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		foreach my $rngType (qw/geneRng CDSRng UTR3Rng UTR5Rng exonRng/) {
			$geneRngLenHsh_ref->{$geneID}{$rngType} = -1;
			if ($geneInfoHsh_ref->{$geneID}{$rngType}) {
				
				#---get the rng length
				my $rngLength = 0;
				for (my $i=0; $i < $#{$geneInfoHsh_ref->{$geneID}{$rngType}}; $i += 2) {
					$rngLength += ${$geneInfoHsh_ref->{$geneID}{$rngType}}[$i+1] - ${$geneInfoHsh_ref->{$geneID}{$rngType}}[$i];
				}
				$geneRngLenHsh_ref->{$geneID}{$rngType} = $rngLength + 1;
			}
		}
	}

	return ($geneRngLenHsh_ref);
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	dependOnSub: currentTime|214
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|71, 6_finishingTasks|142
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 77, 147
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->214
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->214
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->214
		print "=========================================================================\n\n";
	}
}
sub printRRScore {
#....................................................................................................................................................#
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 5_printRRScoreXls|132
#	secondaryAppearInSection: >none
#	input: $RRScoreDataHsh_ref, $geneInfoHsh_ref, $geneRngLenHsh_ref, $resultLogDir
#	output: 
#	toCall: &printRRScore($RRScoreDataHsh_ref, $geneRngLenHsh_ref, $resultLogDir, $geneInfoHsh_ref);
#	calledInLine: 137
#....................................................................................................................................................#
	my ($RRScoreDataHsh_ref, $geneRngLenHsh_ref, $resultLogDir, $geneInfoHsh_ref) = @_;

	#---get ctgry
	open RRSCORELOG , ">", "$resultLogDir/RRScore.xls";

	foreach my $geneID (sort keys %{$RRScoreDataHsh_ref}) {
		my @outputAry = (qw/geneID description location CDSLength UTR3Length/);
		push @outputAry, $_ foreach sort keys %{$RRScoreDataHsh_ref->{$geneID}};
		print RRSCORELOG join "", ((join "\t", (@outputAry)), "\n");
		last;
	}
	
	foreach my $geneID (sort keys %{$RRScoreDataHsh_ref}) {
		my $description = $geneInfoHsh_ref->{$geneID}{'description'};
		my $location = $geneInfoHsh_ref->{$geneID}{'cntg'}.":".${$geneInfoHsh_ref->{$geneID}{'geneRng'}}[0]."-".${$geneInfoHsh_ref->{$geneID}{'geneRng'}}[1];
		my $CDSLength = $geneRngLenHsh_ref->{$geneID}{'CDSRng'};
		my $UTR3Length = $geneRngLenHsh_ref->{$geneID}{'UTR3Rng'};
		my @outputAry = ($geneID, $description, $location, $CDSLength, $UTR3Length);
		push @outputAry, $RRScoreDataHsh_ref->{$geneID}{$_} foreach sort keys %{$RRScoreDataHsh_ref->{$geneID}};
		print RRSCORELOG join "", ((join "\t", (@outputAry)), "\n");
	}
	close RRSCORELOG;
	return ();
}
sub readGFF_oneRNAPerGene {
#....................................................................................................................................................#
#	dependOnSub: currentTime|214
#	appearInSub: >none
#	primaryAppearInSection: 3_retrieveTheStorables|109
#	secondaryAppearInSection: >none
#	input: $gffPath
#	output: $geneInfoHsh_ref
#	toCall: my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);
#	calledInLine: 114
#....................................................................................................................................................#

	my ($gffPath) = @_;

	my $geneInfoHsh_ref = {};
	
	#---read the gff
	my $geneByRNAHsh_ref = {};

	open (GFF, $gffPath);
	print "[".&currentTime()."] Reading: $gffPath\n";#->214
	while (my $theLine = <GFF>) {

		chomp $theLine;
		
		last if $theLine =~ m/^##FASTA/;
		
		if ($theLine !~ m/^\#|^\@/ and $theLine !~ m/\tsupercontig\t/) {

			my ($seq, undef, $geneCategory, $featureStart, $featureEnd, undef, $geneStrd, undef, $dscrptns) = split (/\t/, $theLine);
			
			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));
			
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				
				$geneInfoHsh_ref->{$geneID}{'strnd'} = $geneStrd;
				$geneInfoHsh_ref->{$geneID}{'cntg'} = $seq;
				$geneInfoHsh_ref->{$geneID}{'description'} = uri_unescape($geneName);
				$geneInfoHsh_ref->{$geneID}{'description'} =~ s/\+/ /g;
				@{$geneInfoHsh_ref->{$geneID}{'geneRng'}} = ($featureStart, $featureEnd);

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'CDSRng'}}, ($featureStart, $featureEnd);
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'exonRng'}}, ($featureStart, $featureEnd);
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh_ref->{$RNAID} = $geneID;
				$geneInfoHsh_ref->{$geneID}{'ctgry'} = $geneCategory;
				@{$geneInfoHsh_ref->{$geneID}{'RNARng'}} = ($featureStart, $featureEnd);
				$geneInfoHsh_ref->{$geneID}{'RNAID'} = $RNAID;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close GFF;
	
	#---get the UTR if any
	my $minUTRLength = 10;
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {

		if (exists $geneInfoHsh_ref->{$geneID}{'CDSRng'}) {
			my $exonMin = min(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $exonMax = max(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $CDSMin = min(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});
			my $CDSMax = max(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});

			if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			} else {
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			}
		}
	}
	
	foreach my $rngType (qw/geneRng CDSRng exonRng UTR3Rng UTR5Rng/) {
		foreach my $geneID (keys %{$geneInfoHsh_ref}) {
			if ($geneInfoHsh_ref->{$geneID}{$rngType}) {
				@{$geneInfoHsh_ref->{$geneID}{$rngType}} = sort {$a <=> $b} @{$geneInfoHsh_ref->{$geneID}{$rngType}};
			}
		}
	}
	
	return ($geneInfoHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|71
#	secondaryAppearInSection: >none
#	input: none
#	output: $RNACountStorablePath, $gffPath, $outDir, $riboCountStorablePath
#	toCall: my ($RNACountStorablePath, $riboCountStorablePath, $gffPath, $outDir) = &readParameters();
#	calledInLine: 78
#....................................................................................................................................................#
	
	my ($RNACountStorablePath, $riboCountStorablePath, $gffPath, $outDir);
	
	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/ribosomeReleaseScoreCalculator/";
	
	GetOptions 	("RNACountStorablePath=s"  => \$RNACountStorablePath,
				 "riboCountStorablePath=s"  => \$riboCountStorablePath,
				 "gffPath=s"  => \$gffPath,
				 "outDir:s"  => \$outDir)

	or die		("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($RNACountStorablePath, $riboCountStorablePath, $gffPath) {
		die "Can't read $fileToCheck" if not -s $fileToCheck;
	}

	system "mkdir -p -m 777 $outDir/";
	
	return($RNACountStorablePath, $riboCountStorablePath, $gffPath, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	dependOnSub: currentTime|214
#	appearInSub: calculateRRScore|157
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_getTheRatioAndScore|122
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 185
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->214

	return ();
}

exit;
