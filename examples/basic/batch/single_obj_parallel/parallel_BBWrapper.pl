#!/usr/bin/perl
# use strict;
use warnings;
 
use Data::Dumper;
 
use threads;
use threads::shared;
use Thread::Semaphore;
use Config;
 

my $numberParallelJobs:shared=4;
my $semaphoreBBEval = Thread::Semaphore->new($numberParallelJobs);

## The blackbox executable
my $OSname:shared = "$Config{osname}";
my $BBdotEXE:shared="./bb.exe";
if ( $OSname eq "MSWin32" ) {
    $BBdotEXE="bb.exe"; 
}





if ( ! exists $ARGV[0]) {
	$nameInputFile = "";     
} else {
	$nameInputFile = $ARGV[0];
}

sub OneEvalBB($$$$$){
	my $x = shift;
	my $index = shift;
	my $completedEval = shift;
	my $sema_ref = shift;
	my $output = shift;
	
	
 
	my $input_file_name="";
	#if ($index > 3) {   # arbitrarily reject evaluation with index > 3!!!!! THIS IS NOT A GOOF THING TO DO
	#	$$output="REJECT\n";
    #} else {	
		
    	# write single evaluation input file
		$input_file_name="x.$$.$index.txt";
		open(SINGLE_INPUT_FILE,"> $input_file_name");
		print SINGLE_INPUT_FILE "$x \n";
		close(SINGLE_INPUT_FILE);

	
		# start single evaluation and get output
		open(OUTPUT_I,"$BBdotEXE $input_file_name |") ; # or die "Can't run blackbox executable: $!\n";
		$$output=<OUTPUT_I>;
		close(OUTPUT_I);
	#}


	$$completedEval++;

	# on a une place de libre. Ne pas oublier de libérer le sémaphore même en cas d'erreur
	$$sema_ref->up();

	if ( -e $input_file_name) {
		unlink $input_file_name ;
	}
 
	return;
}

open(LIST_INPUT_FILE,"$ARGV[0]") or die "Can't open input file: $ARGV[0]\n";
@X = <LIST_INPUT_FILE>;
my $started = 0;
my $completedBBEval:shared =0;
my @output:shared;

my @thrList= ();



# démarre tous les jobs 
while ( $started < scalar @X ){
	
	my $x=$X[$started];
	
	$started++;
	
	# avons nous une place de libre ?
	$semaphoreBBEval->down();
 
	# si le sémaphore est a 0, le processus principal va se bloquer en attendant une nouvelle place
	my $thr= threads->create("OneEvalBB", (
		$x,
		$started,
		\$completedBBEval,
		\$semaphoreBBEval,
		\$output[$started-1]
		)
	);
	
	push(@thrList,$thr);

	# détache le job du thread principal, rend la main au thread principal
	$thr->detach();
	
}

 
# attend les derniers jobs
while ( @thrList ){
	sleep(0.01);

	my $r=0;  
	foreach $thr (@thrList)
	{
		if ( $thr->is_running() ) {
			$r=1;
			last;
		}
	}
	last if ($r==0);
} 
print "@output";
