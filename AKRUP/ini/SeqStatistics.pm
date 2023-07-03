package SeqStatistics;

use vars qw( @ISA $CODONS $synsites %synchanges );
use strict;
use Bio::Align::DNAStatistics;
use Bio::Root::Root;

@ISA = qw( Bio::Align::DNAStatistics );

my @t = split '', "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
$CODONS = $Bio::Align::DNAStatistics::CODONS;
$synsites = $Bio::Align::DNAStatistics::synsites;
%synchanges = %Bio::Align::DNAStatistics::synchanges;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
}

sub calc_K_values {
  my ( $self, $aln, $seq1_id, $seq2_id) = @_;
  $self->throw("Needs 3 arguments - an alignment object, and 2 sequence ids")
    if @_!= 4;
  $self->throw ("This calculation needs a Bio::Align::AlignI compatible object, not a [ " . ref($aln) . " ]object") unless $aln->isa('Bio::Align::AlignI');
  my @seqs = (
	      {id => $seq1_id, seq =>($aln->each_seq_with_id($seq1_id))[0]->seq},
	      {id => $seq2_id, seq =>($aln->each_seq_with_id($seq2_id))[0]->seq}
	     ) ;
  if (length($seqs[0]{'seq'}) != length($seqs[1]{'seq'})) {
    $self->throw(" aligned sequences must be of equal length!");
  }
  my $results = [];
  $self->_get_av_ds_dn_da(\@seqs, $results);
  return $results;
}

# Override the method of base class
sub _get_av_ds_dn_da {
  # takes array of hashes of sequence strings and ids   #
  my $self = shift;
  my $seq_ref = shift;
  my $result = shift if @_;
  my @caller = caller(1);
  my @seqarray = @$seq_ref;
  my $bootstrap_score_list;
  #for a multiple alignment considers all pairwise combinations#
  my %dsfor_average = (ds => [], dn => []); 
  for (my $i = 0; $i < scalar @seqarray; $i++) {
    for (my $j = $i +1; $j<scalar @seqarray; $j++ ) {
#      print "comparing $i and $j\n";
      if (length($seqarray[$i]{'seq'}) != length($seqarray[$j]{'seq'})) {
#	$self->warn(" aligned sequences must be of equal length!");
	$self->throw(" aligned sequences must be of equal length!");
	next;
      }
#      print "length is " , length($seqarray[$i]{'seq'}) , "\n";
#      my $syn_site_count = count_syn_sites($seqarray[$i]{'seq'}, $synsites);
#      my $syn_site_count2 = count_syn_sites($seqarray[$j]{'seq'}, $synsites);
#      print "syn 1 is $syn_site_count , syn 2 is $syn_site_count2\n";
      my ($syn_count, $non_syn_count, $aa_count, $r_p_count, $gap_cnt) = analyse_mutations($seqarray[$i]{'seq'}, $seqarray[$j]{'seq'});
#      print "syn sub is $syn_count, non-syn sub is $non_syn_count, aa sub is $aa_count, gap count is $gap_cnt\n";
      #get averages
#      my $av_s_site = ($syn_site_count + $syn_site_count2)/2;
      my $av_s_site = count_av_syn_sites($seqarray[$i]{'seq'}, $seqarray[$j]{'seq'}, $synsites);
      my $av_ns_syn_site = length($seqarray[$i]{'seq'}) - $gap_cnt- $av_s_site;
#      print "av syn site is $av_s_site, av non-syn site is $av_ns_syn_site\n";

      #calculate ps and pn  (p54)
      my $syn_prop = $syn_count / $av_s_site;
      my $nc_prop = $non_syn_count / $av_ns_syn_site;
#      print "syn prop is $syn_prop, non-syn prop is $nc_prop\n";

      #now use jukes/cantor to calculate D_s and D_n, would alter here if needed a different method
      my $d_syn = $self->jk($syn_prop);
      my $d_nc = $self->jk($nc_prop);
#      print "d_syn = $d_syn, d_nc = $d_nc\n";

      #get codon number of seq alignments (remove gaps)
      my $num_aa = (length($seqarray[$i]{'seq'}) - $gap_cnt) / 3;

      #calculate pa
      my $aa_prop = $aa_count / $num_aa;

      #calculate p1, p2, p3
      my $p1 = $r_p_count->[0] / $num_aa;
      my $p2 = $r_p_count->[1] / $num_aa;
      my $p3 = $r_p_count->[2] / $num_aa;

      #use Poisson correction to calcuate D_a
      my $d_aa = $self->pc($aa_prop);

      #JK calculation must succeed for continuation of calculation
      #ret_value = -1 if error
#      next unless $d_nc >= 0 && $d_syn >=0;

      push @{$dsfor_average{'ds'}}, $d_syn;
      push @{$dsfor_average{'dn'}}, $d_nc;

      #if not doing bootstrap, calculate the pairwise comparisin stats
#      if ($caller[3] =~ /calc_KaKs_pair/ || $caller[3] =~ /calc_all_KaKs_pairs/) {
      if ($caller[3] =~ /calc_K_values/) {#now calculate variances assuming large sample
#	my $d_syn_var =  Bio::Align::DNAStatistics::jk_var($syn_prop, length($seqarray[$i]{'seq'})  - $gap_cnt );
#	my $d_nc_var =  Bio::Align::DNAStatistics::jk_var($nc_prop, length ($seqarray[$i]{'seq'}) - $gap_cnt);
#	my $d_aa_var = pc_var($aa_prop, $num_aa);
	#now calculate z_value
#	print "d_syn_var is  $d_syn_var,and d_nc_var is $d_nc_var\n";
#	my $z = ($d_nc - $d_syn) / sqrt($d_syn_var + $d_nc_var);
	#	print "z is $z\n";
#	print "push results\n";
	push @$result , {S => $av_s_site, N=>$av_ns_syn_site, A=>$num_aa,
			 S_d => $syn_count, N_d =>$non_syn_count, N_a=>$aa_count,
			 P_s => $syn_prop, P_n=>$nc_prop, P_a=>$aa_prop,
			 D_s => @{$dsfor_average{'ds'}}[-1],
			 D_n => @{$dsfor_average{'dn'}}[-1], D_a => $d_aa,
#			 D_n_var =>$d_nc_var, D_s_var => $d_syn_var, D_a_var => $d_aa_var,
			 P1 => $p1, P2 => $p2, P3 => $p3,
			 Seq1 => $seqarray[$i]{'id'},
			 Seq2 => $seqarray[$j]{'id'},
#			 z_score => $z,
			};
	$self->warn (" number of mutations too small to justify normal test for  $seqarray[$i]{'id'} and $seqarray[$j]{'id'}\n- use Fisher's exact, or bootstrap a MSA")
	  if ($syn_count < 10 || $non_syn_count < 10 ) && $self->verbose > -1 ;
      }				#endif
    }
  }

  #warn of failure if no results hashes are present
  #will fail if Jukes Cantor has failed for all pairwise combinations
  #$self->warn("calculation failed!") if scalar @$result ==0;

  #return results unless bootstrapping
#  return $result if $caller[3]=~ /calc_all_KaKs/ || $caller[3] =~ /calc_KaKs_pair/; 
  return $result if $caller[3] =~ /calc_K_values/;
  #else if getting average for bootstrap
  return( Bio::Align::DNAStatistics::mean ($dsfor_average{'ds'}),
	  Bio::Align::DNAStatistics::mean ($dsfor_average{'dn'})) ;
}

# Override method of base class
sub analyse_mutations {
#  print ">>> in SeqStatistics::analysie_mutations\n";
  #compares 2 sequences to find the number of synonymous/non synonymous
  # mutations between them
  my ($seq1, $seq2) = @_;
  my %mutator = (2=> {0=>[[1,2], #codon positions to be altered depend on which is the same
			  [2,1]],
		      1=>[[0,2],
			  [2,0]],
		      2=>[[0,1],
			  [1,0]],	},
		 3=> [		#all need to be altered 
		      [0,1,2],
		      [1,0,2],
		      [0,2,1],
		      [1,2,0],
		      [2,0,1],
		      [2,1,0] ],
		);
  my $TOTAL = 0;		# total synonymous changes
  my $TOTAL_n = 0;      	# total non-synonymous changes
  my $TOTAL_a = 0;		# total amino-acid changes
  my @TOTAL_p = ( 0, 0, 0 );	# total uncleotide acid changes in 3 positions of codons
  my $gap_cnt = 0;

  my %input;
  my $seqlen = length($seq1);
  for (my $j=0; $j< $seqlen; $j+=3) {
    $input{'cod1'} = substr($seq1, $j,3);
    $input{'cod2'} = substr($seq2, $j,3);

    #ignore codon if beeing compared with gaps! 
    if ($input{'cod1'} =~ /\-/ || $input{'cod2'} =~ /\-/) {
      $gap_cnt += 3; #just increments once if there is a apair of gaps
      next;
    }

    $TOTAL_p[0]++ unless substr($input{'cod1'}, 0, 1) eq substr($input{'cod2'}, 0, 1);
    $TOTAL_p[1]++ unless substr($input{'cod1'}, 1, 1) eq substr($input{'cod2'}, 1, 1);
    $TOTAL_p[2]++ unless substr($input{'cod1'}, 2, 1) eq substr($input{'cod2'}, 2, 1);

    # different animo acids
    $TOTAL_a++ unless _is_same_aa(\%input);

    my ($diff_cnt, $same) = Bio::Align::DNAStatistics::count_diffs(\%input);

    #ignore if codons are identical
    next if $diff_cnt == 0 ;
    if ($diff_cnt == 1) {
      $TOTAL += $synchanges{$input{'cod1'}}{$input{'cod2'}};
      $TOTAL_n += 1 - $synchanges{$input{'cod1'}}{$input{'cod2'}};
      #print " \nfordiff is 1 , total now $TOTAL, total n now $TOTAL_n\n\n"
    } elsif ($diff_cnt ==2) {
      my $s_cnt = 0;
      my $n_cnt = 0;
      my $tot_muts = 4;
      #will stay 4 unless there are stop codons at intervening point
    OUTER:for my $perm (@{$mutator{'2'}{$same}}) {
	my $altered = $input{'cod1'};
	my $prev= $altered;
	#		print "$prev -> (", $t[$CODONS->{$altered}], ")";
	for my $mut_i (@$perm) { #index of codon mutated
	  substr($altered, $mut_i,1) = substr($input{'cod2'}, $mut_i, 1);
	  if ($t[$CODONS->{$altered}] eq '*') {
	    $tot_muts -=2;
	    #print "changes to stop codon!!\n";
	    next OUTER;
	  } else {
	    $s_cnt += $synchanges{$prev}{$altered};
	    #					print "$altered ->(", $t[$CODONS->{$altered}], ") ";
	  }
	  $prev = $altered;
	}
	#		print "\n";
      }
      if ($tot_muts != 0) {
	$TOTAL += ($s_cnt/($tot_muts/2));
	$TOTAL_n += ($tot_muts - $s_cnt)/ ($tot_muts / 2);
      }
 
    } elsif ($diff_cnt ==3 ) {
      my $s_cnt = 0;
      my $n_cnt = 0;
      my $tot_muts = 18;	#potential number  of mutations
    OUTER: for my $perm (@{$mutator{'3'}}) {
	my $altered = $input{'cod1'};
	my $prev= $altered;
	#	print "$prev -> (", $t[$CODONS->{$altered}], ")";
	for my $mut_i (@$perm) { #index of codon mutated
	  substr($altered, $mut_i,1) = substr($input{'cod2'}, $mut_i, 1);
	  if ($t[$CODONS->{$altered}] eq '*') {
	    $tot_muts -=3;
	    #	print "changes to stop codon!!\n";
	    next OUTER;
						
	  } else {
	    $s_cnt += $synchanges{$prev}{$altered};
	    #			print "$altered ->(", $t[$CODONS->{$altered}], ") ";
	  }
	  $prev = $altered;
	}
	#	print "\n";
			 
      }				#end OUTER loop
      #calculate number of synonymous/non synonymous mutations for that codon
      # and add to total
      if ($tot_muts != 0) {
	$TOTAL += ($s_cnt / ($tot_muts /3));
	$TOTAL_n += 3 - ($s_cnt / ($tot_muts /3));
      }
    }				#endif $diffcnt = 3
  }				#end of sequencetraversal
#  print " there are $TOTAL syn mutations and $TOTAL_n non-syn  mutations $TOTAL_a aa\n";
#  print "<<< get out of SeqStatistics::analysie_mutations\n";
  return ($TOTAL, $TOTAL_n, $TOTAL_a, \@TOTAL_p, $gap_cnt);
}

sub _is_same_aa {
  my ($r_input) = @_;
  my $cod1 = uc $r_input->{'cod1'};
  my $cod2 = uc $r_input->{'cod2'};
  return ( $t[$CODONS->{$cod1}] eq $t[$CODONS->{$cod2}] );
}

sub pc {
  my ($self, $p_a) = @_;
  return -1 * log(1 - $p_a);
}

sub pc_var {
  my ($p, $n) = @_;
  return $p / log((1 - $p) * $n);
}

sub jk {
	my ($self, $p) = @_;
	if ($p >= 0.75) {
#		$self->warn( " Jukes Cantor won't  work -too divergent!");
		return -1;
		}
	return -1 * (3/4) * (log(1 - (4/3) * $p));
}

# Override for debug
sub count_syn_sites {
    #counts the number of possible synonymous changes for sequence
    my ($seq, $synsite) = @_;
    die "not integral number of codons" if length($seq) % 3 != 0;
    my $S = 0;
    for (my $i = 0; $i< length($seq); $i+=3) {
	my $cod = substr($seq, $i, 3);
	next if $cod =~ /\-/;	#deal with alignment gaps
	$S +=  $synsite->{$cod}{'s'};
    }
    #print "S is $S\n";
    return $S;
}

sub count_av_syn_sites {
  # counts the average number of possible synonymous changes
  # for a pair of sequences
  my ($seq1, $seq2, $synsite) = @_;
  die "not integral number of codons" if length($seq1) % 3;
  my $S = 0;
  for (my $i = 0; $i < length($seq1); $i += 3) {
    my $cod1 = substr($seq1, $i, 3);
    my $cod2 = substr($seq2, $i, 3);
    next if ( $cod1 =~ /\-/ || $cod2 =~ /\-/ );
    $S += $synsite->{$cod1}{'s'} + $synsite->{$cod2}{'s'};
  }
  return $S / 2;
}


1;

