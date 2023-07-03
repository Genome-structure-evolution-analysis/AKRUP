
# https://metacpan.org/pod/BioPerl

use strict;

use Bio::SeqIO;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::Seq::EncodedSeq;
use Bio::AlignIO;

use FindBin qw($Bin);
use lib "$Bin";
use SeqStatistics;

my $dirname = $ARGV[0];
my $infile = "$ARGV[1]";
my $cdsseq = "$ARGV[2]";
my $id1 = "$ARGV[3]";
my $id2 = "$ARGV[4]";

# my $output = "./$ARGV[0]/$ARGV[5]";
# open(OUT, ">>" . $output) or die "cannot open outfile $output due to $!.\n";


my %dna_hash;
my $seq = Bio::SeqIO->new(-file => "./$dirname/$cdsseq", '-format' => 'fasta');
while (my $seqobj = $seq->next_seq()) {
    $dna_hash{$seqobj->id} = $seqobj;
}

my $is_prot_aln = Bio::AlignIO->new(-file => "./$dirname/$infile", -format => "CLUSTALW");
my $prot_aln = $is_prot_aln->next_aln();
my $dna_aln = &aa_to_dna_aln($prot_aln, \%dna_hash);
my $stats = new SeqStatistics;
my $result = $stats->calc_all_KaKs_pairs($dna_aln);
my ($Da, $Ds, $Dn, $N, $S, $S_d, $N_d);
for my $an (@$result) {
    for (sort keys %$an) {
        next if /Seq/;
        if ($_ eq "D_n") {$Dn = $an->{$_}};
        if ($_ eq "D_s") {$Ds = $an->{$_}};
        if ($_ eq "S_d") {$S_d = $an->{$_};}
        if ($_ eq "N_d") {$N_d = $an->{$_};}
        if ($_ eq "S") {$S = $an->{$_};}
        if ($_ eq "N") {$N = $an->{$_};}

    }
}

if ($Dn !~ /\d/) {$Dn = -2;}
if ($Ds !~ /\d/) {$Ds = -2;}

# return $id1 . "\t" . $id2 . "\t" . $Dn . "\t" . $Ds . "\n";
print "<";
print $id1 . "\t" . $id2 . "\t" . $Dn . "\t" . $Ds . "\n";
print ">";
# 获取返回值
