#!/usr/bin/env perl
# Copyright [2021] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use warnings;
use strict;
use feature 'say';


my $seq_hash = {};

my @all_lines = ();
my $count = 0;
open(IN,$ARGV[0]);
while(<IN>) {
  my $line = $_;
  chomp($line);
  if($line =~ /\>/) {
    if($count) {
      $count++;
      $all_lines[$count] = $line;
      $count++;
    } else {
      $all_lines[$count] = $line;
      $count++;
    }
  } else {
    if($all_lines[$count]) {
      $all_lines[$count] .= $line;
    } else {
      $all_lines[$count] = $line;
    }
  }

}
close IN;

for(my $i=0; $i<scalar(@all_lines); $i+= 2) {
  my $header = $all_lines[$i];
  my $seq = $all_lines[$i+1];
  unless(length($seq) >= 100) {
    next;
  }
  $seq_hash->{$seq} = $header;
}

foreach my $key (keys(%{$seq_hash})) {
  say $seq_hash->{$key};
  say $key;
}

exit;

