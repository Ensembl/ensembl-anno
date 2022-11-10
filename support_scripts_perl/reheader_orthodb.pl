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

# OrthoDB = >6573_0:002ec4
my $outname=$ARGV[1];

open(IN,$ARGV[0]);
open(OUT, ">${outname}_Reheader.fa.tmp");

while(<IN>) {
  my $line = $_;

  if($line =~ /^>[0-9]+_[0-9]:\w+/) {
    my $header = $&;
    print OUT "$header\n";
  } else {
    print OUT ($line);
  }
}

close IN;
close OUT;
exit;
