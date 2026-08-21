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

#!/usr/bin/perl -w
use strict;

my $current_extension = $ARGV[0];
my $new_extension = $ARGV[1];
my $pre_or_post = $ARGV[2];

my $usage = "quick_rename.pl 'current pre/ext' 'new pre/ext' 'prefix/extension'\n";
die "$usage" unless (@ARGV == 3);

my @Infiles;

print "New: $new_extension    Old: $current_extension\n";

if ($pre_or_post eq "prefix"){

	@Infiles = <$current_extension*>;
}
elsif ($pre_or_post eq "extension"){

	@Infiles = <*$current_extension>;
}
else{

	die "Option for prefix or postfix not specified correctly\n";
}

foreach my $file (@Infiles)
	{
		my $new_file_name = $file;

		$new_file_name =~ s/$current_extension/$new_extension/;
		`mv $file $new_file_name`;

	}

my $folder_contents=`ls -l ./*$new_extension*`;

print $folder_contents;

exit;
