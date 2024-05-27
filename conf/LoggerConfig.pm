package conf::LoggerConfig;

use strict;
use warnings;

use Log::Log4perl;

Log::Log4perl->init(\<<'EOT');
    log4perl.rootLogger = DEBUG, LOGFILE, SCREEN
    log4perl.appender.Screen = Log::Log4perl::Appender::Screen
    log4perl.appender.SCREEN.stderr = 1
    log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
    log4perl.appender.File = Log::Log4perl::Appender::File
    log4perl.appender.File.filename = /path/to/your/logfile.log
    log4perl.appender.File.layout = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.File.layout.ConversionPattern = %d %p %m %n

EOT

our $logger = Log::Log4perl->get_logger();

1;
