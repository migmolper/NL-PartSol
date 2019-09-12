

################################################################################
#
#    This file must be copied in the build directory of gidpost
#    distribution. After its execution, it will create two files:
#       gidpost.h gidpost.c
#    These are the only two files to be used in any c/c++ program
#
################################################################################

set header_h {
    
    #define HDF5
	
    #ifdef HDF5
    #include "hdf5.h"
    #endif
}

set header_c {
    #include <assert.h>
    #include <math.h>
    #include <string.h>
    #include <stdlib.h>
    #include <ctype.h>
    #include <stdarg.h>
    #include <stddef.h>
    #include <string.h>
    #include <stdio.h>
    #include <errno.h>
    #include <zlib.h>

    #include "gidpost.h"
    
    #ifdef WIN32
    #include <windows.h>
    #define snprintf _snprintf
    #else
    #include <pthread.h>
    #endif
}


set files_h_full [list gidpost.h hdf5c.h]

set files_h [list standard.h lookupa.h recycle.h hashtab.h gidpostHash.h gidpostInt.h gidpostHDF5.h]
set files_c [list gidpost.c     gidpostHash.c  gidpostInt.c  hdf5c.c    recycle.c   \
	gidpostfor.c  gidpostHDF5.c  hashtab.c           lookupa.c ]

set fout [open gidpost.h w]

puts $fout $header_h

foreach file $files_h_full { 
    puts $fout "\n"
    puts $fout "/*--------------------------------------------------------------------------------"
    puts $fout "    file: [file tail $file]"
    puts $fout "--------------------------------------------------------------------------------*/"
    puts $fout ""

    set fin [open [file join .. source $file] r]
    set data [read $fin]
    regsub -all -line {\s*#\s*include} $data "\n//#include" data
    puts $fout $data
    close $fin
}

close $fout


set fout [open gidpost.c w]

puts $fout $header_c

foreach file $files_h { 
    puts $fout "\n"
    puts $fout "/*--------------------------------------------------------------------------------"
    puts $fout "    file: [file tail $file]"
    puts $fout "--------------------------------------------------------------------------------*/"
    puts $fout ""

    set fin [open [file join .. source $file] r]
    set data [read $fin]
    regsub -all -line {\s*#\s*include} $data "\n//#include" data
    puts $fout $data
    close $fin
}

set nlines 0
foreach file $files_c {
    
    puts $fout "\n"
    puts $fout "/*--------------------------------------------------------------------------------"
    puts $fout "    file: [file tail $file]"
    puts $fout "--------------------------------------------------------------------------------*/"
    puts $fout ""
    set fin [open [file join .. source $file] r]
    set data [read $fin]
    regsub -all -line {\s*#\s*include} $data "\n//#include" data
    puts $fout $data
    incr nlines [regexp -all {\n} $data]
    close $fin
}

close $fout
exit

