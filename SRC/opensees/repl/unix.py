file_util_commands = r"""
proc echo {args} {puts "$args"}

proc ls { {dir .} } { 
        # Get the current year, because the date format depends on it.

        set thisYear [clock format [clock seconds] -format %Y]
    
        # Walk the files in the given directory, accumulating lines
        # in $retval

        set retval {}
        set sep {}
        # In Tcl older than 8.3 use 'glob [file join $dir *]'
        foreach fileName [lsort [glob -dir $dir *]] {
    
            append retval $sep
            set sep \n
    
            # Get status of the file

            #file stat $fileName stat
            # use 'file lstat' instead: if the file is a symbolic link we don't want info about its target
            file lstat $fileName stat
    
            # Put in one character for file type.  Use - for a plain file.

            set type -
            if { [info exists stat(type)]
                 && [string compare file $stat(type)] } {
                set type [string index $stat(type) 0]
            }
            append retval $type
    
            # Decode $stat(mode) into permissions the way that ls does it.

            foreach { mask pairs } {
                00400 { 00400 r }
                00200 { 00200 w }
                04100 { 04100 s 04000 S 00100 x }
                00040 { 00040 r }
                00020 { 00020 w }
                02010 { 02010 s 02000 S 00010 x }
                00004 { 00004 r }
                00002 { 00002 w }
                01001 { 01001 t 01000 T 00001 x }
            } {
                set value [expr $stat(mode) & $mask]
                set bit -
                foreach { x b } $pairs {
                    if { $value == $x } {
                        set bit $b
                    }
                }
                append retval $bit
            }
    
            # Put in link count, user ID, and size.  Note that the UID
            # will be numeric.  If you know how to back-translate this
            # from Tcl, please feel free to edit it in!

            # LV writes - use file userid and file groupid to convert the numbers back to names.
            #   I don't know what version of Tcl added those commands...

            append retval [format %4d $stat(nlink)] { }

            array set attribs [file attributes $fileName]
            if {[info exists attribs(-owner)]} {
                append retval [format %-8s $attribs(-owner)]
                append retval [format %-8s $attribs(-group)] 
            } else {
                append retval [format %8d $stat(uid)]
                append retval [format %8d $stat(gid)]
            }
            append retval [format %9d $stat(size)]

            # Put in the date.  The current year is formatted differently
            # from prior years.

            set year [clock format $stat(mtime) -format "%Y"]
            if { $year == $thisYear } {
                set modified [clock format $stat(mtime) -format "%h %e %H:%M"]
            } else {
                set modified [clock format $stat(mtime) -format "%h %e  %Y"]
            }
            # glennj: see note below
            append retval { } $modified { }

            # Finally, put in the file name, stripping off the directory.

            append retval [file tail $fileName]

            if {[string compare $stat(type) link] == 0} {
                append retval " -> [file readlink $fileName]"
            }    

            unset stat attribs
    
        }
    
        return $retval
    
}
"""
