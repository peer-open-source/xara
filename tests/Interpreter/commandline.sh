
# run file, then drop into a REPL
opensees -i truss.tcl

# run file 'args.tcl'. Output should be 'a b c'
opensees args.tcl a b c

# This should be a file not found error
cat args.tcl | opensees a b c

#
cat args.tcl | opensees -c 'set argv newargs'

# run command 'set argv 3'. Output should be '3' on stdout
opensees -c 'puts "[set argv 3]"'

# run command 'set argv 3', then read file 'args.tcl', then drop into a REPL.
# The command should override argv and the output should be '3' on stdout
opensees -c 'set argv 3' args.tcl a b c

# this should print out "3"
echo 'set a 3' | python -m opensees -c 'puts $a'
