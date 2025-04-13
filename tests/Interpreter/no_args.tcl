proc range args {
  foreach {start stop step} [switch -exact -- [llength $args] {
      1 {concat 0 $args 1}
      2 {concat   $args 1}
      3 {concat   $args  }
      default {error {wrong # of args: should be "range ?start? stop ?step?"}}
  }] break
  if {$step == 0} {error "cannot create a range when step == 0"}
  set range [list]
  while {$step > 0 ? $start < $stop : $stop < $start} {
      lappend range $start
      incr start $step
  }
  return $range
}

model basic 3 3

set cmds [info commands]
puts "$cmds"

foreach cmd $cmds {
  if {[lsearch -exact {fault quit exit} $cmd] >= 0} continue
  puts "$cmd"
  catch {$cmd} err
}

