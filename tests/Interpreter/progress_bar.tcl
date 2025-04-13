
set delay 80 ; # 0.08 seconds

progress create 200

for {set i 0} {$i <  50} {incr i} {
  after $delay;
  progress update ""
}

for {set i 0} {$i <  50} {incr i} {after $delay; progress update "Part I"}
for {set i 0} {$i <  50} {incr i} {after $delay; progress update "Part II"}
for {set i 0} {$i <  50} {incr i} {after $delay; progress update ""}

puts ""
