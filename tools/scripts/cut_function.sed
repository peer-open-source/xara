#/^OPS_Export *void *\* *\nOPS_[A-z0-9]*(void)/s/.*/&/; tx; b; :x p; n; /^}/s/.*/&/p; Tx;
#/^void *\* *\nOPS_[A-z0-9]*(void)/s/.*/&/; tx; b; :x p; n; /^}/s/.*/&/p; Tx;
/^OPS_Export *void *\* *OPS_[A-z0-9]*(void)/s/.*/&/; tx; b; :x p; n; /^}/s/.*/&/p; Tx;
/^void *\* *OPS_[A-z0-9]*(void)/s/.*/&/; tx; b; :x p; n; /^}/s/.*/&/p; Tx;
/^OPS_Export *void *\* *OPS_[A-z0-9]*()/s/.*/&/; tx; b; :x p; n; /^}/s/.*/&/p; Tx;
/^void *\* *OPS_[A-z0-9]*()/s/.*/&/; tx; b; :x p; n; /^}/s/.*/&/p; Tx;
