&pcm
    NVARS      = 8
    var_names  = 'N', 'mu', 'sig', 'kap',   'V',  'T',     'P',    'ac',
    log_var    =   T,    T,     F,    F,      T,    F,       F,       F,
    lower_bnds = 1.0,  -3.0,     1.2,     0.0, -2.0, 240.,  50000., 0.1, 
    upper_bnds = 4.0,  -1.0,     3.0,     1.2,  1.0, 310., 105000., 1.0,
    PRNT_DBG   = T,
/

