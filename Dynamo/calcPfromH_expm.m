function calcPfromH_expm(t)

global OC;

OC.cache.P{t} = OC.config.expmFunc(-OC.seq.tau(t) * OC.cache.H{t});

