clc();
clear();

y0 = 1          ;
f  = @(t) cos(t);

sd = IntrepidTwilight.AdamantWave.SpaceDiscretization();
sd.rhs    = f;
sd.getAll = @() y0;

evolver = makeEvolver(sd);