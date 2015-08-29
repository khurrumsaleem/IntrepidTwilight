clc();
clear();


sim = IntrepidTwilight.executive.Simulation2();
sim.model.set('module','AdamantWave');
sim.model.set('name','BasicModel');
sim.spacediscretization.set('module','AdamantWave');
sim.spacediscretization.set('name','Quasi2DUpwind');
sim.build();