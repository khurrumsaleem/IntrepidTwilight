function hem = HEM()

    %   Inherit
    hem = IntrepidTwilight.executive.Simulation();

    %   Set module and remove user's ability to do likewise
    hem.setBuildModuleName(struct('spacediscretization','AdamantWave','model','AdamantWave'));
    hem = rmfield(hem,'setBuildModuleName');
    
    
    %   Make name setters and remove user's ability to do likewise
    hem = hem.makeChoice(hem,'spacediscretization',IntrepidTwilight.AdamantWave.spacediscretizations());
    hem = hem.makeChoice(hem,'timediscretization',IntrepidTwilight.TransientStride.timediscretizations());
    hem = hem.makeChoice(hem,'residual',{'residual'});
    hem = hem.makeChoice(hem,'preconditioner',{'Preconditioner'});
    hem = hem.makeChoice(hem,'solver',IntrepidTwilight.TenaciousReduction.solvers());
    hem = hem.makeChoice(hem,'evolver',{'evolver'});
    hem = rmfield(hem,'makeChoice');

    
    %   Set Defaults
    hem.choose.spacediscretization.Quasi2DUpwind();
    hem.choose.timediscretization.ImplicitEuler();
    hem.choose.residual.residual();
    hem.choose.preconditioner.Preconditioner();
    hem.choose.solver.JFNK();
    hem.choose.evolver.evolver();

    %   Bind basic AdamantWave module
    hem.bind(IntrepidTwilight.AdamantWave.BasicModel());



end
% 
% function component = buildComponent(module,object,varargin)
%     component = IntrepidTwilight.executive.build(module,object,varargin{:});
% end
