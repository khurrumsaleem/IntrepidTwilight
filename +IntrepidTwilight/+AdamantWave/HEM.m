function hem = HEM()

    model          = IntrepidTwilight.AdamantWave.BasicModel();
    residual       = [];
    solver         = [];
    preconditioner = [];

    hem.model                = model            ;
    hem.discretization.space = 'Quasi2DUpwind'  ;
    hem.discretization.time  = 'ImplicitEuler'  ;
    hem.solver               = 'JFNK'           ;
    hem.preconditioner
    
    parameters = [];


end