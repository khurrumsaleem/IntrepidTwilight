clc();
clear();

hem = testHEM();

%   Adjust evolver time stuff
dt = 50;
n  = numel(dt);
q{1,n} = [];
Dq{1,n} = [];
t{1,n} = [];
tElapsed(1,n) = 0;

for k = 1:n
    
    %   Set evolver
    hem.set('evolver','time.span'         , [0,4000] );
    hem.set('evolver','time.step.maximum' ,    dt(k) );
    hem.set('evolver','time.step.minimum' ,   1E-12  );
    hem.set('evolver','time.step.goal'    ,   1      );
    hem.set('evolver','saveRate'          ,   0.5    );
    
    %   Run
    tic;
    hem.run();
    tElapsed(k) = toc;
    [q{k},Dq{k},t{k}] = hem.results();
    
end

t = {(0:0.5:4E3).'};

%%
%   Calculate pressures
mass{1,n}     = [];
energy{1,n}   = [];
volume{1,n}   = [];
momentum{1,n} = [];
T{1,n}        = [];
P{1,n}        = [];

for k = 1:n
    mass{k}     = q{k}( (1:14) + 0   ,:);
    energy{k}   = q{k}( (1:14) + 14  ,:);
    volume{k}   = q{k}( (1:14) + 14*2,:);
    momentum{k} = q{k}( (1:14) + 14*3,:);
    rho         = mass{k}   ./ volume{k};
    i           = energy{k} ./ mass{k}  ;
    T{k}        = Temperature(rho,i)    ;
    P{k}        = Pressure(rho,T{k})    ;
end

%% 

figure(6);
mask          = 1:numel(t{1});
ords          = {P{1}(2,mask),P{1}(4,mask),P{1}(6,mask),P{1}(8,mask),P{1}(10,mask),P{1}(12,mask),P{1}(14,mask)}   ;
n             = numel(ords)             ;
args{1,2*n}   = []                      ;
args(1:2:end) = {t{1}(mask)}            ;
args(2:2:end) = ords                    ;
h = plot(args{:},'LineWidth',2);
legend(strcat({'\Delta{t} = '},num2str(dt(:),'%7.2E'),{' [s]'},{';  t_{clock} = '},num2str(tElapsed.','%7.2f'),{' [s]'}));


% %%
% figure(2);
% Tmax = 1.0001*max(T{1}(:));
% Tmin = 0.9999*min(T{1}(:));
% Tplt = linspace(Tmin,Tmax,500);
% [Psat,rhoL,rhoG] = SaturationStateGivenTemperature(Tplt);
% iL = InternalEnergyOne(rhoL,Tplt);
% iG = InternalEnergyOne(rhoG,Tplt);
% for k = 1:numel(P{1}(2,:))
% %     plot(rhoL,iL,'k',mass{1}(:,k)./volume{1}(:,k),energy{1}(:,k)./mass{1}(:,k),'bo-');
%     plot(Tplt,Psat,'k',T{1}(:,k),P{1}(:,k),'bo-');
% %     axis([964.035,964.085,3.8464E5,3.848E5]);
%     title(['Time: ',num2str(t{1}(k),'%9.7e'),' [s]'])
%     drawnow();
%     pause(0.1);
% end
% 
%%
% figure(3);
% 
% for k = 1:numel(P{1}(2,:))
%     plot(mass{1}(:,k)./volume{1}(:,k),energy{1}(:,k)./mass{1}(:,k),'bo-');
%     title(['Time: ',num2str(t{1}(k),'%9.7e'),' [s]'])
%     xlim([40,75])
%     ylim([4.35E5,4.49E5]);
%     drawnow();
%     pause(0.01);
% end
