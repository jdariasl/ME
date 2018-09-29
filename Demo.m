function Demo()

%Demo estimation of Markov-based entropies

% Load Rossler attractor data
load('RosslerAttractor.mat');

% Estimates the Principal Curve using the method: Subspace Constrained Mean Shift, see:
% Ozertem, Umut; Erdogmus, Deniz, "Locally Defined Principal Curves and
% Surfaces", Journal of Machine Learning Research 12 (2011) 241-274

mCentroids = SCMS(x,1);

% Plot attractor and principal curve
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
view(axes1,[322.5 30]);
grid(axes1,'on');
hold(axes1,'all');

plot3(x(:,1),x(:,2),x(:,3),'Parent',axes1,...
    'Color',[0.831372559070587 0.815686285495758 0.7843137383461],...
    'DisplayName','R\"ossler Attractor');
hold on
plot3(mCentroids(:,1),mCentroids(:,2),mCentroids(:,3),'Parent',axes1,...
    'Marker','*','LineStyle','none',...
    'Color',[1 0 0],...
    'DisplayName','Principal Curve');
grid on;
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.151820349761527 0.636913875598085 0.267885532591415 ...
    0.0852272727272727],'Interpreter', 'latex');
%--------------------------------------------------------------------------
% Generatation a time series from a MIX process 
vSignal = MixF(0.5,1000);
iDim = 2;
iTao = 1;
Entropies=Markov_Entropies(vSignal,iDim,iTao);
figure();
plot(vSignal)
title('MIX Process')
%--------------------------------------------------------------------------
disp('----------------------------------------------');
fprintf('MIX process characterization \n')
disp('----------------------------------------------');
fprintf('Markov chain entropy-Shannon: %g \n', Entropies.OriginalDim.EMcs );
fprintf('Markov chain entropy-Renyi 2: %g \n',Entropies.OriginalDim.EMcr );
fprintf('Conditional HMP entropy (Normalized)-Shannon: %g \n', Entropies.OriginalDim.EhmmsN );
fprintf('Conditional HMP entropy (Normalized)-Renyi 2: %g \n', Entropies.OriginalDim.EhmmrN );
fprintf('Recurrence State Entropy-Shannon: %g \n', Entropies.OriginalDim.EMrds );
fprintf('Recurrence State Entropy-Renyi 2: %g \n',Entropies.OriginalDim.EMrdr );
disp('----------------------------------------------');
end


function vSerie = MixF(p,N)
rng('default');
X = sqrt(2)*sin(2*pi*(1:N)./12);
Y = 2*sqrt(3)*rand(1,N) - sqrt(3);
Z = rand(1,N);
Z = Z <= p;
vSerie = (1 - Z).*X + Z.*Y;

end

