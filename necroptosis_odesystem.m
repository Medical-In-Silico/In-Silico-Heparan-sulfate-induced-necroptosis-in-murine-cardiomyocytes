function [ ydot ] = nekroptosis_odesystem(t, y, x, hs)
%% Basic Model of Nekrose
% Inhibition form: s * s / (w2*s + w1*y)
% steady-state: linear
% HS: prescribed function
%% required inputs
% (t, y, x, hs)
%% required parameters: 28
% kinetics:8
% arcweigts: 16
% inhibition: 4
%% best results
% none
x(x <0) = 0;
k = x(1:8);
w = x(9:24);
v = x(25:28);

%p = {'HS','pERK','CytoC','cleavedCasp3','cleavedPARP','ERK','Casp3','PARP','Apoptose','pRIP3','TNFa','RIP3','Nekrose'}; % places
% 6,7,8,12 assumed to be irrelevant

% old proteins
ydot(1) = hs(t); % HS
ydot(2) = (w(3)*k(2)) ... % pERK
    +w(2)*(k(1)*y(1)^w(1))-w(3)*(k(2)*y(2)^w(3));    
ydot(3) = (w(5)*k(3)*v(2)/(v(1)*v(2)+v(1)+v(2))-w(4)*k(2)) ... % Cytochrom C 
    + w(4)*(k(2)*y(2)^w(3)) - w(5)*k(3)*y(3)*y(3)^w(5) / (y(3)+v(1)*(1+y(11)/v(2))); 
ydot(4) = (w(7)*k(4)*v(4)/(v(4)*v(3)+v(4)+v(3)) - w(6)*k(3)*v(2)/(v(2)*v(1)+v(2)+v(1))) ... %cleaved Caspase 3
    + w(6)*k(3)*y(3)*y(3)^w(5) / (y(3)+v(1)*(1+y(11)/v(2))) - w(7)*k(4)*y(4)*y(4)^w(7) / (y(4)+v(3)*(1+y(11)/v(4))); 
ydot(5) = (w(9)*k(5)-w(8)*k(4)*v(4)/(v(4)*v(3)+v(4)+v(3))) ... % cleaved PARP
    + w(8)*k(4)*y(4)*y(4)^w(7) / (y(4)+v(3)*(1+y(11)/v(4))) - w(9)*(k(5)*y(5)^w(9)); 
ydot(6) = 0; % ERK 
ydot(7) = 0; % Casp 3
ydot(8) = 0; % PARP
ydot(9) = -(w(15)*k(5)) + w(15)*k(5)*y(5)^w(9); % Apoptosis
ydot(10) = (w(13)*k(8)-w(11)*k(7)) + w(11)*(k(7)*y(11)^w(10)) - w(13)*(k(8)*y(10)^w(13)); % pRIP3
ydot(11) = (w(10)*k(7)) + w(12)*(k(6)*y(1)^w(1)) - w(10)*(k(7)*y(11)^w(10)); % TNF-alpha
ydot(12) = 0; % RIP3
ydot(13) = -(w(16)*k(7)) + w(16)*k(7)*y(10)^w(13); % Necroptosis
ydot = ydot';
end