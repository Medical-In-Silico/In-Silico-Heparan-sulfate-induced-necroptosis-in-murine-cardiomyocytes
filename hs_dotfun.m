% % % % % % % % % % % % % % % % % % % % % % % % % % % %
% hs_dotfun(t, hs0) 
% is the derivative of a prescribed input function that
% is approximated for the heparan sulfate behavior over time.
%% parameters
% t : time since induction of hs0
% hs0 : initial amount [µg/ml] of heparan sulfate
%% further comments
% The approximated function is basically a Gaussian bell.
% Its DERIVATIVE is calculated by ydot_hs. Used for providing a HS dynamic within
% the ODE system.
% a, b, s, c, t0 are the parameters of the Gaussian bell.
% % % % % % % % % % % % % % % % % % % % % % % % % % % %
function ydot_hs = hs_dotfun(t,hs0)
a = 2.032;
b = 15.27;
s = 0.35;
c = 17.25*s;
t0 = 14 - b;

ydot_hs = -2* hs0/1.95 .* a./c.^2 .* (t+t0) .* exp(-(t+t0).^2/(c.^2));
end