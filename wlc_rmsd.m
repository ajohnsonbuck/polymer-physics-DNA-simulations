function [rmsd] = wlc_rmsd(N,unitlength,p)
% Calculate rmsd end-to-end distance of a worm-like chain polymer with <N> monomers 
% of length <unitlength> and persistence length <p>
% N = number of monomers
% unitlength = contour length per unit of polymer
% p = persistence length

l = unitlength;
L = l*N;

rmsd = sqrt(2*p*L*(1-p/(L)*(1-exp(-L/p))));

end

