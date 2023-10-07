function [rmsd Rg] = fjc_rmsd(N,l,Lp)
%   Calculate rms end-to-end distance <rmsd> and radius of gyration <Rg>
%   For a freely-jointed chain polymer with <N> monomers of per-monomer contour
%   length <l> and persistence length <Lp> (= 1/2 Kuhn length)
%   N = number of monomers in chain
%   l = contour length per monomer
%   Lp = persistence length (NOT Kuhn length)

L = 2*Lp; % Kuhn length
Lc = N*l; % Contour length of polymer

rmsd = sqrt(L*Lc);
Rg = rmsd/sqrt(6);

end

