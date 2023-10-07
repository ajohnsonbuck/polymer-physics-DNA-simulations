function [r r2 rmsd] = swing_arm_simulation_pollack_finitewidth(N, t)
% Sample multiple polymer conformations using freely rotating chain model with electrostatic self-repulsion
% N = number of nucleotides
% t = 1000; % number of conformations to sample

kb = 1.38065E-23;
T = 293.15;
Rtot = cell(2,1);

[Wchain(1) Rtot{1}] = FRC_pollack_ssDNA_finitewidth(N,'n'); % Initial conformation

num = random('unif',0,1,t,1); % Generate list of pseudo-random numbers for Monte Carlo sampling
r = NaN(t,3);

na = size(Rtot{1},1);

for m = 1:t
    if mod(m,100) == 0
        disp(num2str(m));
    end
    [Wchain(2) Rtot{2}] = FRC_pollack_ssDNA_finitewidth(N,'n'); % Candidate next conformation
    ratio = exp(-Wchain(2)/(kb*T))/exp(-Wchain(1)/(kb*T)); % Compare free energy of candidate conformation to present conformation
    if ratio < 1  % Monte Carlo sampling: Decide whether to transition to candidate conformation based on energetics and pseudo-random numbers
        if num(m)<ratio
            Wchain(1) = Wchain(2);
            Rtot{1} = Rtot{2};
            r(m,:) = Rtot{1}(na,:);
        else
            r(m,:) = Rtot{1}(na,:);
        end
    else
        Wchain(1) = Wchain(2);
        Rtot{1} = Rtot{2};
        r(m,:) = Rtot{1}(na,:);
    end
end

r2 = r.^2;
d = sum(r2,2).^0.5;
rmsd = mean(d);

end

