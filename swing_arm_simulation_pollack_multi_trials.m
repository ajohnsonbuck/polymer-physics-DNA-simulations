clear all
close all

% Copyright 2014 Alex Johnson-Buck, The University of Michigan
%
% Performs polymer physics simulations of engineered DNA swinging arm constructs 
% to determine local effective concentration of functional components as described in 
% J. Fu et al. Nature Nanotechnology 9, pages 531–536 (2014).
%
% For flexible polymer simulations (ssDNA), a freely-rotating chain model with electrostatics is used,
% after Meisburger et al. (Lois Pollack) Biopolymers, 2013, “Polyelectrolyte properties of single stranded DNA measured using SAXS and single molecule FRET: Beyond the wormlike chain model
% For rigid polymer simulations (dsDNA), a rigid rod model is used

% User-specified parameters

d = 7; % Distance between anchor points of swinging arm and enzyme anchor site, in nm

% N = 23; % Number of nucleotides in chain
N = 10;
t = 1000; % Number of trials per run
nruns = 1;

la = 0.34*10; % Length of the helix holding WN1 to the DX tile, in nm
%Note: this value assumes a persistence length of 2 nm!
r_cutoff_max = la+2.88; % Value of 2.88 nm is an estimate of the rms length of the 3t+3t linkers, modeling them as a single 6t wormlike chain
r_cutoff_min = la-2.88;

totalhist = zeros(1,301);

volmap = cell(50,1);
for i = 1:size(volmap,1)
    volmap{i} = zeros(50,50);
end

rinclude = zeros(3,0);
rtotal = zeros(3,0);

for m = 1:nruns
    disp('Run:');
    disp(num2str(m));
    [r r2 rmsr] = swing_arm_simulation_pollack_finitewidth(N,t); % Sample t conformations of polymer with N monomers
    r = r'*1E9;
    r2 = r2'*1E18;
    rmsr = rmsr'*1E9;
    rtotal = cat(2,rtotal,r);   
    radj = r;
    radj(1,:) = radj(1,:)+d;
    
    r2adj = sum(radj.^2,1);
    rmsradj = sqrt(r2adj);

    for i = 1:size(radj,2)
        if rmsradj(i) < r_cutoff_max && rmsradj(i) > r_cutoff_min && radj(3,i) > 0
            include(i) = 1;
            rinclude = cat(2,rinclude,radj(:,i)+20);
        else
            include(i) = 0;
        end
    end

    f1 = sum(include)/size(radj,2);
    
    v1 = 4/3*pi*(r_cutoff_max*1E-8)^3;
    v2 = 4/3*pi*(r_cutoff_min*1E-8)^3;
    V = (v1-v2)/2;
    
    c_eff = f1/6.022E23/V; % local effective concentration
    
    
    [A B] = hist(rmsradj(1,:),0:0.1:30);
    
    
    f(m) = f1;
    C(m) = c_eff;

    for p = 1:size(radj,2)
        volmap{round(radj(2,p)+20)}(round(radj(1,p)+20),round(radj(3,p)+20)) = volmap{round(radj(2,p)+20)}(round(radj(1,p)+20),round(radj(3,p)+20)) + 1;
    end

totalhist = totalhist+A;
end

figure(1)  % Plot histogram of distances of cofactor from enzyme anchor point
bar(B,totalhist);
ylabel('counts');
xlabel('Distance from anchor (nm)');
xlim([-4 24]);
hold on
ymaxval = get(gca,'YLim');
ymaxval = ymaxval(2);
plot(r_cutoff_min*[1 1],[0 ymaxval],'r');
plot(r_cutoff_max*[1 1],[0 ymaxval],'r');
hold off

volmap2 = (volmap{18}'+volmap{19}'+volmap{20}'+volmap{21}'+volmap{22}');
volmap2 = volmap2/max(max(volmap2));

figure(2)
contourf(volmap2,1000,'LineStyle', 'none');

colormap(jet);
axis square;
xlim([5 50]);
ylim([5 50]);

volmap3 = zeros(size(volmap,1),size(volmap,1),0);
for i = 1:size(volmap,1);
    volmap3 = cat(3,volmap3,volmap{i});
end

[gx gz gy] = gradient(volmap3);

for i = 1:size(rinclude,2)
    rinclude2 = round(rinclude(:,i));
    rx = rinclude2(1);
    ry = rinclude2(2);
    rz = rinclude2(3);
    F(i) = 1.38E-23*298*(sqrt(((gx(rx,rz,ry))^2+(gy(rx,rz,rz))^2+(gz(rx,rz,ry))^2))/1E-9)/volmap3(rx,rz,ry);
end

Fmean = mean(F); % Assumes that all bound states contribute equally to the mean force felt by the duplex.

xt = 3.8E-9; % Distance to transition state, in m, estimated by worm-like chain model and assuming 4.5 base pairs have unwound in t.s.
% xt = 1E-9;

kfactor = exp(mean(F)*xt/(1.38E-23*298)); % xt is the distance to the transition state
