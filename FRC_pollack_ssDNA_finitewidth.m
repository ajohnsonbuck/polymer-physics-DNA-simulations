function [Wchain Rtot] = FRC_pollack_ssDNA_finitewidth(N,plotyn)
% Freely rotating chain polymer model with electrostatic self-repulsion,
% based on Meisburger et al. Biopolymers. 2013 Dec; 99(12)
% Inputs:
% N is the number of nucleotides in the chain.
% plotyn is a string specifying whether to plot the conformation of each
% chain generated.  'y' will result in plotting; any other argument will
% yield no plotting.

l = 3.69; % Virtual bond length, in Angstroms
% l = 1;
theta = 57.8; % Valence angle, in degrees
theta = theta*pi/180;
% theta = pi/100;
f = 0.229; % Charge renormalization factor, determined empirically
K = 1/13.4; % in A-1
a = 5.6; % in A
d = a; % in A

kb = 1.38065E-23; % in J/K
lb = 7.14; % Bjerrum length, in Angstroms
T = 293.15; % in K


%Convert all length values to meters;
l = l*10^-10;
K = K*10^10;
a = a*10^-10;
d = a;
lb = lb*10^-10;

na = 2*N; % na is the number of virtual atoms 
nb = na-1; % nb is the number of virtual bonds 
nv = na-1; % the number of valence angles
nd = na-1; % the number of dihedral angles

Z = zeros(na,1);
for j = 1:na
if mod(j,2)==0
    Z(j,1) = -1;
end
end

Z2 = Z*Z';

for j = 1:na % Remove self-interactions (diagonal elements in Z^2)
    k = j;
    Z2(j,k) = 0;
end

goodchain = 0;

tries = 1;

while goodchain == 0

dihed = random('unif',0,2*pi,nd,1);

R = zeros(na,3);

R(1,:) = [0 0 l];
R(2,:) = [l*sin(theta)*cos(dihed(1)) l*sin(theta)*sin(dihed(1)) l*cos(theta)];

for n = 3:na
    rotaxis = cross([0 0 1],R(n-1,:));
    rotaxis = rotaxis/norm(rotaxis);
    u = rotaxis(1);
    v = rotaxis(2);
    w = rotaxis(3);
    th = acos(dot(R(n-1,:),[0 0 1])/l);
    rotmat = [u^2+(1-u^2)*cos(th) u*v*(1-cos(th))-w*sin(th) u*w*(1-cos(th))+v*sin(th); 
              u*v*(1-cos(th))+w*sin(th) v^2+(1-v^2)*cos(th) v*w*(1-cos(th))-u*sin(th);
              u*w*(1-cos(th))-v*sin(th) v*w*(1-cos(th))+u*sin(th) w^2+(1-w^2)*cos(th)];
    R(n,:) = (rotmat*[l*sin(theta)*cos(dihed(n-1)) l*sin(theta)*sin(dihed(n-1)) l*cos(theta)]')';
end

Rtot = NaN(size(R,1),3);

for i = 1:size(R,1)
Rtot(i,:) = sum(R(1:i,:),1);
end

% if min(Rtot(:,3)) < 0   % Prevent chain from crossing 2-helix tile

goodchain = 1;

for i = 1:size(Rtot,1)
    if goodchain == 1
        if Rtot(i,3) < 0
            if Rtot(i,3) > -2E-9
                if abs(Rtot(i,2)) < 2.5E-9
                    goodchain = 0;
                    %     disp('Chain crosses z = 0 plane');
                    tries = tries+1;
                end
            end
        end
    end
end

if goodchain == 1
    D = pdist(Rtot);
    D = squareform(D);
    for j = 1:size(D,1)
        for k = 1:size(D,2)
            if j == k || abs(j-k) < 2
                D(j,k) = 50;
            end
        end
    end
    if min(min(D)) < a
        goodchain = 0;
%         disp('Chain collides with itself');
        tries = tries+1;
    else
        goodchain = 1;
    end
end
end


if ischar(plotyn) == 1
if strcmp(plotyn,'y') == 1
plot3(Rtot(:,1),Rtot(:,2),Rtot(:,3),'o-','MarkerFaceColor','r','MarkerEdgeColor','k','Color','k');
axis('vis3d');
% axis('square');
axis('equal');
grid;
xlabel('X (m)','FontSize',16);
ylabel('Y (m)','FontSize',16);
zlabel('Z (m)','FontSize',16);
set(gca,'LineWidth',3,'FontSize',16);
end
end

% disp(num2str(tries));
% max(Rtot,[],1)

W = zeros(na,na);

for j = 1:na
    for k = 1:na
        W(j,k) = kb*T*lb*(Z2(j,k))/(1-K*d/2)^2*exp(-K*(D(j,k)-d))/D(j,k);
    end
end

Wchain = f^2*0.5*sum(sum(W));

end

