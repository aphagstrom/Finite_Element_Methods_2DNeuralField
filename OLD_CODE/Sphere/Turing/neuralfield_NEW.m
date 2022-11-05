clc
clear

%Input
    %e - edges/connectivity of the nodes, ne is number of connections.
    %first column is e(first point (1), mesh k).
    %p - points/nodal cooredinates, np is number of points n the mesh.
    %p(x-coordinate (1), point k), p(y-(2), point k)
    %t - triangle cooredinates, nt is number of triangles
plt1=load('sphere750.mat')
%load('BRAIN.mat')

t = plt1.T'; 
p = plt1.X';
%t = nfv.faces'; %use for neuroimaging data BRAIN.MAT
%p = nfv.vertices'; %use for neuroimaging data BRAIN.MAT


%AREAS OF TRIANGLES
%number of element triangles
% ne=size(param.t,2);
ne=size(t,2);

%number of points
np=size(p,2);

%global matrix
K=sparse(np,np);
M=sparse(np,np);

x=p(1,:);
y=p(2,:);
z=p(3,:);

dett = @(J)sqrt(abs(det(J*J')));
L = [   1/3        1/3               1/3;       
        0.47014206 0.05971587 0.47014206;
        0.47014206 0.47014206 0.05971587;
        0.05971587 0.47014206 0.47014206;
        0.10128651 0.79742699 0.10128651;
        0.10128651 0.10128651 0.79742699;
        0.79742699 0.10128651 0.10128651;
        ];
%weights
W = [   0.225
        0.13239415
        0.13239415
        0.13239415
        0.12593918
        0.12593918
        0.12593918];

for i=1:ne
    i1=t(1,i);
    i2=t(2,i);
    i3=t(3,i);
    
    I=[i1 i2 i3];
    
    %vertex points
    
    x1=p(1,i1); %vertex 1 xcoordinate
    x2=p(1,i2); %vertex 2 xcoordinate
    x3=p(1,i3); %vertex 3 xcoordinate
    
    y1=p(2,i1); %vertex 1 ycoordinate
    y2=p(2,i2); %vertex 2 ycoordinate
    y3=p(2,i3); %vertex 3 ycoordinate
    
    z1=p(3,i1); %vertex 1 zcoordinate
    z2=p(3,i2); %vertex 2 zcoordinate
    z3=p(3,i3); %vertex 3 zcoordinate
    
    J=[x2-x1  y2-y1  z2-z1;
       x3-x1  y3-y1  z3-z1];
    
    Jinv=pinv(J);
    
    for h=1:7
    
        %(J^_(-1) *  NablaHat((N_i)))
    
        g1=Jinv*[-1;-1];
        g2=Jinv*[1;0];
        g3=Jinv*[0;1];
    
        g=[g1 g2 g3];
  

        %A Matrix

        for l=1:3
            for j=1:3
                K(I(:,l),I(:,j)) = K(I(:,l),I(:,j)) ...
                    + W(h)*abs(dett(J)) * dot(g(:,l),g(:,j));
                
                %Mass Matrix
                M(I(:,l),I(:,j)) = M(I(:,l),I(:,j)) ...
                    + W(h)*abs(dett(J)) * L(h,l) * L(h,j);
            end
        end

    end  
    
end



nt=10;%Number of time steps

%parameters that yield Turing pattern on Sphere
v=5.78758;
eta_0=8.12;
Delta=0.5;
alfa=5;
kappaV=0.262442;
kappaS=12;
tau=1;
maxiter=200;


param.p=p;
param.t=t;
param.v = v;
param.eta_0 = eta_0;
param.Delta = Delta;
param.alfa = alfa;
param.kappaV = kappaV;
param.kappaS = kappaS;
param.tau = tau;
param.M=M;
param.K=K;
param.np=np;
param.nt=nt;
param.maxiter=maxiter;


%perturbations
F0 = [0.1,0.1]';
options = optimoptions('fsolve','Display','iter'); % Option to display output
FF = fsolve(@(F)find_steady_state_NEW(F,param),F0,options); % Call solver

nterms=3;
kc=4.49022271163; %eigenvalue
phi_vec=[0,2*pi/3,4*pi/3];
c_vec=[1,1,1];
[k1,k2,init] =initial_conds(nterms,kc,phi_vec,c_vec,x,y);

%init=exp(-(20*x.^2+20*z.^2));
init=exp(-(20*x.^2+20*z.^2)) + exp(-(20*x.^2+20*y.^2));
perturb = real(0.1*init)'; 


%initial conditions

U=zeros(8*np,nt);

R = FF(1)*ones(size(perturb)) + perturb;
V = FF(2)*ones(size(perturb)) + perturb;
psi = zeros(np,1);
A1 = perturb;
A2 = perturb;
A3 = perturb;
g = perturb;
pp = perturb;

U(:,1) = [R; V; psi; A1; A2; A3; g; pp];

Dxxyy = @(x)-M\(K*x);

%[V,D] = eigs(K,M,2562); %used to find eigenvalues


param.cholM = chol(M);
U0 = [R; V; psi; A1; A2; A3; g; pp];
U = U0;
dt = 0.1;
options = odeset('OutputFcn',@odetpbar,'RelTol',1e-3,'AbsTol',1e-3);
for i=1:5000
    tspan = 0:dt:1;
    U0 = U(:,end);
    % [tout,U] = ode45(@(t,U)rhs(t,U,param,M,K),tspan,U0,options);
    [tout,Unew] = ode45(@(t,U)ff(U,param),tspan,U0,options);
    Unew = Unew';
    U = [U Unew];
    for tt=1:length(tout)
        clf
        plt = trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),Unew(np+1:2*np,tt));
        shading interp
        axis equal
         colorbar
         %caxis([-0.1 1.5])
        %caxis([0,1]) 
        drawnow
    end
end
    
%%
%MAKE VIDEO after outputting data
for i=20000:24645
   clf
   plt = trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),U(1:np,i));
   axis tight
    shading interp
   colorbar
   caxis([0 1.5])
   % title(num2str(tt));
   %view(0,90)
   F(i) = getframe(gcf) ;
  drawnow
end


 % create the video writer with 1 fps
  writerObj = VideoWriter('TuringSPHERE2.avi');
  writerObj.FrameRate = 5;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
%ax=gca();
for i=20000:length(F)
    %convert the image to a frame
    %set(ax, 'XLimMode', 'auto', 'YLimMode', 'auto');
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


%%
%% redraw
figure;
trimesh(plt1.T,x,y,z)
hold on
axis tight
% Rewatch at faster speed:
W=pi*U(1:np,:)+1i*U(np+1:2*np,:);
Z=(1-conj(W))./(1+conj(W));

for tt=1:200:600
    clf
    plt = trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),abs(Z(1:np,tt)));
    shading interp
    axis equal
    colorbar
    %title(num2str(tt));
    %view(0,90)
    drawnow
end
xlabel('x')
ylabel('y')
zlabel('z')
hold on
plot3(p(1,1),p(2,1),p(3,1),'-o','MarkerSize',10,...
'MarkerEdgeColor','black', 'MarkerFaceColor','red')






