clc
clear

dx = (pi-(-pi))./50;   %0.05;
x = -pi:dx:pi ;      %-1:dx:1;
y =  -pi:dx:pi  ;    %-1:dx:1;
[X,Y] = meshgrid(x,y);
DT = DelaunayTri(X(:),Y(:));
p = DT.X';
t = DT.Triangulation';

%AREAS OF TRIANGLES
%number of element triangles
ne=size(t,2);

%number of points
np=size(p,2);

p=[p; zeros(1,np)];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify M and K to incorporate periodic boundary conditions
% This code will only work for the mesh obtained from DelaunayTri above.
tol = 1e-13;
ileft = find(abs(x+1)<tol);
iright = find(abs(x-1)<tol);
ibottom = find(abs(y+1)<tol);
itop = find(abs(y-1)<tol);
K0 = K;
M0 = M;
for i = 1:length(ileft)
    K(ileft(i),:)  = K(ileft(i),:)  + K0(iright(i),:);
    K(iright(i),:) = K(iright(i),:) + K0(ileft(i),:);
    
    K(ibottom(i),:) = K(ibottom(i),:) + K0(itop(i),:);
    K(itop(i),:)    = K(itop(i),:)    + K0(ibottom(i),:);
    
    M(ileft(i),:)  = M(ileft(i),:)  + M0(iright(i),:);
    M(iright(i),:) = M(iright(i),:) + M0(ileft(i),:);
    
    M(ibottom(i),:) = M(ibottom(i),:) + M0(itop(i),:);
    M(itop(i),:)    = M(itop(i),:)    + M0(ibottom(i),:);
    
end

itl = intersect(itop,ileft);
itr = intersect(itop,iright);
ibl = intersect(ibottom,ileft);
ibr = intersect(ibottom,iright);

K(itl,:) = K(itl,:) + K0(ibr,:);
K(itr,:) = K(itr,:) + K0(ibl,:);
K(ibr,:) = K(ibr,:) + K0(itl,:);
K(ibl,:) = K(ibl,:) + K0(itr,:);

M(itl,:) = M(itl,:) + M0(ibr,:);
M(itr,:) = M(itr,:) + M0(ibl,:);
M(ibr,:) = M(ibr,:) + M0(itl,:);
M(ibl,:) = M(ibl,:) + M0(itr,:);

B = zeros(np,length(union(iright,itop)));
for i=1:length(ileft)
    B([ileft(i) iright(i)],i) = [1 -1];
    if i<length(ileft)
        B([ibottom(i) itop(i)],length(ileft)+i) = [1 -1];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%--- Time parameters ---%%%
dt=0.001;%0.1;          % Size of time step
nt=500;%10;          %Number of time steps

eta_0=0.4;
Delta=0.5;
v=0.336756;
kappaV=1.3207;
kappaS=12;
tau=1;


%v = 0.172923
%eta_0 = 0.4
%Delta = 0.5
% = 5.0
%kappaV = 1.28773
%kappaS = 12.0
%tau = 1
%tau_a = 0.5
%delta = 0.0


maxiter=200;

param.B = B;
param.v = v;
param.eta_0 = eta_0;
param.Delta = Delta;
param.alfa = alfa;
param.kappaV = kappaV;
param.kappaS = kappaS;
param.tau = tau;
param.delta = delta;
param.M=M;
param.K=K;
param.np=np;
param.nt=nt;
param.maxiter=maxiter;
%perturbations
F0 = [1.0; 1.2]';
options = optimoptions('fsolve','Display','iter'); % Option to display output
FF = fsolve(@(F)find_steady_state_NEW(F,param),F0,options); % Call solver


nterms=3;
% kc=1;
% kc=10;
kc = 6*pi;
kc2 = 2;
phi_vec=[0,2*pi/3,4*pi/3];
c_vec=[1,1,1];
[k1,k2,init] =initial_conds(nterms,kc,phi_vec,c_vec,x,y);
[k1,k2,init2] =initial_conds(nterms,kc2,phi_vec,c_vec,x,y);

% init=exp(-(20*x.^2+20*z.^2));
%init=exp(-(20*x.^2+20*z.^2)) + exp(-(20*x.^2+20*y.^2));

% perturb = real(0.001*init)';
% perturb = real(0.1*init)';
perturb = real(0.01*(init2))';

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

% Dxxyy = @(x)-M\(K*x);

U0 = [R; V; psi; A1; A2; A3; g; pp];
U = U0;
dt = 0.1;
options = odeset('OutputFcn',@odetpbar,'RelTol',1e-3,'AbsTol',1e-3);
for i=1:5000
    tspan = (0:dt:1);
    U0 = U(:,end);
    %[tout,Unew] = ode45(@(t,U)ff(U,param),tspan,U0,options);
    [tout,Unew] = ode45(@(t,U)ff_periodic(U,param),tspan,U0,options);
    Unew = Unew';
    U = [U Unew];
    for tt=1:length(tout)
        clf
        plt = trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),Unew(1:np,tt));
        shading interp
        %title([num2str(tt),'  ',num2str(max(abs(Unew(:,tt))))]);
        colorbar
        axis tight
         %caxis([-.006 0.006])
        view(0,90)
        %F(i) = getframe(gcf) 
        drawnow
    end
end
    
%%
% create the video writer with 1 fps
  writerObj = VideoWriter('TuringHopf.avi');
  writerObj.FrameRate = 1;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
%ax=gca();
for i=1:length(F)
    %convert the image to a frame
    %set(ax, 'XLimMode', 'auto', 'YLimMode', 'auto');
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);



%% redraw
figure;
trimesh(DT,x,y)
hold on
axis tight

% Rewatch at faster speed:
for tt=1:200:500
    clf
    plt = trisurf(t(1:3,:)',p(1,:),p(2,:),p(3,:),U(1:np,tt));
    %shading interp
    axis equal
    %title(num2str(tt));
    view(0,90)
    drawnow
end
xlabel('x')
ylabel('y')
hold on
plot(p(1,1677),p(2,1677),'-o','MarkerSize',5,...
'MarkerEdgeColor','black', 'MarkerFaceColor','red')
