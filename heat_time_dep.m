% %a(u,v)= int(gradu * gradv)dx + int(vf)dx
%where f=1
%A=int(gradu * gradv)dx=Sum(phi_i *phi_j)*a
%b=int(vf)dx


dbclear if error

%create geometric model
model = createpde(1);

g = decsg(gd,sf,ns);
geometryFromEdges(model,g);
generateMesh(model,'GeometricOrder','quadratic');
[p,e,t] = meshToPet(model.Mesh);



%Input
    %e - edges/connectivity of the nodes, ne is number of connections.
    %first column is e(first point (1), mesh k).
    %p - points/nodal cooredinates, np is number of points n the mesh.
    %p(x-coordinate (1), point k), p(y-coordinate(2), point k)
    %t - triangle cooredinates, nt is number of triangles
    
    ib =unique([e(1,:),e(2,:)]);
 
%AREAS OF TRIANGLES
%number of element triangles
ne=size(t,2);

%number of points
np=size(p,2);


%global matrix
K=zeros(np);
M=zeros(np);
F=zeros(np,1);


      
   
for i=1:ne       
    i1=t(1,i); 
    i2=t(2,i);
    i3=t(3,i);
    i4=t(4,i); 
    i5=t(5,i);
    i6=t(6,i);
    
    I=[i1 i2 i3 i4 i5 i6];
    %%
   %vertex points
   
    x1=p(1,i1); %vertex 1 xcoordinate
    x2=p(1,i2); %vertex 2 xcoordinate
    x3=p(1,i3); %vertex 3 xcoordinate
  
   
    y1=p(2,i1); %vertex 1 ycoordinate
    y2=p(2,i2); %vertex 2 ycoordinate
    y3=p(2,i3); %vertex 3 ycoordinate
  
%%

    J=[x1-x3  y1-y3;  
    x2-x3 y2-y3];         
    Jinv=inv(J);
    
%Gauss Points

 for h=1:7
     
   
     L=[
         
        1/3        1/3                1/3;
        
        0.47014206 0.05971587  0.47014206;
      
        0.47014206 0.47014206  0.05971587;
      
        0.05971587 0.47014206 0.47014206;
        
        0.10128651  0.79742699 0.10128651;
        
        0.10128651 0.10128651 0.79742699;
        
        0.79742699 0.10128651 0.10128651;
      ];
 
  L1=L(h,1);
  L2=L(h,2);
  L3=L(h,3);
 
  
   N=[L1*(2.*L1-ones(1,1));
      L2*(2.*L2-ones(1,1));
      L3*(2.*L3-ones(1,1));
      4.*L2*L3;
      4.*L1*L3;
      4.*L1*L2];

  
  %%
    
    g1(:,h)=J\[4.*L1-ones(1,1);zeros(1,1)];
    g2(:,h)=J\[zeros(1,1);4*L2-ones(1,1)];
    g3(:,h)=J\[-4.*L3+ones(1,1);-4.*L3+ones(1,1)];
    g4(:,h)=J\[-4.*L2;4.*L3-4.*L2];
    g5(:,h)=J\[4.*L3-4.*L1;-4.*L1];
    g6(:,h)=J\[4.*L2;4.*L1];

    g=[g1(:,h) g2(:,h) g3(:,h) g4(:,h) g5(:,h) g6(:,h)];
 
    %K Matrix
    
    %weights 
    W=[0.225        0.225         0.225
        
      0.13239415    0.13239415    0.13239415 
      0.13239415    0.13239415    0.13239415 
      0.13239415    0.13239415    0.13239415 
      
      0.12593918    0.12593918    0.12593918
      0.12593918    0.12593918    0.12593918
      0.12593918    0.12593918    0.12593918
      ];
    
     for i=1:6
       for j=1:6
           for k=1:1
    K(I(:,i),I(:,j))=K(I(:,i),I(:,j))...
    +W(h,k)*det(J) * dot(g(:,i),g(:,j));
           end
       end
     end
  
   
     %M Matrix
   for i=1:6
       for j=1:6
           for k=1:1
    M(I(:,i),I(:,j))=M(I(:,i),I(:,j))...
    +W(h,k)*det(J) * dot(N(i),N(j));
           end
       end   
   end
   

 end
 
   
       
      
    %F matrix
  
    F(I(:,1))=F(i1)+det(J)/3;
    F(I(:,2))=F(i2)+det(J)/3;
    F(I(:,3))=F(i3)+det(J)/3;
    F(I(:,4))=F(i4)+det(J)/3;
    F(I(:,5))=F(i5)+det(J)/3;
    F(I(:,6))=F(i6)+det(J)/3;

 end


%%%--- Initial conditions ---%%%


%%%--- Time parameters ---%%%
dt=1;          % Size of time step
nt=100;
    
r=zeros(2,np)+25;
u=zeros(2,np);
for i=1:1133
   if p(1,i).^2+p(2,i).^2<r(1,i).^2
      u(1,i)=r(1,i).^2-p(1,i).^2-p(2,i).^2;
   else 
      u(1,i)=zeros(1,1); 
   end
end

u(2,:)=[];



F1=zeros(np,nt);
F1(:,1)=ones(1133,1);

T1=zeros(np,nt);
T1(:,1)=u';

for tt=1:nt-1

A= M+1/2.*dt.*(K);

B=(M-1/2.*dt.*(K))*T1(:,tt)...
    
   +dt.*(1/2.*F1(:,tt+1)+1/2.*F1(:,tt));

  
   
interior= ones(np,1);
interior(ib)=0;
interior=logical(interior);

A=A(interior,interior);
B=B(interior);

T1(interior,tt+1)=inv(A)*B;
end

for i=1:nt
trisurf(t(1:3,:)',p(1,:),p(2,:),T1(:,i));

axis([-25 25 -25 25 0 1000]);
drawnow


pause(0.02)
end







