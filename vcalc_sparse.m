clear all 
%K=exp(5*rand(200));  %debug
K=load('kout_exp_1','-ascii');
K=log(K);
K=K-mean(mean(K))+1;   % adjust mean(log(K)), although it makes no difference in hillslope case.
K=exp(K);
nn=size(K); t1=cputime;  L=100; H=25;
dx=L/nn(2); dy=H/nn(1);
% iteration parameters: 1<om<2 acceleration, pointwise tolerance, max iterations
hriv=0; I=2.5e-4;  por=0.25;
Infilt=I/dy;
% make an intial head map guess that goes uniformly from 1 to 0, L to R.
%h=repmat(nn(1)-1:-1:0,nn(1),1)*(1/(nn-1))*(hhi-hlo)+hlo;
h=zeros(size(K));
for k=1:nn-2
    x0(k*(nn-2)-(nn-2)+1:k*(nn-2))=h(1,k+1);
end
x0=x0';  % this is a vector of the uniform guess heads in the eventual Ax=b
figure(10), imagesc(log(K)); axis tight; axis equal;
hold on;
%K=exp(K);   % exponentiate OS-logNormal K values
%get K harmonic means between nodes
Kiave=ones(nn(1),nn(2)-1);
Kjave=ones(nn(1)-1,nn(2));

for  i=1:nn(1)-1;for j=1:nn(2);
  Kiave(i,j)=2/(1/K(i,j)+1/K(i+1,j));
    end;end;
for  i=1:nn(1);for j=1:nn(2)-1;
  Kjave(i,j)=2/(1/K(i,j)+1/K(i,j+1));
    end;end;

% construct sparse coefficient matrix by storing only diagonals
A=zeros(nn(1)*nn(2)+nn(1),5);
%where are the diags?:
dis=[-nn(1) -1 0 1 nn(1)];
nntot=nn(1)*nn(2);
k=1
for j=1:nn(2); for i=1:nn(1);
    if i==1
        Kup=0; 
    else
        Kup=Kiave(i-1,j)/dy^2;
    end
    
    if i==nn(1) 
        Kdown=0; 
    else
        Kdown=Kiave(i,j)/dy^2;
    end
    
    if j==1
        Kleft=0;
    else
        Kleft=Kjave(i,j-1)/dx^2;
    end
        
    if j==nn(2)  
        Kright=K(i,j)/dx^2;  % for constant head at river
    else
        Kright=Kjave(i,j)/dx^2;
    end
        
    A(k,3)=-(Kup+Kdown+Kright+Kleft);
    if k>1 A(k-1,2)=Kup; end
    A(k+1,4)=Kdown;
    if k>nn(1) A(k-nn(1),1)=Kleft; end
    A(k+nn(1),5)=Kright;
    
    k=k+1;
 end; end; 

% add the no-flow boundary at top (i=1) and bottom (i=nn) to main 3 diag terms
%build the known vector.
b=zeros(nn(1)*nn(2),1);

b(nntot-nn(1)+1:nntot)=hriv*(1/dx^2)*Kjave(1:nn(1),nn(2));  %river on Right
for j=1:nn(2)
  b(nn(1)*(j-1)+1)=b(nn(1)*(j-1)+1)-I/dy;
end

%A=A;b=b;
P=spdiags(A,dis,nntot,nntot);  %assemble sparse matrix P from diags in A

%clear A; 
x=inv(P)*b;
% Solve implicit with iteration.  These use L and U for conditioning.
% Any of these three solvers seem to work about the same.  Uncomment 
% the following line and one of the next three.  VERY stable solutions
% but slow. And two big pre-conditioning matrices.
%[L,U] = ilu(P,struct('type','ilutp','droptol',1e-6));
%x=cgs(P,b,1e-12,200,L,U,x0);
%x=gmres(P,b,10,1e-12,200,L,U,x0);
%x=bicgstab(P,b,1e-12,200,L,U,x0);

% This seems to work quite well, using inc. Cholesky conditioner and pcg.
%M = ichol(P, struct('type','ict','droptol',1e-3));
%x=bicgstab(P,b,1e-12,200,M,M',x0);
%x=pcg(P,b,1e-12,200,M,M',x0);
t1=cputime-t1

blah=reshape(x,nn(1),nn(2));   % put solved head vector x back into 2-d map
h=blah;

contour(h,40,'black');   %place head contours on top of ln(K) map
v=3*ones(nn(1),nn(2),2);vmag=ones(nn-1);
%calculate velocity by darcy's law between nodes (normalize by distance 1/nn)
for i=1:nn(1)-1;for j=1:nn(2)-1;
    v(i,j,2)=(1/dy)*Kiave(i,j)*(h(i+1,j)-h(i,j));    
    v(i,j,1)=(-1/dx)*Kjave(i,j)*(h(i,j+1)-h(i,j));  %(X- or j-velocity)
%    vmag(i,j)=sqrt(v(i,j,1)^2+v(i,j,2)^2);
    end;end;
%  Add BC stuff
v(1:nn(1),nn(2),1)=(-1/dx)*K(:,nn(2)).*(hriv-h(:,nn(2)));  % x-vel on left bdy
v(nn(1),1:nn(2)-1,1)=(-1/dx)*Kjave(nn(1),1:nn(2)-1).*(h(nn(1),2:nn(2))-h(nn(1),1:nn(2)-1));
v(nn(1),:,2)=0;
v(1:nn(1)-1,nn(2),2)=(1/dy)*Kiave(1:nn(1)-1,nn(2)).*(h(2:nn(1),nn(2))-h(1:nn(1)-1,nn(2)-1));
v=v/por;  vmag=sqrt(v(:,:,1).^2+v(:,:,2).^2);
figure(3), imagesc(log(vmag))  % color map of velocity magnitude.
hold on
axis equal; axis tight
contour(h,40,'black');   %place head contours on top 
hold off

% comment this out typically - it saves K, h, and v subfields for other
% programs
%xmax=800;ymax=400;blah1=K(1:xmax,1:ymax);blah2=h(1:xmax,1:ymax),blah3=v(1:xmax-1,1:ymax-1,:);
%save('Velocities_fractal.mat','blah1','blah2','blah3');
save('Velocities_exp.mat','K','h','v','L','H','I','por');

