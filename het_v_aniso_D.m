clear all
% vertrecharge is vertical velocity BEFORE hitting aquifer, 
% x-length, height: lengths in m, time in days  
Dmol=2.5e-3; atrans=0.0; along=0.0;  % atrans is isotropic dispersivity (m), 
% along is extra longitudinal
% For same flux, because v dot n = v_z(ztop) is the flux, and v_z(ztop) = -I*ztop
%I=vertrecharge/H;
eps=1e-8;
load('Velocities_exp.mat')  %  'K','h','v','L','H','I','por'
nn=size(K); dx=L/nn(2); dz=H/nn(1);
vz=zeros(nn(1)+1,nn(2)+1); vx=zeros(nn(1)+1,nn(2)+1);
vx(2:nn(1)+1,2:nn(2)+1) = v(:,:,1); 
vz(2:nn(1)+1,2:nn(2)+1) = v(:,:,2); 
vz(1,:)=-I;

meantime=H*por/I;
maxtime=10*meantime;
vmax=(I/H/por)*sqrt(L^2+H^2); Dmax=Dmol+atrans*vmax;  % For D=D(v)
Peclet=I*H/Dmol/por
dt=.5*por/I;
maxdtadv=max( dx/max(max(abs(vx))), dz/max(max(abs(vz))));
nadvect=1+ceil(dt/maxdtadv)    % How many sub-steps for advection
nadvect=max([5 nadvect]);

ntsteps=floor(maxtime/dt); maxtime=dt*ntsteps
stats=zeros(ntsteps,4);
Nage=floor(maxtime/dt);
N=100;             % number of particles input per timestep
ds_x=L/N;  ds_z = I*dt/por;  % Horiz and vert separations
ds=ds_x*ds_z;            % sep=L/NUP (Mike - might be wrong for this 2-D prob? 
X=zeros(floor(N*maxtime/dt/20),Nage+7);     % Reserve X(:,1:7) for x,y,z,D,streamtube #,birth time,death time
dead=zeros(ntsteps*N,Nage+7); ensmeans=zeros(ntsteps,1);
size(X)
gridpeclet=L^2*H/(N*H*Dmol)
% ----------------------
Ntot=1; deadtot=0;
a=2; %factor multiplying diffusion distance - distance bigger than this do not react
currtime=0;
tic
for k=1:ntsteps
% Age existing particles
    X(1:Ntot,9:Nage+7)=X(1:Ntot,8:Nage+6);
    X(1:Ntot,8)=0;
    if k<2 
        Ntot=0;
    end
% Insert particles at upper bndry (z=0)
%    X(Ntot+1:Ntot+N,1)=L*rand(N,1);       % Random particle entry points May not be divergence free!
    X(Ntot+1:Ntot+N,1)=(L/(N+1))*[1:N]';  % Evenly spaced points (debugging?)
    X(Ntot+1:Ntot+N,3)=H-eps;     % Could randomize over depth for 0<time<dt
    X(Ntot+1:Ntot+N,6)=currtime; % Birthdate: Also could make correspond to above random start depth
    X(Ntot+1:Ntot+N,8)=1/dt;     % Ditto
    X(Ntot+1:Ntot+N,5)=[1:N]';   % last minute hack to add "streamtube number"
    Ntot=Ntot+N;
    
%  Figure out particle velocities and dispersion coeff.
     cellvj=1+floor(1+(1/dx)*X(1:Ntot,1));
     cellvi=1+floor(1+(1/dz)*(H-X(1:Ntot,3))); 
     lindexA=(nn(1)+1)*(-1+cellvj)+cellvi;  % gets linear index of each particle's cell
     cellposx=(1/dx)*X(1:Ntot,1)-cellvj+2;
     cellposz=(1/dz)*(H-X(1:Ntot,3))-cellvi+2;
     vxparts=cellposx.*vx(lindexA)+(1-cellposx).*vx(-(nn(1)+1)+lindexA);
     vzparts=cellposz.*vz(lindexA)+(1-cellposz).*vz(-1+lindexA);
     X(1:Ntot,4)=Dmol+atrans*sqrt(vxparts.^2+vzparts.^2);
     Dmax=max(X(1:Ntot,4));
     
% Find nearby particles to mix with
if(Dmax>0)
    dist=a*sqrt(2*(2*Dmax*dt));
    [idx, r]=rangesearch(X(1:Ntot,1:3),X(1:Ntot,1:3),dist,'BucketSize',50);
    
    Nclose = sum(cellfun('length', idx));
    row = zeros(1, Nclose);
    col = zeros(1, Nclose);
    val = zeros(1, Nclose);

    start = 1;

    for ii = 1 : Ntot
        finish=start-1+length(idx{ii});
%        row(start:finish)=ii*ones(1,length(idx{ii}));
        row(start:finish)=ii;
        col(start:finish)=idx{ii};
        val(start:finish)=r{ii};
        D(start:finish)=0.5*X(ii,4)+0.5*X(col(start:finish),4);  % Average D between ii particle and others
        start=finish+1;
    end

    clear idx r
    
%     calculate this before sparsing to avoid exponentiating the zero entries

%     D is formed such that the entries of val that will end up in the i^th
%     column of Sdist will correspond to the average diffusion coefficient for the 
%     i^th and all jth particles

    val =  (1 ./ (8 * pi * D * dt)) .* exp(-(val.^2) ./ (8 * D * dt));  
    
    Pmat = sparse(row, col, val);

    clear row col val
    
%     normalize the columns of Pmat to ensure they sum to 1
    Pmat = Pmat * spdiags(1./(sum(Pmat)'), 0, Ntot, Ntot);

%     this is the explicit matrix from the AWR Accuracy paper
    Emat = (speye(Ntot) - (1 / 2) * spdiags(Pmat * ones(Ntot, 1), 0, Ntot, Ntot)) + (1 / 2) * Pmat;
    
    clear Pmat
    
%     conduct the mass-transfer
%     temp = Emat * X(1 : Ntot, 6 : Nage + 5);
    X(1 : Ntot, 8 : Nage + 7) = Emat * X(1 : Ntot, 8 : Nage + 7);
    
    clear Emat
      
%  RW particles (along-atrans)|v| in the direction of transport
  if(along>atrans)
%     vx=(I/H)*X(1:Ntot,1);
%     vz=(-I/H)*X(1:Ntot,3);
     vmag=sqrt(vx.^2+vz.^2);
%     tantheta=vx./vz; 
     cospart=vx./vmag; sinpart=vz./vmag; 
     X(1:Ntot,1)=X(1:Ntot,1)+sqrt(2*(along-atrans).*vmag*dt).*cospart.*randn(Ntot,1);
%   X(1:Ntot,2)=X(:,2)+sqrt(2*D*dt)*randn(Ntot,1);
     X(1:Ntot,3)=X(1:Ntot,3)+sqrt(2*(along-atrans).*vmag*dt).*sinpart.*randn(Ntot,1);
     apple=find(X(1:Ntot,1)>L);    %  Reflect particles gone too far
     X(apple,1)=2*L-X(apple,1);
     apple=find(X(1:Ntot,3)<0);
     X(apple,3)=-X(apple,3);
  end
end
%  advect particles  May need to refine this step using forward Euler, or
%  do RK4
   numdeaddt=0;
 %  nadvect=5;   % This many sub-steps to get better 
                % advection & death resolution if needed


                
     for jadv=1:nadvect     
     currtime=currtime+dt/nadvect;

     cellvj=1+floor(1+(1/dx)*X(1:Ntot,1));
     cellvi=1+floor(1+(1/dz)*(H-X(1:Ntot,3))); 
%     vx=(I/H)*X(:,1);
%     vz=-(I/H)*X(:,3);
%     lindexA=zeros(NA,1);
     lindexA=(nn(1)+1)*(-1+cellvj)+cellvi;  % gets linear index of each particle's cell
     cellposx=(1/dx)*X(1:Ntot,1)-cellvj+2;
     cellposz=(1/dz)*(H-X(1:Ntot,3))-cellvi+2;
     vxparts=cellposx.*vx(lindexA)+(1-cellposx).*vx(-(nn(1)+1)+lindexA);
     vzparts=cellposz.*vz(lindexA)+(1-cellposz).*vz(-1+lindexA);
     X(1:Ntot,1)=X(1:Ntot,1)+(dt/nadvect)*vxparts;
     X(1:Ntot,3)=X(1:Ntot,3)+(dt/nadvect)*vzparts;
     X(1:Ntot,4)=Dmol+atrans*sqrt(vxparts.^2+vzparts.^2);

%     X(:,1)=X(:,1)+(dt/nadvect)*vx;
%     X(:,2)=X(:,2)+(dt/nadvect)*vy;
%     X(:,3)=X(:,3)+(dt/nadvect)*vz;

%  Convert particles to death, keeping info
     deadapple=find(X(:,1)>L);
     X(deadapple,7)=currtime;
     numdead=length(deadapple);
     if(numdead>=1)
         dead(deadtot+1:deadtot+numdead,:)=X(deadapple,:);
         X(deadapple,:)=[];
         Ntot=Ntot-numdead;
         deadtot=deadtot+numdead;
         numdeaddt=numdeaddt+numdead;
     end
   end
   
%plot age pdfs, including total mixed
   plottime=dt*[1:Nage]';
   pdfsnow=dead(deadtot-numdeaddt+1:deadtot,8:end)';
   pdfsnow(pdfsnow<1e-20)=1e-20;
   avepdf=(1/numdeaddt)*sum(pdfsnow,2);
if(numdeaddt>0)
   figure(2)
   plot(plottime,pdfsnow);
   axis([0 maxtime 0 1/dt]);
   hold on
   plot(plottime,avepdf,'-+k','LineWidth',2);
   plot(plottime,(I/por)*exp(-plottime*I/por),'-r','LineWidth',1)
   xlabel('Age (d)');ylabel('Age pdf');
   drawnow
   hold off
%pause
   meanage=sum(plottime.*avepdf)*dt;
   varage=sum((plottime-meanage).^2.*avepdf)*dt;

   figure(3)
   semilogy(plottime,pdfsnow);
   axis([0 maxtime 1e-9 2/dt]);
   hold on
   loglog(plottime,avepdf,'-+k','LineWidth',2);
   semilogy(plottime,(I/H/por)*exp(-plottime*I/H/por),'-r','LineWidth',1)
   xlabel('Age (d)');ylabel('Age pdf');
   title(['Mean age of ensemble exiting particles = ',num2str(sum(plottime.*avepdf)*dt)]); 
   drawnow
   hold off

   figure(4)
   semilogy(plottime,pdfsnow);
   axis([0 maxtime*I/H/por 1e-5 10]);
   hold on
   loglog(plottime*I/H/por,(H*por/I)*avepdf,'-+k','LineWidth',2);
   semilogy(plottime*I/H/por,(H*por/I)*(I/H/por)*exp(-plottime*I/H/por),'-r','LineWidth',1)
   xlabel('Dimensionless Age (tI/H\phi)');ylabel('Age pdf');
   title(['Mean dimensionless age of ensemble exiting particles = ',num2str(sum(plottime.*avepdf)*dt*I/H/por)]); 
   drawnow
   hold off

   
   randsample=randperm(numdeaddt,ceil(numdeaddt/10));
   randpdfs=pdfsnow(:,randsample);
   avepdf=(1/ceil(numdeaddt/10))*sum(randpdfs,2);

   figure(5)
   semilogy(plottime,randpdfs);
   axis([0 maxtime*I/H/por 1e-5 10]);
   hold on
   loglog(plottime*I/H/por,(H*por/I)*avepdf,'-+k','LineWidth',2);
   semilogy(plottime*I/H/por,(H*por/I)*(I/H/por)*exp(-plottime*I/H/por),'-r','LineWidth',1)
   xlabel('Dimensionless Age (tI/H)');ylabel('Age pdf');
   title(['Mean dimensionless age of ensemble exiting particles = ',num2str(sum(plottime.*avepdf)*dt*I/H/por)]); 
   drawnow
   hold off
  
   currtime
end   
   

   ensmeans(k,1)=sum(plottime.*avepdf)*dt;
   minsize=min(X(1:Ntot,2)); 
   diffdiam=2*sqrt(2*X(1:Ntot,4)*dt);
   minsize=.000001; 
   %diffdiam=2*sqrt(2*X(1:Ntot,2).*(currtime-X(1:Ntot,4)))
   %diffdiam=1e-8;
   
   figure(1)
  ptplot=scatter(X(1:Ntot,1),X(1:Ntot,3),diffdiam);
  daspect([1 1 1]);
  axis([0 L 0 H]);
  title(['Diameter is approximate diffusion distance per dt=',num2str(dt)]);
  xlabel('X (m)');ylabel('Z (m)');
    ax=gca;
    currentunits=ax.Units;
    ax.Position;
    ax.Units='Points';
    axpos=ax.Position;
    pba=ax.PlotBoxAspectRatio;
    ax.Units=currentunits;
    scalemarker=min([axpos(3),pba(1)/pba(2)*axpos(4)])/diff(xlim);
    ultsize=max(((diffdiam.*scalemarker).^2),1);
    set(ptplot,'SizeData',ultsize);
  drawnow

 end     %________ End of time loop ________________
toc
 
% Change particle plot to not use every point:

  figure(1)
  for ll=1:2:N
      plotmask=X(:,5)==ll;
      diffdiam=2*sqrt(2*X(plotmask,4)*dt);
      minsize=.000001; 
      ptplot=scatter(X(plotmask,1),X(plotmask,3),diffdiam);
  daspect([1 1 1]);
  axis([0 L 0 H]);
  title(['Diameter is approximate diffusion distance per dt=',num2str(dt)]);
  xlabel('X (m)');ylabel('Z (m)');
    ax=gca;
    currentunits=ax.Units;
    ax.Position;
    ax.Units='Points';
    axpos=ax.Position;
    pba=ax.PlotBoxAspectRatio;
    ax.Units=currentunits;
    scalemarker=min([axpos(3),pba(1)/pba(2)*axpos(4)])/diff(xlim);
    ultsize=max(((diffdiam.*scalemarker).^2),1);
    set(ptplot,'SizeData',ultsize);
    hold on
  end  
  hold off
  drawnow

 