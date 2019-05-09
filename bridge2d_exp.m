% Edited by Mark M Meerschaert on 14 April 2011
% Based on code by Hans-Peter Scheffler circa 2008
%
% Simulation code for fractal K field with operator scaling
% Code can produce a normal or stable (heavy tailed) log K field
% Field has stationary increments, like fractional Brownian motion
% Hurst index in x1 direction is H/a_1
% Hurst index in x2 direction is H/a_2
% IMPORTANT NOTE: The code requires 0 < H < a_i !!! 
%
% Input variable definitions
% alpha = 2 for normal noise, 0<alpha<2 for stable noise (heavy tail)
% H = Hurst index applied to filter
% A = Scale for noise variables (standard deviation if alpha=2)
% M = size of OUTPUT field (internal fields are 2M+1 by 2M+1)
% a1 = power law index of filter for x1 coordinate
% a2 = power law index of filter for x2 coordinate
% v1 = direction of x1 coordinate axes in radians 
% v2 = direction of x2 coordinate axes in radians 
% C1 = correlation length for x1 coordinate 
% C2 = correlation length for x1 coordinate  
% rho = 2 (see formula for filter phi)
% mu = sample standard deviation of log K field
% sigma = sample standard deviation of log K field
% nuse  = # of conditioning points to use for each generated point
clear all
% statistical parameters
 alpha = 2; %(2 for normal, 0<alpha<2 for stable)
 H = 0.4;   A = 1;  
 a1 = .9;    a2 = 2-a1;   a3=1;
 C1 =  1;    C2 =3;   C3=1; 
 rho = 2;   mu= 0;  sigma = .00000001; nuse=100;
 lam1=1; lam2=.1*lam1;  lammag=sqrt(lam1^2+lam2^2) % for exponential
 lammag=0.5/C1+0.5/C2;
 % angles associated with the scaling eigenvectors
 angle1=0*pi/36; angle2=18*pi/36;

 % Geometry of the output field
 origin=[0,0]; max=[25,100]; m = [32,128]; dx=(max-origin)./m 
 nreals=1;

 %unit vectors in the scaling eigenvector directions
 theta1=[cos(angle1),sin(angle1)]; theta2=[cos(angle2),sin(angle2)];
 
 % check to see if there is conditioning data
 if exist('conditioning.txt')==2;
  docondition=1
  cpts1 = load('conditioning.txt','-ascii'); cpts1(:,1)=ceil((cpts1(:,1)-origin(1))/dx(1));  
  cpts1(:,2)=ceil((cpts1(:,2)-origin(2))/dx(2)); 
  %temp=cpts1(:,1); cpts1(:,1)=cpts1(:,2); cpts(:,2)=temp;
  cpts1=sortrows(cpts1, [1 2]);  cpts=cpts1;
  % consolidate repeated points in a single cell by averaging the conditioning points
  ncp=length(cpts(:,1))
  k=ncp;
  while k>1;
    loc=cpts(k,1:2);
    count=1;
    value=cpts(k,3);
    while k-count>0 & cpts(k-count,1)==loc(1) & cpts(k-count,2)==loc(2);
        value=value+cpts(k-count,3);
        count=count+1;
    end;
    cpts(k,3)=value/count;
    if count>1 ;
        cpts(k-1:k-count+1,:)=[];
    end;
    k=k-count;
  end;
  ncp=length(cpts(:,1))
 else; docondition=0;  % no conditioning data
 end;    
 
  tic;
 
 % Make the Fourier kernel phi 
 
 Q=a1+a2;
        
 k1 = repmat(-m(1):m(1),2*m(2)+1,1);
 k2 = repmat(-m(2):m(2),2*m(1)+1,1);
 %k1 = k1*(A/m(1));  
 k1=k1';
 %k2=k2*(A/m(2));
 k1=k1/m(1)/dx(1); k2=k2/m(2)/dx(2);

        %d3 = ones(2*m+1,2*m+1,2*m+1);
        %n=1;
        %for j=-m:m;
        %d3 (:,:,n)= d3(:,:,n).*j;
        %n=n+1;
        %end
        
        %d2 = d1(:,:,1); d2=d2';
        %d2 = repmat(d2,[1,1,2*m+1]);
        
        % phi is the FFT filter
        % for exponential
        
        phi = lammag./(lammag^2/4 + (C1*(abs(k1*theta1(1) + k2*theta1(2)))).^2 +...
               (C2*(abs(k1*theta2(1) + k2*theta2(2)))).^2 ) ; 

%        phi(m(1)+1,m(2)+1)=m(1)*m(2)/dx(1)/dx(2);
%end
        phi=ifftshift(phi);

        cov=real(fftn(phi.*conj(phi)))/m(1)/m(2);
        cov=fftshift(cov);
        
 % interpolate the trend surface of the conditioning points
 if (docondition ==1); 
     trendc=krigs2d(cpts,cov,m,nuse); 
     %trendc= trend2d(cpts,cov,m,nuse); 
 figure(5),imagesc(trendc); hold on; plot(cpts(:,2),cpts(:,1),'o'); hold off;
 end;

 % generate nreals unconditional fields
 % simulate FFT of noise field using IID normal (stable if alpha<2)
 % This would be a good point to loop and make multiple realizations
 
for n=1:nreals 
 Z = randn(2*m(1)+1,2*m(2)+1);    Z = Z* ((A/m(1)/m(2))^(1/alpha));
 Z1 = randn(2*m(1)+1,2*m(2)+1);   Z1 = Z1* ((A/m(1)/m(2))^(1/alpha));
 Z = complex(Z,Z1);        
 Z(m(1),m(2)) = 0;   Z(m(1)+1,m(2))=0;   Z(m(1),m(2)+1)=0;    Z(m(1)+1,m(2)+1)=0;  
        %Z(m,m,m+1)=0;   Z(m+1,m,m+1)=0; Z(m,m+1,m+1)=0;  Z(m+1,m+1,m+1)=0;

         K = phi.*Z; % FFT of log K field obtained by multiplying FFT of noise field Z with FFT filter phi
       
        clear Z Z1;

        K = real(fftn(K)); % Invert the FFT to get the simulated field
%        Y(1:2:end)=-Y(1:2:end); % Correct for shifting FFT integral from [-A,A] to [0,2*A] by (-1)^{i+j} 

        %clear X;
        K=K(m(1)/2:3*m(1)/2-1,m(2)/2:3*m(2)/2-1);
        s=mean(std(K)); % average stdev for the simulated field (mean is zero)
        K=(sigma/s)*(K-mean(mean(K)))+mu; % convert to mean mu, standard deviation sigma
 
        if(docondition==1)
          newpts=cpts;   % make a copy of the conditioning points list so I can use the location data
          for k=1:ncp;
            newpts(k,3)=K(cpts(k,1),cpts(k,2));   % replace the random values with the new unconditioned field
          end;
        %   get trend surface of uncond. field
         % trendnew=trend2d(newpts,cov,m,nuse);
          trendnew=krigs2d(newpts,cov,m,nuse);
        %   make the new conditioned field
          K=trendc+(K-trendnew);
        end;  
 
        % plot the field
        figure(3), imagesc(K);
        axis tight; axis equal;
        pause(0.1)
          
        K=exp(K);    
        toc;
        % make some conditioning points for testing
        ncp=100; 
        condwrite=rand(ncp,2); condwrite(:,1)=m(1)*condwrite(:,1); condwrite(:,2)=m(2)*condwrite(:,2);
        for k=1:ncp
        condwrite(k,3)=K(ceil(condwrite(k,1)),ceil(condwrite(k,2))); 
        end
        condwrite(:,1)=dx(1)*condwrite(:,1); condwrite(:,2)=dx(2)*condwrite(:,2);
        %if(docondition==0); save('conditioning.txt','condwrite','-ascii'); end; 
        fname=strcat('kout_exp_',int2str(n))
        
        save (fname, 'K', '-ascii' );
end;        