function [fest,spike,imu,festfine,sigmaminus]= TVRegSpikeRemove( data, iter, u0, scale, ep, dx,nanti, plotflag, diagflag )
% u = TVRegDiff( data, iter, alph, u0, scale, ep, dx, plotflag, diagflag );
% u = tvdiff( e, dx, iter, ep, alph );
% Rick Chartrand (rickc@lanl.gov), Apr. 10, 2011
% Please cite Rick Chartrand, "Numerical differentiation of noisy,
% nonsmooth data," ISRN Applied Mathematics, Vol. 2011, Article ID 164564, 
% 2011. 
%
% Inputs:  (First three required; omitting the final N parameters for N < 7
%           or passing in [] results in default values being used.) 
%       data        Vector of data to be differentiated.
%
%       iter        Number of iterations to run the main loop.  A stopping
%                   condition based on the norm of the gradient vector g
%                   below would be an easy modification.  No default value.
%
%       alph        Regularization parameter.  This is the main parameter
%                   to fiddle with.  Start by varying by orders of
%                   magnitude until reasonable results are obtained.  A
%                   value to the nearest power of 10 is usally adequate.
%                   No default value.  Higher values increase
%                   regularization strenght and improve conditioning.
%
%       u0          Initialization of the iteration.  Default value is the
%                   naive derivative (without scaling), of appropriate
%                   length (this being different for the two methods).
%                   Although the solution is theoretically independent of
%                   the intialization, a poor choice can exacerbate
%                   conditioning issues when the linear system is solved.
%
%       scale       'large' or 'small' (case insensitive).  Default is
%                   'small'.  'small' has somewhat better boundary
%                   behavior, but becomes unwieldly for data larger than
%                   1000 entries or so.  'large' has simpler numerics but
%                   is more efficient for large-scale problems.  'large' is
%                   more readily modified for higher-order derivatives,
%                   since the implicit differentiation matrix is square.
%
%       ep          Parameter for avoiding division by zero.  Default value
%                   is 1e-6.  Results should not be very sensitive to the
%                   value.  Larger values improve conditioning and
%                   therefore speed, while smaller values give more
%                   accurate results with sharper jumps.
%
%       dx          Grid spacing, used in the definition of the derivative
%                   operators.  Default is the reciprocal of the data size.
%
%       plotflag    Flag whether to display plot at each iteration.
%                   Default is 1 (yes).  Useful, but adds significant
%                   running time.
%
%       diagflag    Flag whether to display diagnostics at each
%                   iteration.  Default is 1 (yes).  Useful for diagnosing
%                   preconditioning problems.  When tolerance is not met,
%                   an early iterate being best is more worrying than a
%                   large relative residual.
%                   
% Output:
%
%       u           Estimate of the regularized derivative of data.  Due to
%                   different grid assumptions, length( u ) = 
%                   length( data ) + 1 if scale = 'small', otherwise
%                   length( u ) = length( data ).

% Copyright notice:
% Copyright 2010. Los Alamos National Security, LLC. This material
% was produced under U.S. Government contract DE-AC52-06NA25396 for
% Los Alamos National Laboratory, which is operated by Los Alamos
% National Security, LLC, for the U.S. Department of Energy. The
% Government is granted for, itself and others acting on its
% behalf, a paid-up, nonexclusive, irrevocable worldwide license in
% this material to reproduce, prepare derivative works, and perform
% publicly and display publicly. Beginning five (5) years after
% (March 31, 2011) permission to assert copyright was obtained,
% subject to additional five-year worldwide renewals, the
% Government is granted for itself and others acting on its behalf
% a paid-up, nonexclusive, irrevocable worldwide license in this
% material to reproduce, prepare derivative works, distribute
% copies to the public, perform publicly and display publicly, and
% to permit others to do so. NEITHER THE UNITED STATES NOR THE
% UNITED STATES DEPARTMENT OF ENERGY, NOR LOS ALAMOS NATIONAL
% SECURITY, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY,
% EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
% RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF
% ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR
% REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED
% RIGHTS. 

% BSD License notice:
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met: 
% 
%      Redistributions of source code must rsigmaminusin the above
%      copyright notice, this list of conditions and the following
%      disclaimer.  
%      Redistributions in binary form must reproduce the above
%      copyright notice, this list of conditions and the following
%      disclaimer in the documentation and/or other materials
%      provided with the distribution. 
%      Neither the name of Los Alamos National Security nor the names of its
%      contributors may be used to endorse or promote products
%      derived from this software without specific prior written
%      permission. 
%  
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
% CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
% AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE. 

% code starts here
% Make sure we have a column vector.


data = data( : );
% Get the data size.
n = length( data );

% Default checking. (u0 is done separately within each method.)
if nargin < 9 || isempty( diagflag )
    diagflag = 1;
end
if nargin < 8 || isempty( plotflag )
    plotflag = 1;
end
if nargin < 7 || isempty( dx )
    dx = 1 / n;
end
if nargin < 6 || isempty( ep )
    ep = 1e-6;
end
if nargin < 5 || isempty( scale )
    scale = 'small';
end

% nanti=1;
expvar=1;
dd=abs(diff(data));
sigmaplus=sum(dd.*(1-dd/max(dd)))/sum(1-abs(dd)/max(abs(dd)));
alph=((sigmaplus*sigmaplus)^expvar)*dx^(nanti-1);
falph=ones(length(data),1);
falph(:)=alph;
pinf=0;

thresh=1;
threshrm=1;
taumin=1;

        % Construct antidifferentiation operator and its adjoint.
	str='v';
	for ii=1:nanti
		str=['cumsum(' str ')*dx']	;
	end
	AN=str2func(['@(v,dx)' str]);


	str='w';
	for ii=1:nanti
		str=['AT(' str ')*dx'];
	end


        if nargin < 4 || isempty( u0 )
%       u0 = [ 0; diff( data +0.1*rand(size(data)))/dx ];
        u0 = [ 0; diff( data )/dx ];
	    i=1;
	    while i<nanti 
	    
		u0=[0; diff( u0 )/dx;];
        i=i+1;
        end
        
             if plotflag
                subplot(2,1,1);
                plot(1:n,AN(u0,dx),'-r',1:n,data,'-b');
                subplot(2,1,2);
                plot( u0, 'ok' ), drawnow;
            end

        end
f=AN(u0,dx);
u = u0;
% res=zeros(size(u));

eta=exp(-(sum(diff(data))^2)/(2*sigmaplus^2));%   + 0.5;
%         eta=(eta-0.5)*2;

tau=taumin+eta*(thresh-taumin);
%                  tau=taumin+1/(1-eta);
phi=1;
chi=1;
        
beta=0.2;
sigmaminus=sigmaplus;

function calcDen8()
    
    
            nsr=sigmaplus/(max(AN(u,dx))-min(AN(u,dx)));
    
            wdu=(1-((abs(D*u))/max(abs(D*u))));
            wu=(1-((abs(u))/max(abs(u))));
        
                
%         ddp=[0;abs(diff(data)).*(1-abs(u(2:end))/max(abs(u(2:end)))).*wdu(2:end)]...
%             +[abs(diff(data)).*(1-abs(u(1:end-1))/max(abs(u(1:end-1)))).*wdu(1:end-1);0];
                
%         ddp=[sigmaminus;abs(diff(data)).*(1-abs(u(2:end))/max(abs(u(2:end)))+nsr)]/2+...
%             [abs(diff(data)).*(1-abs(u(1:end-1))/max(abs(u(1:end-1)))+nsr);sigmaminus]/2;
        
                ddp=[sigmaminus;abs(diff(data)).*(wdu(1:end-1)+nsr)]/2+...
            [abs(diff(data)).*(wdu(2:end)+nsr);sigmaminus]/2;


%         den8 = (abs((res)/(sigmaplus)).*(1-wdu))/(sum(1-wdu));

        
               dd2=[abs(diff(data));abs(diff(data(end-1:end)))];
        if exist('res')
            rsmax=(sum(abs((res)/(sigmaminus)).*(1-wdu))/(sum(1-wdu)));
            rsm=(sum(abs((res(~spike))/(sigmaminus)).*(1-wdu(~spike)))/(sum(1-wdu(~spike))));
        else
            rsm=(sum(abs((dd2)/(sigmaminus)).*(1-wdu(~spike)))/(sum(1-wdu(~spike))));
            rsp=(mean(abs((dd2)/(sigmaminus)).*(1-wdu(~spike)))/(sum(1-wdu(~spike))));
            rsmax=rsm;
        end
        if exist('res') 
%             if ~isempty(newnoise)
            ares=abs(res);
            
                    p=sum(ares<sigmaminus)/n;
            

            if sum(spike)~=0
                threshc=mean(ares(~spike)/sigmaplus)+std(ares(spike)/sigmaplus);
                threshrm=mean(ares(spike)/sigmaplus)-3*std(ares(spike)/sigmaplus);
            else
                threshc=1;
            end
%             plotDensity(ares(spike)/sigmaplus);
            thresh=max(thresh,threshc);
            
            threshrm=max(threshrm,threshc);
%             threshrm=min(abs(res(newnoise)/sigmaplus))/2
%             end
        else
            ares=sigmaplus;
              p=sum(ddp<sigmaplus)/n;
            threshc=1;
            thresh=2;
            threshrm=1;
        end
        
        tau=(1-eta)*taumin+(eta)*(thresh);
        sigmamax=min(sigmaplus,sigmaminus);
        neff=sum(~spike);
        eta=exp(-((sigmaminus-sigmamax)^2)/(2*sigmamax^2))*exp(-(abs(rsp)^2)*neff);
        
        qq=(rsm/rsmax)^(1/1);

        e8=((1-nsr^(tau/thresh))/(sqrt(1+eta)*(rsmax*(1-qq)+qq*rsm)));
        e8=(1-nsr^(1/rsm));
%         e8=(1+nsr^(tau/thresh))/(max(rsm,1));
%         e8=-1+1/eta;
%         e8=1;
        
%  
%         den8=ep^(threshrm)+(exp(-(abs(ddp-max(ddp)))/std(ddp))).^e8;
% %         den8(den8<ep^(e8))=ep^(e8);
        if exist('res')
            sigmalocal=filter(b,a,ddp);
        else
            sigmalocal=ddp;
        end
%             
%             pp=exp(-((abs(sigmalocal-max(sigmalocal)).^2)/var(sigmalocal))-(abs(u)/var(abs(u))));
%             pp=exp(-(abs(u.^2)/var(abs(u))));
%             ms=max(ddp(ddp<(mean(ddp)+3*std(ddp))));
%             jj=(1-abs(u)/max(abs(u)));
% 
%             kk=(abs(dd2)/max(abs(dd2)));
%             pp=exp(-((abs(sigmalocal-max(sigmalocal)).^2)/(ms)^2)-((u.^2)/var(u))).*kk.*jj;
%             pp=pp/max(pp);
% %             plot(pp);
%             ind=sigmalocal>sigmaplus&abs(u)<std(u);
% 
%     %         den8(ind)=den8(ind)+sigmalocal(ind);
% %             den8=(1-pp).*den8+pp.*(sigmalocal+sigmaplus);
%         den8=abs(u)/max(abs(u))+ep+sigmalocal/max(AN(u,dx));
        
        ee=mean(ddp(ddp<=(1-nsr^(1/rsmax))*(mean(ddp(~spike))-std(ddp(ddp<mean(ddp(~spike)))))));
        
        mean(ddp<mean(ddp));
        mmu=Inf;
        while std(ddp(ddp<mmu))>mean(ddp(ddp<mmu))/3
        mmu=mean(ddp(ddp<mmu));
        end
        ee = mean(ddp<mmu);
%         ee=mean(ddp(ddp<=-log((1-(p^(1+rsmax))))/mean(ddp)));
        
        
        
%         plotDensity(ddp);hold on;y=ylim();h=plot([ee ee],y,'-g');set(h,'Linewidth',1);hold off;
          den8=nsr.^((threshc*sigmaplus/rsmax)./ddp + std(u.^2)./(max(abs(u))-abs(u)).^2)+ee;%+(exp(-(abs(ddp-max(ddp)))/std(ddp))).^e8;
    ddp2=ddp;
    ind=ddp2<sigmalocal;
    ddp2(ind)=sigmalocal(ind);
     ff=std(u*dx)*sum(u>std(u)/2)/(max(AN(u,dx))-min(AN(u,dx)));
    uu=[u(2:end);u(end)]+(([wu(2:end);wu(end)]).^(2)).*[u(1);u(1:end-1)];
%     AAU = AN(AN(u,dx)-min(AN(u,dx)),1);
%     AAU(end)/(2*sqrt(2)*sigmaminus*n);
    ef=max((rsmax+1)*abs(diff(data))/(2*sigmaminus))-std((rsmax+1)*abs(diff(data))/(2*sigmaminus));
    
    ef=exp(-abs(ef).^(sqrt(rsmax)/2));
     den8=nsr.^((((0.5+(1/threshc))*sigmaplus/(1+rsmax))^1/2)./((sigmalocal+ddp).^1/2)+ mean(abs(u))./(max(abs(uu))-abs(uu)+nsr))...
         +exp(sqrt(2+rsmax)*log(nsr)*(abs([u(2:end);u(end)])/min(abs(u(u~=0)))).^(1/2))*(1/sigmaplus^(1+sqrt(nsr)))...
         + ef/sigmaplus + ee;%+(exp(-(abs(ddp-max(ddp)))/std(ddp))).^e8;
%      den8=1;
%     
%             den8=nsr.^((threshc*sigmaplus)./sigmalocal + mean((u))./(max(abs(u))-abs(u)+sigmaplus));%+(exp(-(abs(ddp-max(ddp)))/std(ddp))).^e8;

%           plot(1:n,log(den8),1:n,sigmalocal );
if plotflag
    plot(1:n,log(den8),1:n,ddp)
        1;
end
    
end



rs=0;
rsp=1;

 rsm=1;
 rsmax=1;

% Different methods for small- and large-scale problems.

     spike=false(size(data));
      

	ATN=str2func(['@(w,dx)' str]);

        %AT2 = @(w) AT(AT(w));
        % Construct differentiation matrix.
        c = ones( n, 1 );
        D = spdiags( [ -c c ], [ 0 1 ], n, n ) / dx;
        D( n, n ) = 0;
        
        Db = spdiags( [ -c c ], [ 1 0 ], n, n ) / dx;
        Db( n, n ) = 0;
        
        clear c
        DT = D';
        % Since Au( 0 ) = 0, we need to adjust.
        ofst=data( 1 );
        data = data - ofst;
        datap=data;
        % Default initialization is naive derivative.

        u = u0;
        uprev=u0;
        % Precompute.
        ATNd = ATN( data ,dx);
        
        weight=ones(size(data));
%         thresh=6;
%         weight=linspace(1,0,n)';
            NRG=@(z) alph*sum(abs(diff(z)*dx))+0.5*sum((weight.*(AN( z,dx )-data)*dx).^2);
            NRGP=@(z) alph*weight.*[0; abs(diff(z)*dx)]+0.5*(weight.*(AN( z,dx )-data)*dx).^2;
            E0=NRG(u);
            E=E0;
        % Main loop.
        ii=1;
        iii=1;
        tic;
        
        [mu,imu]=max(u);
        su=std(u);
        detect=true;
        gamma=0;
        abort=false;
        while ii<=iter && ~abort
            if toc>30
              break;
            end
            % Diagonal matrix of weights, for linearizing E-L equation.
        windowSize = 3;
        
            ep=sigmaplus/(max(AN(u,dx))-min(AN(u,dx)));
            ep=ep^2;
            

            b = (1/windowSize)*ones(1,windowSize);
            a = 1;
%             epu=filter(b,a,abs(D * u))+(ep*max(abs(u)))^2;
                aDu=abs(D*u);
                
                epu=filter(b,a,abs(D * u))+(ep*max(abs(u)))^2;
                
                epu=min(aDu(aDu~=0));
                epu=std(aDu(aDu<=sqrt(ep)*max(aDu)));
                
%             Q = spdiags( 1./ (sqrt( (D * u) .^2 ) + epu), 0, n, n );
            if iii>1
            p=(1-eta^2)+(eta^2)*pinf;
            p=1;
%             peff=(1-peff)*p+peff;
            
%             cf(99);plot(peff);
            else
                p=1;
                peff=1;
            end
            Q = spdiags( (abs( (D * u) ) + epu ).^(p-2), 0, n, n );
            % Linearized diffusion matrix, also approximation of Hessian.
            try
            L = dx * DT * Q * D;
            catch e
                disp(e);
            end
            % Gradient of functional.
            %
%             gamma=((1-alph^iii)*exp((1/iter)^2))*exp(-(iii/iter)^2)+alph^iii
                dd=abs(diff(data(~spike)));
                sigmaplus=sum(dd.*(1-dd/max(dd)))/sum(1-abs(dd)/max(abs(dd)));
%                 res=(datap-AN(u,dx));
                if iii>1
                    falph=abs(res);
                    windowSize = 3;

                    b = (1/windowSize)*ones(1,windowSize);
                    a = 1;
                    falph = filter(b,a,falph);
%                     sigmaminus=mean(falph)
                    falph=falph;
%                     falph=((falph.*falph).^expvar)*dx^(nanti-1);
                    falph=alph*exp(-((max(abs(res))-abs(res)).^2)/var(abs(res)));
%                     falph=(falph+alph);
                    mfalph=mean(falph);
                    falph(:)=alph;
                    
                    
          if plotflag
                subplot(2,1,1);
                h=plot(1:n,data,'-b',1:n,data,'.b',1:n,AN(u,dx),'-r',find(spike),data(spike),'og',imu,data(imu),'vb',1:n,min(data)+falph*(max(data)-min(data))/max(falph));
                set(h(3),'linewidth',2);
                legend('Regularized Estimate','Data','Data','Spikes');
                subplot(2,1,2);
                plot(1:n, u, 'ok',imu,mu,'vb'), drawnow;
            end
                    
%                     falph(:)=alph;
%              del0=(1-(abs(sigmaminus-sigmaplus)/max(sigmaminus,sigmaplus)))^2;
%                 del0=((1-eta^2)+(1-chi))/2 +eta^2;
                           

                else
                    del0=1;
                end
                
%                  del=min(del*del,1)     


             if iii>1
                gamma=sigmaplus;
%                 del=1-del

%                     gamma=sigmaplus*(sigmaplus*(1-phi)*(1-eta))+(sigmaplus*(eta^2));
                    gamma=sigmaplus*((1-eta)^2)+sigmaminus*sigmaminus*eta;
                    
             else
                 gamma=sigmaplus*((1-eta)^2)+sigmaminus*sigmaminus*eta;
             end
             
             gamma=(sigmaplus^1.5)*eta+sigmaplus*sigmaplus*(1-eta);
             gamma = sigmaminus*sigmaminus*p;
%              del=del0
             gamma=0;
             del=del0;
             del=1;
             
%              gamma=0
%                 if gamma==0
%                     disp('its over!!!')
%                 end
% %              del=del
% %              gamma=0;
% %            gamma=alph^(iii/iter)
%             delta=0;
            if iii>1
                res=AN( u ,dx)-data;
                ep1=std(abs(res(~spike)));
                b = (1/windowSize)*ones(1,windowSize);
                a = 1;
                sigmalocal=filter(b,a,abs(res));
                sigmaplocal=filter(b,a,dd.*(1-dd/max(dd)));
                ep2=filter(b,a,abs(res));
                ep2=ep*min(ep2(ep2~=0));
%                 ep2(ep2>ep1)=ep1;
            else
                ep2=ep;
            end
            
          
            den3=1;
            
            epd=min(abs(u(u~=0)))/max(abs(u));
            abu=abs(filter(b,a, u));
            mxu=max(abu);
            abDu=abs(filter(b,a,D * u));
            

            mxDu=max(abDu);
            epDu=min(abDu(abDu~=0));
            
            den4=(abu/mxu);
            
            den5=exp(-(((abDu).^2))/(2*var(abDu)));
            den6=1-exp(-(((mxu-abu).^2))/(2*var(mxu-abu)));
            
%            delta=abs[data(1:end)-[data(2:end); data(end)] +]
%             den5=den5/max(den5);
            den2=(den4.*den5.*den6.*abs( (AN( u ,dx)-data)));
            den2=den2+min(dd(dd~=0));
            
            den2=abs( (AN( u ,dx)-data))+min(dd(dd~=0));
            
%             min(dd(dd~=0))
            if iii>1
%                 den8=sigmalocal0/sigmaplus+ep;
%                 plot(den8);
            else                 
                calcDen8()
            end
            
            if iii>1
%               denu=filter(b,a,abs(res));
%               denu=denu/sigmaplus;
%               denu=denu/max(denu);
%               denu=1-denu+max(min(abs(res(res~=0))),ep*max(AN(u,dx)));
              
            else
                denu=1;
            end
            
            
%             cf(777);plot(1:n,den4,'-g',1:n,den5,'-r',1:n,den6,'-c');legend('den4','den5','den6');
            %%mulitply den2 by something which is 1 most of the time and
            %%very small at the maximum of u or possibly Du
            
            
            
%              den2=(abs( weight.*(AN( u ,dx)-data))+ep2);
%             weight=weight;
% %             weight(weight==0)=0;
%            wdu=(exp(-((mean(abs(D*u))-abs(D*u)).^2)/var(abs(D*u))));
%         den2=(wdu+ep).*(abs( weight.*(AN( u ,dx)-data))+ep2);
            
%            rs=(sum(res(~spike).*wdu(~spike))/(sum(wdu(~spike))/sum(~spike)));

            if diagflag

                    fprintf('sigmaplus= %f, sigmaminus = %f, phi = %f, chi = %f , rsp = %f, rsm = %f\n',[sigmaplus sigmaminus phi chi rsp rsm])
                    fprintf('alph = %f, gamma = %f, eta = %f tau = %f \n',[alph gamma eta tau])

            end

            w2=weight;
            g = gamma*(ATN( w2.*(AN( u ,dx)-data)./den2,dx ))+del.*(ATN( w2.*(AN( u ,dx)-data)./(den8),dx ));
            g= g + falph .* (L * u).*denu;
            % Build preconditioner.
            c = cumsum( n : -1 : 1 ).';
            B = alph * L + spdiags( c( end : -1 : 1 ), 0, n, n );
            %droptol = 1.0e-2;
            try
                R = ichol( B, struct('droptol',1.0e-3) );
            catch
                alpha = max(sum(abs(B),2)./diag(B))-2;
                R = ichol( B, struct('droptol',1.0e-3,'diagcomp',alpha) );
            end
            % Prepare to solve linear equation.
            tol = nsr;
            maxit = int16(sqrt(n))*n;
            
            hs = @(x) ( alph *( L * x).*denu + gamma*(ATN( w2.*(AN( x ,dx))./den2,dx )) + del.*ATN( w2.*AN( x,dx )./den8,dx ) );
            
            if diagflag
                [s,flag,relres,ip] = pcg( hs , -g, tol, maxit, R', R );
                
                if flag==0
                    disp(['pcg converged at iteration ' int2str(ip) ' to a solution with relative residual ' num2str(relres) '.']);
                else
                    
          
                    disp(flag);
                    abort=true;
                    
                    s = pcg( hs , -g, tol, maxit, R', R );
                end
                
 
                
                fprintf( 'iteration %2d: relative change = %.3e, gradient norm = %.3e, gradient elm max = %.3e, dF =  %.3e\n', ii, norm( s ) / norm( u ), norm( g ),max(abs(g)),NRG(u+s)-NRG(u));
            else
               [s,flag,relres,ip] = pcg( hs , -g, tol, maxit, R', R );
                
            
            end
            
            if norm( s ) / norm( u ) < (nsr)^(3/4)
                abort=true;
            end
            % Update current solution
            u = u + s;
            
%             uprev=u;
            
            if ~diagflag
                hess = ( alph * L * (u) + ATN( AN( u,dx ),dx ) );
            end
            E=NRG(u);
             res=(data-AN(u,dx));
            
            % Display plot.
%             exp(-(res(res>0)/(thresh*sigma)));

%             weight(res>thresh*sigma)=0;

          if plotflag
                subplot(2,1,1);
                h=plot(1:n,data,'-b',1:n,data,'.b',1:n,AN(u,dx),'-r',find(spike),data(spike),'og',imu,data(imu),'vb',1:n,min(data)+falph*(max(data)-min(data))/max(falph));
                set(h(3),'linewidth',2);
                legend('Regularized Estimate','Data','Data','Spikes');
                subplot(2,1,2);
                plot(1:n, u, 'ok',imu,mu,'vb'), drawnow;
          end
         wdu=(1-((abs(D*u))/max(abs(D*u))));
           rsp=(mean(abs(res(~spike))/(sigmaplus).*(1-wdu(~spike)))/(sum(1-wdu(~spike))));
          
        sigmamax=min(sigmaplus,sigmaminus);
        neff=sum(~spike);
        eta=exp(-((sigmaminus-sigmamax)^2)/(2*sigmamax^2))*exp(-(abs(rsp)^2)*neff);
        nsr=(1/(((max(AN(u,dx))-min(AN(u,dx)))/sigmaminus)));
        pinf= (1-nsr)^rsm;

%         eta=(eta-0.5)*2;
        tau=(1-eta)*taumin+(eta)*(thresh);
        
        ares=abs(res);
        ares(spike)=mean(ares(~spike));
        sigmalocal=filter(b,a,ares);
        sigmalocal0=sigmalocal;
        peff=exp(-(((sigmalocal-sigmaplus).^2)/(sigmaplus^2)));
        
        sigmalocal(sigmalocal>(mean(sigmalocal)+3*std(sigmalocal)))=(mean(sigmalocal)+3*std(sigmalocal));
        sigmalocal(sigmalocal<sigmaplus)=sigmaplus;
        
        wdu=(1-((abs(D*u))/max(abs(D*u))));
        rsmax=(sum(abs((res)/(sigmaminus)).*(1-wdu))/(sum(1-wdu)));
%         thresh=max(thresh,2*rsmax);
        rsm=(sum(abs((res(~spike))/(sigmaminus)).*(1-wdu(~spike)))/(sum(1-wdu(~spike))));

        neff=sum(~spike);
        eta=exp(-((sigmaminus-sigmamax)^2)/(2*sigmamax^2))*exp(-(abs(rsp)^2)/n);
        tau=(1-eta)*taumin+(eta)*(thresh);
        sigmamax=min(sigmaplus,sigmaminus);
        
%                     weight(spike)=sqrt(exp(-((abs(res(spike))-sigmaplus).^2)./(2*(tau*sigmalocal(spike)).^2)));
%             cf(88);plot(weight)
        

        
%         sigmalocal(spike)=mean(sigmalocal(~spike));
        
%         cf(888);plot(1:n,tau*sigmalocal,1:n,abs(res))
%                      tau=taumin+1/(1-eta);
%          indnew=~spike&abs(res)>tau*sigmaplus;
%         tau=3;
%          indnew=~spike&(res>tau*sigmalocal|-res>max(1,tau^2)*sigmalocal);
    
         indnew=~spike&(res>tau*sigmaplus|-res>max(1,tau^2)*sigmaplus);

         indrm=spike&(res<threshrm*sigmaminus|-res>max(1,threshrm^2)*sigmaminus);
         indnem(indrm)=false;
         rem = sum(indrm);
         
        if sum(spike)~=0 && detect
            newold=indnew|spike;
            if spike(end) || spike(1)
                if ~spike(1)
                 uclean=interp1([find(~spike);n],[u(~spike);u(end-1)],1:n);
                elseif ~spike(end)
                     uclean=interp1([1;find(~spike)],[u(2);u(~spike)],1:n);
                else
                    uclean=interp1([1;find(~spike);n],[u(2);u(~spike);u(end)],1:n);
                end
            else
             uclean=interp1(find(~spike),u(~spike),1:n);
            end
%             uclean(newold)=0;
try
          [mu,imu]=findpeaks(uclean);
catch e
    disp(e)
end
%             if ~isempty(find(imu==24))
%                 disp('24 is in the race');
%             end
            ind=mu>std(uclean(abs(uclean)<1.3236*std(uclean)/sqrt(2)));%&mu;%>=beta*max(uclean)&(data(imu)-(data(imu-1)+u(imu-1)*dx)<tau*sigmaplus | data(imu)-(data(imu+1)-u(imu)*dx)<tau*sigmaplus)'&(~spike(imu+1) | ~spike(imu-1))';
            mu=mu(ind);
           imu=imu(ind);
          
        else
            try
            [mu,imu]=findpeaks(u);
            catch e
                disp(e)
            end
             ind=mu>std(u(abs(u)<1.3236*std(u)/sqrt(2)));%&mu>=beta*max(u)&(data(imu)-data(imu-1)<tau*sigmaplus | data(imu)-data(imu+1)<tau*sigmaplus);

            mu=mu(ind);
           imu=imu(ind);
           newold=spike;
        end

     jj=1;
       while jj<=length(imu)


%            if indnew(imu(jj))
%                 indnew(imu(jj))=false;
%            end
%             if imu(jj)<n-1 && newold(imu(jj)+1) && ~indnew(imu(jj)+2)
%                 kk=imu(jj)+1;
%                 while kk<n
%                     if newold(kk) 
% %                         if data(kk)<data(kk+1)  || u(kk+1)*(dx^nanti)>-tau*sigmaplus 
%                         if  u(kk+1)*(dx^nanti)>-tau*sigmaplus 
%                             spike(kk)=false;
%                             indnew(kk)=false;
%                             weight(kk)=1;
%                         end
%                     else
%                         kk=n;
%                     end
%                     kk=kk+1;
%                 end
%             end
%             if imu(jj)>2 && newold(imu(jj)-1) && ~indnew(imu(jj)-2)
%                 kk=imu(jj)-1;
%                 while kk>1
%                     if kk==24
%                         disp('this better not pass');
%                     end
%                     if newold(kk)
% %                         if data(kk)<data(kk-1) || u(kk-1)*(dx^nanti)<tau*sigmaplus
%                         if  u(kk-1)*(dx^nanti)<tau*sigmaplus
%                             spike(kk)=false;
%                             indnew(kk)=false;
%                             weight(kk)=1;
%                         end
%                     else
%                         kk=1;
%                     end
%                     kk=kk-1;
%                 end
% 
%             end
        jj=jj+1;
       end
        newnoise=find(indnew);
        pred=AN(u,dx);
        if (~isempty(newnoise)||sum(indrm)~=0)  && detect
        ii=0;

%             if ~isempty(find(newnoise==24))
%                 disp('not today!')
%             end
            spike(newnoise)=true;
            weight(newnoise)=0;
            
            spike(indrm)=false;
            weight(indrm)=1;
            
            
            if plotflag
                subplot(2,1,1);
                h=plot(1:n,data,'-b',1:n,data,'.b',1:n,AN(u,dx),'-r',find(spike),data(spike),'og',imu,data(imu),'vb',1:n,min(data)+falph*(max(data)-min(data))/max(falph));
                set(h(3),'linewidth',2);
                legend('Regularized Estimate','Data','Data','Spikes');
                subplot(2,1,2);
                plot(1:n, u, 'ok',imu,mu,'vb'), drawnow;
            end
            
            datap=AN(u,dx);
            
            for jj=1:n
                if spike(jj) && jj~=1 && jj~=n
                    st=jj-1;
                    jj=jj+1;
                    while jj<n && spike(jj) 
                        jj=jj+1;
                    end
                    en=jj;

                      datap(st:en)=linspace(datap(st),datap(en),en-st+1);

                end
                

            end

%             datap(~noise)=pred(~noise);
            
            u0p = [ u(1); diff( datap )/dx ];
%             i = 1;
%             while i<nanti 
% 
%             u0p=[0; diff( u0p )/dx;];
%               i=i+1;
%             end

%         alph=((sigmaplus*sigmaplus)^expvar)*dx^(nanti-1)*((abs(sigmaminus-sigmaplus)/max(sigmaminus,sigmaplus))^2)*(max(abs(u*dx))/max(abs(diff(data(~spike)))))^2;
        

        
            if plotflag
                subplot(2,1,1);
                h=plot(1:n,data,'-b',1:n,data,'.b',1:n,AN(u0p,dx),'-r',find(spike),data(spike),'og',imu,data(imu),'vb',1:n,min(data)+falph*(max(data)-min(data))/max(falph));
                set(h(3),'LineWidth',2);
                legend('Regularized Estimate','Data','Spikes');
                subplot(2,1,2);
                plot(1:n, u, 'ok',imu,mu,'vb',1:n,u0p,'ob' ), drawnow;
            end
            
                        
            u=u0p;
        else
            if diagflag
                disp('no new noise');
            end
%             expvar=1;
            if ii==iter && detect
                ii=0;
                expvar=0.75;
                detect=false;
                ufine=u;
            end
            


                        % Display plot.
            if plotflag
                subplot(2,1,1);
                h=plot(1:n,data,'-b',1:n,data,'.b',1:n,AN(u,dx),'-r',find(spike),data(spike),'og',imu,data(imu),'vb',1:n,min(data)+falph*(max(data)-min(data))/max(falph));
                set(h(3),'linewidth',2);
                subplot(2,1,2);
                plot( 1:n,u, 'ok',imu,mu,'vb' ), drawnow;
            end
        end
        
        
        
        sigmaminus=mean(abs(res(~spike)));
        
  
                            
        dd=abs(diff(data(~spike)));
        
        sigmaplus=sum(dd.*(1-dd/max(dd)))/sum(1-abs(dd)/max(abs(dd)));

%         ddp2=;

        
        phi=(max(abs(u*dx))/max(abs(diff(data(~spike)))));
        chi=(abs(sigmaminus-sigmaplus)/max(sigmaminus,sigmaplus));
        
%         wdu=(1-((abs(D*u))/max(abs(D*u))));
% 
%         
%                 
% %         ddp=[0;abs(diff(data)).*(1-abs(u(2:end))/max(abs(u(2:end)))).*wdu(2:end)]...
% %             +[abs(diff(data)).*(1-abs(u(1:end-1))/max(abs(u(1:end-1)))).*wdu(1:end-1);0];
%                 
%         ddp=[0;abs(diff(data)).*(1-abs(u(2:end))/max(abs(u(2:end))))]...
%             +[abs(diff(data)).*(1-abs(u(1:end-1))/max(abs(u(1:end-1))));0];
% 
% 
% %         den8 = (abs((res)/(sigmaplus)).*(1-wdu))/(sum(1-wdu));
% %         den8(den8<sqrt(ep))=sqrt(ep);
%         den8=((((abs(u))/max(abs(u)))).^2);
% %         den8(den8<ep)=ep;
% 
% 
% %         den8(den8==0)=min(den8(den8~=0));
%         
%         
% %         
% %         if (sum(abs((res(~spike))/(sigmaminus)).*(1-wdu(~spike)))/(sum(1-wdu(~spike))))>rsm
% %             break;
% %         else
% %             rsm=(sum(abs((res(~spike))/(sigmaminus)).*(1-wdu(~spike)))/(sum(1-wdu(~spike))));
% %         end


%         rsm=(sum(abs((res(~spike))/(sigmaminus)).*(1-wdu(~spike)))/(sum(1-wdu(~spike))));
%         thresh=max(thresh,2*rsmax);
%         tau=(1-eta)*taumin+(eta)*(thresh);
%         sigmamax=min(sigmaplus,sigmaminus);
%         eta=exp(-((sigmaminus-sigmamax)^2)/(2*sigmamax^2))*exp(-abs(rsp-1)^2);
%         
%         qq=(rsm/rsmax)^(1/1);
%         e8=((1+nsr^((1-eta)*tau/thresh))/(sqrt(1+eta)*(rsmax*(1-qq)+qq*rsm)));
% %         e8=-1+1/eta;
% %         e8=1;
%         sigmalocal=filter(b,a,abs(res));
%  
%         den8=(exp(-(abs(ddp-max(ddp)))/std(ddp))+ep).^(e8);
% %         den8(den8<ep^(e8))=ep^(e8);
%         
%         pp=exp(-((abs(sigmalocal-max(sigmalocal)).^2)/var(sigmalocal))-(abs(u)/var(abs(u))));
%         pp=exp(-(abs(u.^2)/var(abs(u))));
%         ms=max(ddp(ddp<(mean(ddp)+3*std(ddp))));
%         jj=(1-abs(u)/max(abs(u)));
%          dd2=[abs(diff(data));abs(diff(data(end-1:end)))];
%         kk=(abs(dd2)/max(abs(dd2)));
%         pp=exp(-((abs(sigmalocal-max(sigmalocal)).^2)/(ms)^2)-((u.^2)/var(u))).*kk.*jj;
%         pp=pp/max(pp);
%         plot(pp);
%         ind=sigmalocal>sigmaplus&abs(u)<std(u);
%         
% %         den8(ind)=den8(ind)+sigmalocal(ind);
%         den8=(1-pp).*den8+pp.*(sigmalocal+sigmaplus);
%         
% %         den8=(1-(abs(ddp-max(ddp)))/std(ddp).^(e8));
% %         den8=(abs(ddp)/max(ddp)).^(e8);
%         
%        
%         plot(1:n,log(den8));
%         legend('den8','pp','sigmalocal')
        %                   rs2=(mean(abs(res(~spike))));

        

%             alph=sigmaplus*((1-eta)^2)*phi + sigmaplus*sigmaplus*phi*chi+sigmaminus*sigmaminus*eta;
        calcDen8();
        alph=sigmaplus*sigmaplus;
        
        sigmaminus=mean(abs(res(~spike)));
         
%          tau=taumin+1/(1-eta);
%         sigma=sigmaminus
%         sigma=mean(abs(res(~noise)));
%         alph=(mean(abs(res(~spike)))^expvar)*dx^(nanti-1);

        ii=ii+1;
        iii=iii+1;
        end
        



fest=zeros(length(data),1+nanti);
fest(:,1+nanti)=u;
ii=nanti;
while ii>0;
fest(:,ii)=cumsum(fest(:,ii+1))*dx;
ii=ii-1;
end
fest(:,1)=fest(:,1)+ofst;


festfine=zeros(length(data),1+nanti);
if exist('ufine','var')==0;
    ufine=u;
end
festfine(:,1+nanti)=ufine;
ii=nanti;
while ii>0;
festfine(:,ii)=cumsum(festfine(:,ii+1))*dx;
ii=ii-1;
end
festfine(:,1)=festfine(:,1)+ofst;

            if plotflag
                subplot(2,1,1);
                h=plot(1:n,data,'-b',1:n,data,'.b',1:n,AN(u,dx),'-r',find(spike),data(spike),'og');
                set(h(3),'linewidth',2);
                subplot(2,1,2);
                plot(1:n, u, 'ok',1:n,u,'ob' ), drawnow;
            end


end

% function w = chop( v )
% w = v( 2 : end );
