function y = activationModelDrift6(A,dact,m,n,rho,ton,zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,t)
    global tmax driftmod fitting signalact driftmod2 tend driftdermod sigma fdrift tsignalact
    %start with drift
    tdriftswitch=ton;

    y=zeros(size(t));
    maxi=find(t==tend);
try    
if ~isempty(maxi) || length(t)>=1
        
    tin=t;
    if ~isempty(maxi) 
        t=t(1:maxi(1));
    end
    
    %drift + Hill function
    ind=t>ton&t<=tmax;
    dact0=dact;
    dact=dact*(tmax-ton);
     x=t(ind)-ton;
    %correct for numerical blunders
    if sum(ind)~=0 && ~isempty(x) && x(1)<0
        x(1)=0;
    end
    xn=x.^n;
    y(ind)=y(ind)+A*xn./(xn+(dact)^n)+m*cumtrapz(x,xn./(xn+(dact^n)));
    
    Atot=y(find(ind,1,'last'));
    
    t1=tmax;
    x0=(1-rho)*(t1-ton);
    t0=ton+x0;
    
    hill=(x0^n)/(dact^n+x0^n);
    xrho=x(x<=x0);
    
    if isempty(xrho)
        xrho=[0 x0];
    elseif xrho(end)~=x0
        xrho=[xrho;x0];
    end
    y0=A*hill+m*trapz(xrho,(xrho.^n)./((xrho.^n)+(dact^n)));
    y0p=(A*n*(x0^(n-1))*(dact^n)/(((x0^n)+(dact^n))^2))+m*((x0^n)/((x0^n)+(dact^n)));
    
    
    y1=Atot;
    y1p=0;

   
%   ind2=t>=t1&t<=t2;

     herm1=@(x) y0+(y0p*x)+((3*y1-3*y0-y1p-2*y0p)*x.^2)+((-2*y1+2*y0+y0p+y1p)*x.^3);
     
%     herm2=@(x) y1+(y1p*x)+((3*y2-3*y1-y2p-2*y1p)*x.^2)+((-2*y2+2*y1+y1p+y2p)*x.^3);

         ind1=find(t>=t0&t<=t1);
         
         if length(ind1)>1
            y(ind1)=herm1((t(ind1)-t0)/(t1-t0));
         end
    

    
%     if fitting && length(t)==length(signal)
%         cf(66);plot(t,signal,t,y);
%     end
if ~isempty(maxi)
    y(1:maxi(1))=y(1:maxi(1))+driftmod2(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,t);
else
     y=y+driftmod2(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,t);
end
    if length(maxi)>1
        ind=tin(maxi(1)+1:maxi(2))>ton&tin(maxi(1)+1:maxi(2))<=tmax;
        ind=find(ind)+maxi(1);
        
        x=tin(ind)-ton;
        xn=x.^n;
        denom=(xn+dact^n).^2;
        
        y(ind)=A*n*(x.^(n-1)).*((dact^n)./denom);
        y(ind)=y(ind)+m*(x.^n)./(xn+dact^n);
        
        y(maxi(1)+1:maxi(2))=y(maxi(1)+1:maxi(2))+driftdermod(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,tin(maxi(1)+1:maxi(2)));
        
        
                
        ind1=find(t>=t0&t<=t1)+maxi(1);

        
        herm1p=@(x) y0p+ (2*(3*y1-3*y0-y1p-2*y0p)*x)+(3*(-2*y1+2*y0+y0p+y1p)*x.^2);
       if length(ind1)>1
            y(ind1)=herm1p((tin(ind1)-t0)/(t1-t0))/(t1-t0);
       end
        
        if length(maxi)>2
            y(maxi(2)+1:maxi(3))=driftmod2(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,tin(maxi(2)+1:maxi(3)));


            if length(maxi)>3
                 y(maxi(3)+1:maxi(4))=driftdermod(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,tin(maxi(3)+1:maxi(4)));
            end
            
        if length(maxi)>4
%             y(maxi(4)+1:end)=abs(driftdermod(zA1,0,zm1,0,ztau1,Inf,ztdriftswitch,tin(maxi(4)+1:end))-driftdermod(0,zA2,0,zm2,Inf,ztau2,ztdriftswitch,tin(maxi(4)+1:end)));
             y(maxi(4)+1:end)=abs(driftmod2(zA1,0,zm1,0,ztau1,Inf,ztdriftswitch,0,tin(maxi(4)+1:end))-driftmod2(0,zA2,0,zm2,Inf,ztau2,ztdriftswitch,0,tin(maxi(4)+1:end)));

            y(maxi(4)+1:maxi(4)+find(tin(maxi(4)+1:end)>ztdriftswitch,1)-1)=0;
%             y(maxi(4)+1+find(tin(maxi(4)+1:end)>tmax,1):end)=sqrt(y(maxi(4)+1+find(tin(maxi(4)+1:end)>tmax,1):end));
        end

        end
    end
end
catch e
    getReport(e);
    
end
if sum(isnan(y))~=0 || sum(imag(y)~=0)~=0
    disp('nana')
end

% if ~isempty(maxi) 
%     plot(t,y(1:maxi(1)),t0,y(find(t>=t0,1)),'gx',tsignalact,signalact);
%     
% end
end

