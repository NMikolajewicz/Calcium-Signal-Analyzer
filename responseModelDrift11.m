function y = responseModelDrift11(Aoff,Aon,dact,dde,mact,nact,nde,rhoact,rhode,tdecay,ton,zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,t)
    global tmax fitting signal noise dff dfb driftmod driftmod2 driftdermod tend sigma fdrift indactmin inddrift inddriftclean indactminclean
    global wint
    %     t=varargin{end};
    %allocate y
    try
    if ~isempty(regexp(lastwarn,'starting'))
        disp(lastwarn)
    end
    y=zeros(size(t));
    maxi=find(t==tend);
    
if ~isempty(maxi) || length(t)>=1
        
    tin=t;
    if ~isempty(maxi)
        t=t(1:maxi(1));
    end
    

    %linear + Hill function (activation)
    dact=dact*(tdecay);

    ind=t>=ton&t<(ton+tdecay);
    xact=(t(ind)-ton);
    %correct for numerical blunders
    if sum(ind)~=0 && length(xact)~=0 && xact(1)<0
        xact(1)=0;
    end
    if isempty(xact)
        xact=[0,tdecay];
        xnact=xact.^nact;
        hill=xnact./(xnact+(dact^nact));
        
    else
        if xact(1)==0
            xnact=xact.^nact;
            hill=xnact./(xnact+(dact^nact));

            y(ind)=Aon*hill+mact*cumtrapz(xact,hill);
        else
            xact=[0;xact];
            xnact=xact.^nact;
            hill=xnact./(xnact+(dact^nact));
            try
            y((find(ind,1)-1):find(ind,1,'last'))=Aon*hill+mact*cumtrapz(xact,hill);
            catch e
                getReport(e);
            end
        end
        
     
    end
    
    hillend=(tdecay^nact)/(dact^nact+tdecay^nact);
    A=Aon*hillend+mact*trapz(xact,hill);

    t1=ton+tdecay;
    
    %linear + Hill function
    ind=t>=(t1);
    dde0=dde;
    dde=dde*(tend-ton-tdecay);
    idecay=find(t>=ton+tdecay,1);
    ide=find(t>=ton+tdecay+dde,1);
    
    if ~isempty(ide) && ~(ide>=idecay+1)
        dde=t(idecay+2)-ton-tdecay;
    end


    xde=t(ind)-t1;
    %correct for numerical blunders
    if sum(ind)~=0 && ~isempty(xde) && xde(1)<0
        x(1)=0;
    end
    
    
%     c=Aon*(tdecay^nact/(tdecay^nact+dact^nact))-Aoff;
%     y(ind)=y(ind)+Aoff*(1-x./(x+(dde)^nde))+c;

    
%     
%     Aoff0=Aoff;
%     xnde=xde.^nde;
%     hillde=1-(xnde./(xnde+(dde^nde)));
%     intlinde=dde/trapz(xde,hillde);
%     mde=-A*(1-Aoff0)/intlinde;
    
    

    
        Aoff0=Aoff;
        Aoff=Aoff*A;
        c=A-Aoff;


      if isempty(xde)
        xde=[0,tend];
        xnde=xde.^nde;
        hillde=1-(xnde./(xnde+(dde^nde)));
        intlinde=trapz(xde,hillde);
        mde=-A*(1-Aoff0*(1-hillde(end)))/intlinde;
        

     
        
        
     else
        if xde(1)==0
            
            xnde=xde.^nde;
            hillde=1-(xnde./(xnde+(dde^nde)));
            intlinde=trapz(xde,hillde);
            
            mde=-A*(1-Aoff0*(1-hillde(end)))/intlinde;

        
            linde=(mde*cumtrapz(xde,hillde));
            y(ind)=Aoff*(hillde)+linde+c;
            
        else
            xde=[0;xde];
            xnde=xde.^nde;
            hillde=1-(xnde./(xnde+(dde^nde)));
            intlinde=trapz(xde,hillde);
            mde=-A*(1-Aoff0*(1-hillde(end)))/intlinde;

        
            linde=(mde*cumtrapz(xde,hillde));
            y((find(ind,1)-1):find(ind,1,'last'))=Aoff*(hillde)+linde+c;
%             y((find(ind,1)-1):find(ind,1,'last'))=Aon*hill+mact*cumtrapz(xact,hill);
        end
        
     
    end
    
    

    x0=(1-rhoact)*tdecay;
    t0=ton+x0;
    
    hill=(x0^nact)/(dact^nact+x0^nact);
    xrhoact=xact(xact<=x0);
    
    if isempty(xrhoact)
        xrhoact=[0 x0];
    elseif xrhoact(end)~=x0
        xrhoact=[xrhoact;x0];
    end
    y0=Aon*hill+mact*trapz(xrhoact,(xrhoact.^nact)./((xrhoact.^nact)+(dact^nact)));
    y0p=(Aon*nact*(x0^(nact-1))*(dact^nact)/(((x0^nact)+(dact^nact))^2))+mact*((x0^nact)/((x0^nact)+(dact^nact)));
    
    
    y1=A;
    y1p=0;

    x2=(2*dde)*rhode;
    t2=t1+x2;
    hill=1-(x2^nde)/(dde^nde+x2^nde);
    xrhode=xde(xde<=x2);
    
    if isempty(xrhode)
        xrhode=[0 x2];
    elseif xrhode(end)~=x2
        xrhode=[xrhode;x2];
    end
    
    y2=Aoff*hill+mde*trapz(xrhode,1-((xrhode.^nde)./((dde^nde)+(xrhode.^nde))))+c;
    y2p=-(Aoff*nde*(x2^(nde-1))*(dde^nde)/(((x2^nde)+(dde^nde))^2))+mde*(dde^nde)/((x2^nde)+(dde^nde));
    
%     
    y0p=y0p*(t1-t0);
    y2p=y2p*(t2-t1);
    

    
    ind1=find(t>=t0&t<=t1);
    ind2=find(t>=t1&t<=t2);

     herm1=@(x) y0+(y0p*x)+((3*y1-3*y0-y1p-2*y0p)*x.^2)+((-2*y1+2*y0+y0p+y1p)*x.^3);
     
    herm2=@(x) y1+(y1p*x)+((3*y2-3*y1-y2p-2*y1p)*x.^2)+((-2*y2+2*y1+y1p+y2p)*x.^3);

     if length(ind1)>1
        y(ind1)=herm1((t(ind1)-t0)/(t1-t0));
     end
     if length(ind2)>1
        y(ind2)=herm2((t(ind2)-t1)/(t2-t1));
     end

% if length(t)==length(signal) && fitting
%     if sum(imag(y))~=0
%         disp('complex value')
%     end
%         cf(66);plot(t(1:maxi(1)),y(1:maxi(1)));
% end


y(y<0)=0;
if ~isempty(maxi)
    idwindle=find(y(1:maxi(1))>=sigma/2,1,'last');
    y(1:maxi(1))=y(1:maxi(1))+driftmod2(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,t);
else
    idwindle=find(y>=sigma/2,1,'last');
    y=y+driftmod2(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,t);
end

if length(maxi)>1&&fitting
    
        ind=find(t>ton&t<t1)+maxi(1);
        
        x=tin(ind)-ton;
        xn=x.^nact;
        denom=(xn+dact^nact).^2;
        
        y(ind)=Aon*nact*(x.^(nact-1)).*((dact^nact)./denom);
        y(ind)=y(ind)+mact*(xn)./(xn+dact^nact);
        
        
        ind=find(t>=t1)+maxi(1);
        
        x=tin(ind)-t1;
        xn=x.^nde;
        denom=(xn+dde^nde).^2;
        

        
        y(ind)=-Aoff*nde*(x.^(nde-1)).*((dde^nde)./denom);
        y(ind)=y(ind)+mde*(dde^nde)./(xn+dde^nde);
        
        ind1=find(t>=t0&t<=t1)+maxi(1);
        ind2=find(t>=t1&t<=t2)+maxi(1);

        
        herm1p=@(x) y0p+ (2*(3*y1-3*y0-y1p-2*y0p)*x)+(3*(-2*y1+2*y0+y0p+y1p)*x.^2);
        herm2p=@(x) y1p+ (2*(3*y2-3*y1-y2p-2*y1p)*x)+(3*(-2*y2+2*y1+y1p+y2p)*x.^2);
        
        if length(ind1)>1
        y(ind1)=herm1p((tin(ind1)-t0)/(t1-t0))/(t1-t0);
        end
        if length(ind2)>1
            y(ind2)=herm2p((tin(ind2)-t1)/(t2-t1))/(t2-t1);
        end
%         y(maxi(1)+1:maxi(2))=y(maxi(1)+1:maxi(2))+driftdermod(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,tin(maxi(1)+1:maxi(2)));
         y(find(t>t1&y(maxi(1)+1:maxi(2))>0)+maxi(1))=0;
        if length(maxi)>2 
            y(maxi(2)+1:maxi(3))=driftmod2(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,tin(maxi(2)+1:maxi(3)));
    %         d=y(maxi(1)+1:maxi(2))-fdrift;
    %         indup=find(d>0);
    %         y(indup+maxi(1))=fdrift(indup)+10*d(indup);


            if length(maxi)>3
                 y(maxi(3)+1:maxi(4))=driftdermod(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,tin(maxi(3)+1:maxi(4)));
            end

            if length(maxi)>4
                y(maxi(4)+1:maxi(5))=y([1:indactminclean,min(inddriftclean,length(t)-1):length(t)]);
                y(maxi(4)+1:maxi(5))=abs(y(maxi(4)+1:maxi(5))-driftmod2(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,tin(maxi(4)+1:maxi(5))));
            end

            if length(maxi)>5
%                 cuckoo=abs(driftmod2(0,zA2,0,zm2,Inf,ztau2,ztdriftswitch,0,tin(maxi(5)+1:maxi(6))));
                
                
                if y(maxi(3)+1)<0
                    d1=(driftmod2(zA1,0,zm1,0,ztau1,Inf,ztdriftswitch,0,tin(maxi(5)+1:maxi(6))));
                    [m1,i1]=min(d1);
                    r1=(abs(m1-max(d1)));
                    d1(1:i1)=m1;
                    tcurr=tin(maxi(5)+1:maxi(6));
                    d1=abs(d1-m1);
                else
                    d1=(driftmod2(zA1,0,zm1,0,ztau1,Inf,ztdriftswitch,0,tin(maxi(5)+1:maxi(6))));
                    [m1,i1]=max(d1);
                    r1=(abs(m1-min(d1)));
                    d1(1:i1)=m1;
                    tcurr=tin(maxi(5)+1:maxi(6));
                    d1=abs(d1-m1);
                    
                end
%                 d1(tcurr>=ztdriftswitch)=m1;
                

%                 y(maxi(5)+1:maxi(6))=driftdermod(0,zA2,0,zm2,Inf,ztau2,ztdriftswitch,tcurr)-driftdermod(zA1,0,zm1,0,ztau1,Inf,ztdriftswitch,tcurr);
                y(maxi(5)+1:maxi(6))=driftmod2(0,zA2,0,zm2,Inf,ztau2,ztdriftswitch,0,tcurr)-driftmod2(zA1,0,zm1,0,ztau1,Inf,ztdriftswitch,0,tcurr);

                y(maxi(5)+1:maxi(5)+find(tcurr>ztdriftswitch,1)-1)=0;
                %                 try
%                 y(maxi(5)+idwindle:maxi(6))=0;
%                 catch e
%                     disp(e);
%                 end
            end
%             if length(maxi)>6
%                 yy=y(1:maxi(1))-driftmod2(zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,tin(maxi(6)+1:maxi(7)));
%                 yy(yy>0)=0;
%                 y(maxi(6)+1:maxi(7))=abs(yy)^2;
%             end
        end
    end
end

    catch e
        getReport(e)
    end
% 
if sum(isnan(y))>0 || sum(~isfinite(y))>0
  disp('nann brad')  
end
% if ~isempty(maxi) && length(t)==length(signal)
%     figure(2)
%     if dde<t(end)-ton-tdecay
%         plot(t,y(1:maxi(1)),t,signal,ton+tdecay+dde,y(find(t>=ton+tdecay+dde,1)),'og',ton+dact,y(find(t>=ton+dact,1)),'og',ton+tdecay,y(find(t>=ton+tdecay,1)),'or');
%         hold on
%         if exist('t0')
%         plot(t0,y(find(t>=t0,1)),'xg',t2,y(find(t>=t2,1)),'xg');
%         else
%             disp('no t0')
%         end
%         
% %         plot(ton+tdecay+xde,linde,ton+tdecay+xde,Aoff*(hillde));
% %         if exist('idwindle') &&  exist('cuckoo')
% %             plot(t(idwindle),y(idwindle),'ob');
% %             plot(tin(maxi(5)+1:maxi(6)),cuckoo)
% %         end
%         hold off
%     else
%         plot(t,y(1:maxi(1)),t,signal,ton+dact,y(find(t>=ton+dact,1)),'og',ton+tdecay,y(find(t>=ton+tdecay,1)),'or');
%    
%     end
%     hold on
%     plot(ton,y(find(t>=ton,1)),'vc')
%     plot(tin(maxi(5)+1:end),y(maxi(5)+1:end));
%     drawnow;
%     hold off
% end
% %     
end

