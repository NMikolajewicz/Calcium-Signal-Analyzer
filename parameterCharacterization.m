 function out = parameterCharacterization(t,f)
global signalact tsignalact tsignal tmax noise fitting signal dfb dff plt tdrift fdrift driftmod driftmod2 tdriftswitch tend driftdermod driftdermod2 sigma indactmin inddrift inddriftclean indactminclean
global wint

%estimate instrumental/experimental drift in the first 8 seconds
%in order to characterize the noise of setup and determine if there
%a response > 4*StdDev within a reasonable time window (25% of the total recording)
ind=t<8;
%lets rezero around the first datapoint:
offset=f(1);
f=f-f(1);
signal=f;
tsignal=t;




function plotEmpty(ttl)
   figure(1);
%     set(gcf, 'Position', get(0, 'Screensize'));
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    h=plot(t,f+offset,'.b',t,fufine(:,1)+offset,'-r',t,f(spike)+offset,'xc');
    set(h(2),'linewidth',2);
    set(gca,'FontSize',18);
    xlabel('Time (s)');
    ylabel('Fluorescence (A.U.)');
    legend('Data','TV Data Estimate');
    set(gcf,'Color','w');
    ht=title(ttl);
    set(ht,'FontSize',8);
end

function cornerLegend(hleg)
    set(hleg,'fontsize',10);
    posleg=get(hleg,'Position');
    set(hleg,'Position',[0 1-posleg(4) posleg(3) posleg(4)]);
    
    set(hleg,'units','pixels')
    posleg=get(hleg,'Position');
    set(hleg,'Position',[2 posleg(2)-2 posleg(3) posleg(4)]);
    
end

function annotateActivation()
        x=x_to_norm(actfit.ton,actfit.ton);
        x(x>1)=1;
        x(x<0)=0;
        y=ylim();
        y=y_to_norm(y(2),actsol(find(t>=actfit.ton,1)));
        y(y>1)=1;
        y(y<0)=0;
        annotation('textarrow',x,y,'String',['t_{switch} =  ' num2str(actfit.ton) ' s'])  
end

function annotateResponse()
    

        x=x_to_norm(tdec,tdec);
        y=y_to_norm(rmax,r5050);
        annotation('line',x,y,'LineStyle',':','Color','green');

        %           x=x_to_norm(respfit.ton,respfit.ton);
        %         ton=respfit.ton;
        %         y=y_to_norm(0.5*fmax,respfit(ton));
        %         annotation('textarrow',x,y,'String',['t_{on} =  ' num2str(respfit.ton) ' s'])  
        %         
        x=x_to_norm(out.ton,tdec);
        y=ylim();

        y=y_to_norm(rmax,rmax);
        annotation('textarrow',x,y,'String','') 

        x=x_to_norm(out.ton,out.ton);
        y=ylim();
        y=y_to_norm(y(2),f(find(t==out.ton,1)));
        annotation('textarrow',x,y,'String',{['t_{act} =  ' num2str(out.ton) ' s'];['Amplitude = ' num2str(out.amplitude) ' A.U.']},'HorizontalAlignment','left')  

        %         
        if isfield(out,'decaytime')

            idecay0=idecay0-1;
            idecay=find(t>=decay,1);
            x=x_to_norm(decay,decay);
            y=ylim;
            y=y_to_norm(y(2),ffit(idecay-idecay0));

            annotation('textarrow',x,y,'String',{['t_{0,decay} = ' num2str(decay) ' s'] ;['  \tau_{decay} =  ' num2str(out.decaytime) ' s']},'HorizontalAlignment','right');

            mt=min(decay+out.decaytime,t(end));

            idecay2=find(t>=(mt),1);         
            x=x_to_norm(mt,mt);
            y=ylim;
            y=y_to_norm(ffit(idecay-idecay0)+0.75*(y(2)-ffit(idecay-idecay0)),ffit(idecay2-idecay0));
            annotation('textarrow',x,y,'String','');

            x=x_to_norm(decay,mt);
            y=ylim;
            y=y_to_norm(ffit(idecay-idecay0)+0.75*(y(2)-ffit(idecay-idecay0)),ffit(idecay-idecay0)+0.75*(y(2)-ffit(idecay-idecay0)));
            annotation('line',x,y);
            
        end
end



% cf(1);plot(t,f);drawnow;

tic;
try
[fu,spike,imu,fufine,noise] = TVRegSpikeRemove(f, 3, 0,0);
catch e
    disp(e)
    out=[];
    return;
end

if toc > 15
    disp('long one');
   dt=datestr(now,'mmmm_dd_yyyy_HH.MM.SS.FFF');
   save(['long_one.' dt '.mat'],'t','f');
end
    
   dt=mean(diff(t));
tend=t(find(~spike,1,'last'));

if isempty(imu)
    out=[];

        disp('No Discernible Response: Derivative did not yield and activation estimates');
    return;
end

yp = fu(:,2);

%if monotonically decreasing function or derivative
if sum(yp>0)<6 
    

    
    out=[];
    disp('No Discernible Response: Monotonically decreasing drift.');
    return;

end

% hh
% yp = TVRegDiff(f, iters, 0.001 , [], 'small', 1e-6, 0.01, 0, 0 );
% %beginning of estimate cannot be trusted
% yp(1:iters)=yp(iters);
% 
% % cf(55);plot(yp);
% % yp=yp(1:end-1);

%find the first significant upstroke
ndrift=2;
% =find(yp(ndrift:end-1)>0&(yp(ndrift:end-1)-yp(ndrift-1:end-2)>var(yp)|(yp(ndrift:end-1)>yp(ndrift-1:end-2)&yp(ndrift:end-1)==max(yp(ndrift:end-1)))),1)+ndrift;
ipmax2=find(yp(ndrift:end-1)>0&(yp(ndrift:end-1)-yp(ndrift-1:end-2)>std(yp)|(yp(ndrift:end-1)>yp(ndrift-1:end-2)&yp(ndrift:end-1)==max(yp(ndrift:end-1)))),1)+ndrift;
ipmax2=find(yp>std(yp),1);
ipmax2=imu(1);



% [ypmax,ipmax]=max(yp);


[ypmin,ipmin]=min(abs(yp(ipmax2+1:end)));
ipmin2=find(abs(yp(ipmax2+1:end))<noise/2,1);
if ~isempty(ipmin2)
    ipmin=min(ipmin,ipmin2);
end

ipmin=ipmin+ipmax2-1;


yp1=yp(1:ipmax2);
yp2=yp(ipmin:end);
%only look at positive derivatives
yp0=find(yp(1:ipmax2)<0,1,'last');
%if we crossed zero get rid of those points before the crossing
if ~isempty(yp0) && yp0<ipmax2
    yp1=yp1(yp0:end);
else
    yp0=2;
end
%shift derivative scale to baseline
yp1=yp1-min(yp1);

ypmax1=yp1(end);


% yp2=yp2-max(yp2);
yp2=abs(yp2-yp2(end));
ypmin2=min(yp2);



% indactmin=find(yp1<var(yp1),1,'last')+yp0;

indactmin=int16(find(yp(yp0:ipmax2)/yp(ipmax2)<std(yp(1:yp0)/yp(ipmax2)),1,'last'))+(yp0-1);
if isempty(indactmin)
    indactmin=int16(find(yp(yp0:ipmax2)/yp(ipmax2)<std(yp(yp0:ipmax2)/yp(ipmax2)),1,'last'))+(yp0-1);
    if isempty(indactmin)
        indactmin=2;
    end
else
    indactmin=min(indactmin,int16(find(yp(yp0:ipmax2)<var(yp(yp0:ipmax2)),1,'last'))+(yp0-1));
end


nn=3;
 inddrift=[];
while nn>=0 && isempty(inddrift);
    inddrift=find((yp2-mean(yp2))>nn*std(yp2),1,'last')+ipmin;
    nn=nn-1;
end

if isempty(inddrift)
    inddrift=find(abs(yp2-mean(yp2))>std(yp2),1,'last')+ipmin;
    inddrift=min(find(~spike,2,'last'));
end

% cf(99);
% subplot(3,1,1);
% plot(yp);
% subplot(3,1,2);
% plot(yp1);
% subplot(3,1,3);
% plot(yp2);
try

    [fmax,indmax]=max(fu(:,2));
    mxu=zeros(size(imu));
    mnu=zeros(size(imu));
    di=zeros(size(imu));
    for ii=1:length(imu)
        i0n=find(fu(1:imu(ii),2)<=0,1,'last');
        if isempty(i0n)
            i0n=1;
        end
        
        mnu(ii)=fu(i0n,1);
        i0x=find(fu(imu(ii):end,2)<=0,1)+imu(ii)-1;
        
        if ~isempty(i0x)
            mxu(ii)=fu(i0x,1);

             mnu(ii)=mnu(ii)+mean(fu([max(i0n-1,1),i0n],2))*(i0x-i0n)/(length(t)*dt);
            di(ii)=(i0x-i0n);
           
        else
            mxu(ii)=-Inf;
        end
    end
catch e
    disp(e)
end
%     [~,ii]=max((mxu-mnu));
if sum(isfinite(mxu))>2
    mxu=mxu(isfinite(mxu));
    mnu=mnu(isfinite(mxu));
    di=di(isfinite(mxu));
    ii=[];
    nn=3;
    try
    mu=(mxu-mnu)./di;
    nu=(mxu-mnu);
    mu(nu<mean(nu)-std(nu))=0;
    
    while nn<=3&&nn>0&&isempty(ii)
        ii=find(mu>mean(mu)-nn*std(mu)&mu~=0,1);
        nn=nn-1;
    end
    if isempty(ii)
        ii=length(imu);
    end
    catch e
        disp(e)
    end
else
    [~,ii]=max((mxu-mnu));    
end
    try
    indmax2=find(fu(imu(ii)+1:end,2)/length(t)<=-noise,1)+imu(ii)+1;
   
%     [fmax,indmax2]=max(fu(imu(ii):indmax2,1));
    catch e
    disp(e)
    end
    try
%     indmax=imu(ii)+indmax2-1;
indmax=indmax2;
    catch e
        getReport(e)
        disp(e)
    end
    fmax=fu(indmax,1);
    tmax=t(indmax+1);
    
    indmax0=indmax;

    fact=f(~spike(1:indmax));
    tact=t(~spike(1:indmax));
    [fmin,indmin]=min(fact);
    
%     indact2=find(fact(tact>t(indmin))-fmin>=0.1*(fmax-fmin),1)+indmin;
    
    %this only works with little noise...needs to mad more robust
%     indact=find(fact(tact>t(indmin))-fmin>=0.3*(fmax-fmin),1)+indmin;
%     
%  if ~isempty(indact)
%      indactmin=max(int16(indactmin),int16(0.8*indact));
%  end

driftind=[1:max(2,indactmin),min(inddrift,length(t)-1):length(t)];
driftind=driftind(~spike(driftind));

tdrift=t(driftind);
fdrift=f(driftind);

actdriftind=[1:indmax,min(inddrift,length(t)-1):length(t)];

  sigma=mean(abs(fu(~spike(1:indactmin),1)-f(~spike(1:indactmin))));
  

    tact2dirty=t(actdriftind);
    ypact=fufine(actdriftind,2)/(length(t)*dt);
     wyp=exp(-((f(actdriftind)-fufine(actdriftind,1)).^4)/(sigma^4));
    wyp(indactmin:indmax)=1;
    
actdriftind=actdriftind(~spike(actdriftind));

    fact2=f(actdriftind);
    tact2=t(actdriftind);



tactmin=t(indactmin);


  ypdrift=fu(driftind,2)/(length(t)*dt);


  ypdriftmin=min(min(ypdrift),0);
  ypdriftmax=max(max(ypdrift),0);
  
  
maxA=2*(max(f)-min(f));
  
lact=[-maxA,-maxA,-Inf,-Inf,0,0,0,-Inf];
uact=[maxA,maxA,Inf,Inf,t(end),t(end),t(end),Inf];


  sigma=mean(abs(fu(~spike(1:indactmin),1)-f(~spike(1:indactmin))));
  

ld=[-maxA,-maxA,ypdriftmin,ypdriftmin,0,0,0,-Inf];
ud=[maxA,maxA,ypdriftmax,ypdriftmax,t(end),t(end),t(end),Inf];

% if ypdrift(1)<0;
%     ud(1)=sigma;
%     u(1)=sigma;
% else
%     ld(1)=-sigma;
%     l(1)=-sigma;
% end



tdriftswitch=t(ipmax2);
tdriftswitch0=t(ipmax2);


inddriftclean=find(~spike(inddrift:end),1)+(inddrift-1-sum(spike(1:inddrift-1+find(~spike(inddrift:end),1))));

driftmod=@(A1,A2,m1,m2,tau1,tau2,z,x) [A1*(1-exp(-x(x<=tdriftswitch)/tau1))+m1*x(x<=tdriftswitch)+z;A1*(1-exp(-x(x>tdriftswitch)/tau1))+A2*(1-exp(-(x(x>tdriftswitch)-tdriftswitch)/tau2))+m2*(x(x>tdriftswitch)-tdriftswitch)+m1*x(x>tdriftswitch)+z];

driftmod2=@(A1,A2,m1,m2,tau1,tau2,tdriftswitch,z,x) [A1*(1-exp(-x(x<=tdriftswitch)/tau1))+m1*x(x<=tdriftswitch)+z;A1*(1-exp(-x(x>tdriftswitch)/tau1))+A2*(1-exp(-(x(x>tdriftswitch)-tdriftswitch)/tau2))+m2*(x(x>tdriftswitch)-tdriftswitch)+m1*x(x>tdriftswitch)+z];


driftdermod=@(A1,A2,m1,m2,tau1,tau2,tdriftswitch,x) [(A1/tau1)*(exp(-x(x<=tdriftswitch)/tau1))+m1;(A1/tau1)*(exp(-x(x>tdriftswitch)/tau1))+(A2/tau2)*(exp(-(x(x>tdriftswitch)-tdriftswitch)/tau2))+m2+m1];
d2mod=@(A1,A2,m1,m2,tau1,tau2,tdriftswitch,z,x)[driftmod2(A1,A2,m1,m2,tau1,tau2,tdriftswitch,z,x(1:floor(length(x)/2)));driftdermod(A1,A2,m1,m2,tau1,tau2,tdriftswitch,x(floor(length(x)/2)+1:end))];

indactminclean=indactmin-sum(spike(1:indactmin));

sact=[(f(end)-f(1)),(f(inddrift)-f(indactmin))+ypdrift(indactminclean)*((t(inddrift)-tdriftswitch))+(f(end)-f(1))*(exp(-t(inddrift)/t(indactmin))),ypdrift(indactminclean),ypdrift(end), t(indactmin),(t(inddrift)-tdriftswitch)/5,tdriftswitch0,0];

sd=[(f(end)-f(1)),(f(inddrift)-f(indactmin))+ypdrift(indactminclean)*((t(inddrift)-tdriftswitch)),ypdrift(indactminclean),ypdrift(end), t(indactmin),(t(inddrift)-tdriftswitch)/5,tdriftswitch0,0];


try 
 
    driftfit=fit(tdrift,fdrift,driftmod2,'Startpoint',sact,'Lower',lact,'Upper',uact);
    
%     
%     
%     cf(99);
%     subplot(2,1,1);
%     plot(t,f,'-b',t(spike),f(spike),'xg',t,fufine(:,1),'-k');
%     subplot(2,1,2);
%     % plot(t,f,'-b',tdrift,fdrift,'or',t(indactmin:ipmax2),f(indactmin:ipmax2),'og',t,driftfit(t),'-k',t,fufine(:,1),'-b',t,fu(:,1),'-c');
%     plot(t,f,'-b',tdrift,fdrift,'or',t(indactmin:ipmax2),f(indactmin:ipmax2),'og',t,driftfit(t),'-k');
%     ylim([min(f),max(f)]);
    
        c=differentiate(driftfit,tdrift)-ypdrift;
        d=double(ones(size(c)));

        cr=find((c(1:end-1)>0&c(2:end)<0)|(c(1:end-1)<0&c(2:end)>0),1,'last');
        
        d(1:indactminclean)=exp(-((abs(c(1:indactminclean))).^2)/(mean(abs(c(1:indactminclean))))^2);
        d(indactminclean+1:end)=exp(-((abs(c(indactminclean+1:end))).^2)/(mean(abs(c(indactminclean+1:end))))^2);


        d(indactminclean+1:end)=d(indactminclean+1:end)*sum(d(1:indactminclean))/sum(d(indactminclean+1:end));
%     dadj(indactminclean+1:end)=d(indactminclean+1:end)*sum(d(1:indactminclean))/sum(d(indactminclean+1:end));
    dadj=d;
    d2fit=fit([tdrift;tdrift],[fdrift;ypdrift],d2mod,'Startpoint',sd,'Lower',ld,'Upper',ud,'Weights',[dadj;dadj]);
catch e
    disp('not enough drift')
    disp(e)
end


df=d2fit;



  c=driftdermod(df.A1,df.A2,df.m1,df.m2,df.tau1,df.tau2,df.tdriftswitch,tdrift)-ypdrift;
        d=double(ones(size(c)));

        cr=find((c(1:end-1)>0&c(2:end)<0)|(c(1:end-1)<0&c(2:end)>0),1,'last');
        
        d(1:indactminclean)=exp(-((abs(c(1:indactminclean))).^2)/(mean(abs(c(1:indactminclean))))^2);
        d(indactminclean+1:end)=exp(-((abs(c(indactminclean+1:end))).^2)/(mean(abs(c(indactminclean+1:end))))^2);


        d(indactminclean+1:end)=d(indactminclean+1:end)*sum(d(1:indactminclean))/sum(d(indactminclean+1:end));
%     dadj(indactminclean+1:end)=d(indactminclean+1:end)*sum(d(1:indactminclean))/sum(d(indactminclean+1:end));
    dadj=d;

% if d2fit.A2>0
%     ld(2)=-sigma;
% else
%     ud(2)=sigma;
% end

% cf(77);plot(driftderfit,tdrift,ypdrift);
hold on


% cf(66);plot(t,driftmod2(df.A1,df.A2,df.m1,df.m2,df.tau1,df.tau2,df.tdriftswitch,df.z,t),t,f,tdrift,fdrift,'ob');
% ylim([min(fdrift) max(fdrift)]);
d2=d2fit([tdrift;tdrift]);



fd=d2(1:length(tdrift));
fpd=d2(length(tdrift)+1:end);

% cf(188);plot(tdrift,fd,'-r',tdrift,fdrift,'-b',tdrift,fpd,':r',tdrift,ypdrift,':b');
noise=mean(abs(fu(~spike)-f(~spike)));

driftsol=driftmod2(df.A1,df.A2,df.m1,df.m2,df.tau1,df.tau2,df.tdriftswitch,df.z,t);

if ~isempty(plt) && plt
    cf(1);

%     set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    h=plot(t(~spike),f(~spike)+offset,'.b',t(spike),f(spike)+offset,'xc',tdrift,fdrift+offset,'.r',t(indactmin:ipmax2),f(indactmin:ipmax2)+offset,'.g',t,fufine(:,1)+offset,'-b',t,driftsol+offset,'-k');
    ylim([min(f),max(f)]+offset);
    xlim([min(t),max(t)]);
    set(gca,'FontSize',18);
    xlabel('Time (s)');
    ylabel('Fluorescence (A.U.)');
    set(h,'linewidth',2);
    set(gcf,'Color','w');
    hlabel=h;
    hleg=legend(hlabel,'Data','Noise Spike','Estimated Drift','Estimated Activation','TV Data Estimate','Initial Drift Fitting','Location','NorthWest');
    
    cornerLegend(hleg);

    
end



%if the estimated "response" is not significantly different
%from noise (4*StdDev) then there is no discernable reponse
if sum(diff(find(f-driftsol>sqrt(2)*noise))==1)<6
        cf(99);

    plot(t,f+offset);
    
    out=[];
    disp('Predicted Response not Discernable from Noise.');
    return;
else
    
%there is a significant response!
    out=struct();
%first we will estimate the time of receptor activation (ton).
%for this purpose we will truncate the signal when it reaches its
%maximal value and fit the beginning of the signal with a 
%piecewise function which switches from being a simple linear
%to the sum of the same linear function and a hill function.


     minn=1;
    maxn=log(realmax)/(2*log(2*t(end)));
   maxm=abs(max(f)-min(f))/t(end);
 maxm=max(abs(fu(:,2)/(dt*length(t))));

    dff=[f(2:end)-f(1:end-1);0];
    dfb=[0;f(1:end-1)-f(2:end)];
    
% driftmod=@(A1,A2,m,tau1,tau2,tdriftswitch,z,x) [A1*(1-exp(-x(x<=tdriftswitch)/tau1))+m*x(x<=tdriftswitch)+z;A1*(1-exp(-x(x>tdriftswitch)/tau1))+A2*(1-exp(-(x(x>tdriftswitch)-ton)/tau2))+m*x(x>tdriftswitch)+z];
% 
% driftmod2=@(A1,A2,m1,m2,tau1,tau2,tdriftswitch,z,x) [A1*(1-exp(-x(x<=tdriftswitch)/tau1))+m1*x(x<=tdriftswitch)+z;A1*(1-exp(-x(x>tdriftswitch)/tau1))+A2*(1-exp(-(x(x>tdriftswitch)-tdriftswitch)/tau2))+m2*(x(x>tdriftswitch)-tdriftswitch)+m1*x(x>tdriftswitch)+z];
% 
% 
% driftdermod=@(A1,A2,m1,m2,tau1,tau2,tdriftswitch,x) [(A1/tau1)*(exp(-x(x<=tdriftswitch)/tau1))+m1;(A1/tau1)*(exp(-x(x>tdriftswitch)/tau1))+(A2/tau2)*(exp(-(x(x>tdriftswitch)-tdriftswitch)/tau2))+m2+m1];
% 

    actfittype=fittype('activationModelDrift8(A,dact,mact,nact,rho,ton,zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,x)');
    %put some bounds on the fitting based on the fit
    %for the initial drift and some data features
  mind=1/double(indmax-indactmin);
  
      [fmin2,imin2]=min(fu(indactmin:indmax,1));
     amp=(fmax-fmin2);
     amp=max(f(1:indmax)-driftfit(t(1:indmax)));
   

%   estimate some good initial guesses based on the data
    indhalf=find(fact(tact>t(indmin))-fmin>=0.5*(fmax-fmin),1)+indmin;
    ton0=(0.25*t(indactmin)+0.75*t(ipmax2));
    ton0=t(indactmin);
    
    drift1=driftmod2(driftfit.A1,0,driftfit.m1,0,driftfit.tau1,Inf,driftfit.tdriftswitch,driftfit.z,t);
 drift1p=driftdermod(driftfit.A1,0,driftfit.m1,0,driftfit.tau1,Inf,driftfit.tdriftswitch,t);

    amp=max([fu(indactmin,1);fu(indactmin:indmax,1)]-[drift1(indactmin);drift1(indactmin:indmax)]);
    if amp<0
        out=[];
        return;
    end

    id0=find(([fu(indactmin,1);fu(indactmin:indmax-1,1)]-[drift1(indactmin);drift1(indactmin:indmax-1)])<0.5*amp&(fu(indactmin:indmax,1)-drift1(indactmin:indmax))>=0.5*amp)+double(indactmin)-1;
     d0=t(id0);
     d0=mean(d0-ton0)/(tmax-ton0);
%      xconv=find(yp(indactmin:ipmax2)-yp((indactmin:ipmax2)-1)<0,1)+indactmin-1;
%      n0=-0.5/log(((t(xconv)-t(indactmin))/(tmax-ton0))/d0);
    im0=find(fu(1:indmax,2)>0,1,'last');
    
     m0=yp(im0)/(length(t)*dt);
%      fphalf=fu(id0,2)/(length(t)*dt);
     fphalf=max(fu(indactmin:indmax,2)/(length(t)*dt));
        n0=4*(fphalf-drift1p(id0)-m0/2)*d0*(tmax-ton0);
        try
        if n0>maxn/2||n0<0
            disp(n0)
        end
        catch e
            getReport(e);
            disp(e)
        end
        xact=(ton0:dt:tmax)-ton0;
         dact0=d0*(xact(end));
        hillact=(xact.^n0)./((xact.^n0)+((dact0)^n0));
        
        hillpup=@(x,n,beta) (n*(beta^n)*(x^(n-1)))/(((x^n)+(beta^n))^2);
        m0=yp(im0)/(length(t)*dt)-amp*hillpup(xact(end-1),n0,d0*(tmax-ton0))-drift1p(im0);
        
        m0=max(0,m0);
        
        n0=4*(fphalf-drift1p(id0)-m0/2)*d0*(tmax-ton0);
        try
        if n0>maxn/2||n0<0
            disp(n0)
        end
        catch e
            getReport(e);
            disp(e)
        end
        
        a0=amp-m0*trapz(hillact,xact);
        
     lact=[max(a0/2,amp/2),mind,0,0,0,0,ld];
    uact=[Inf,1-mind,maxm*10,maxn,1,t(ipmax2),ud];
    
    A0=amp;
    
    
   
    iact=find(t>=ton0,1);
    xact=t(iact:indmax)-t(iact);
    hillact=(xact.^n0)./((xact.^n0)+((dact0)^n0));
    
    idact0=find(t>ton0+d0,1)-iact;
    idact0=max(idact0,1);
  
    
     Hillpact0=m0*hillact+n0*A0*(xact.^(n0-1))*(d0^n0)./(((xact.^n0)+d0^n0).^2);
   
    try
      uu=Hillpact0(idact0:end)-fu(iact+idact0-1:indmax,2)/(length(t)*dt)+drift1p(iact+idact0-1:indmax);
    catch e
        getReport(e)
        disp(e)
    end
    uu=abs(uu)+abs([uu(1);uu(1:end-1)])+abs([uu(2);uu(2:end)]);
    uu=uu/3;
%     cf(); plot(exp(-uu));
    [~,muu]=min(uu);
    
     rho0=1-(muu+idact0-1)/((indmax-iact+1));
     
     
    
           sact=[amp,d0,m0,n0,rho0,ton0,df.A1,df.A2,df.m1,df.m2,df.tau1,df.tau2,tmax+(t(end)-tmax)/double(2),df.z];
 
    snr=((max(f)-min(f))/noise);
        
        
    %fit the data
    if length(tact)>=8
    fitting=1;
    tsignalact=tact;
    signalact=fact;
    try
        w=ones(size(tact2));
        w(1:indactminclean)=d(1:indactminclean);
        w(1:indactminclean)=w(1:indactminclean)./(sum(w(1:indactminclean))/double(indactminclean));
        w(length(tact2)-(length(d)-indactminclean-1):end)=d(indactminclean+1:end);
        w(length(tact2)-(length(d)-indactminclean-1):end)=w(length(tact2)-(length(d)-indactminclean-1):end)./(sum(w(length(tact2)-(length(d)-indactminclean-1):end))/double((length(d)-indactminclean-1)));

        dadj=sum(w)*dadj/sum(dadj);

        wint=abs(fmin2+amp-fu(end,1))/abs(fmin2-fu(end,1));

        actfit=fit([tact2;tact2dirty;tdrift;tdrift;tact2],[fact2;ypact;fdrift;ypdrift;zeros(size(tact2))],actfittype,'Lower',lact,'Upper',uact,'StartPoint',sact,'Weight',[snr*snr*w;snr*wyp;dadj;dadj;wint*mean(w)*ones(size(tact2))],'Robust','Bisquare','Display','off','Algorithm','Trust-Region');


    catch e
        disp('error in activation fit')
        disp(e)
    end
    
%     if ~isempty(regexp(lastwarn,'starting'))
%         disp(lastwarn)
%     end
    
    
    signal=f(~spike);
    fitting=0;
    %amplitube of the response
%     out.amplitude=actfit(tact(end))-actfit(actfit.ton);
    
    %time of onset of response
%     indon= find(t-actfit.ton<0,1,'last');
    indon=indactmin;
    sigma=mean(abs(fu(~spike,1)-f(~spike)));
%     if fu(indon+1)-fu(indon)>4*sigma
%         out.ton=indon;
%     else
        actsol=actfit(t);
        out.amplitude=actsol(find(t>=tact(end),1))-actsol(find(t<=actfit.ton,1,'last'));
        
        

    if ~isempty(plt) && plt

        figure(1);
        %         subplot(2,1,1);
        hold on
        h=plot(t,actsol+offset,'LineStyle','-','Color',[1 0.5 0]);
        hlabel=[hlabel; h];
        set(h,'linewidth',2)
        hleg=legend(hlabel,'Data','Noise Spike','Estimated Drift','Estimated Activation','TV Data Estimate','Initial Drift Fitting','Activation Fitting','Location','NorthWest');
        
        hh=findall(gca,'DisplayName','Initial Drift Fitting');
        if ~isempty(hh)
            set(hh,'LineStyle',':')
        end
        
        cornerLegend(hleg);
        

        
        ylim([min(f),max(f)]+offset);

        annotateActivation();

        hold off
    end

        
        if abs(actsol(end)-fu(end,1))>0.2*amp
            disp('out in the boonies')
        end
%     if plt
%         figure(1);
%         clf(1);
%         hold on
%         plot(actfit,tact,fact);
%         x=x_to_norm(actfit.ton,actfit.ton);
%         ton=actfit.ton;
%         y=y_to_norm(0.5*fmax,actfit(ton));
%         annotation('textarrow',x,y,'String',['t_{on} =  ' num2str(actfit.ton) ' s'])  
%         hold off
%     end
    fitting=1;
    respfittype=fittype('responseModelDrift11(Aoff,Aon,dact,dde,mact,nact,nde,rhoact,rhode,tdecay,ton,zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,x)');
    
%     driftmod2=@(ton,A1,A2,m1,m2,tau1,tau2,z,x) [A1*(1-exp(-x(x<=ton)/tau1))+m1*x(x<=ton)+z;A1*(1-exp(-x(x>ton)/tau1))+A2*(1-exp(-(x(x>ton)-ton)/tau2))+m2*(x(x>ton)-ton)+m1*x(x>ton)+z];
    
    
    fitting=0;
    %find the indices where the fluorescence crosses the half max
    

    actdrift1=driftmod2(actfit.zA1,0,actfit.zm1,0,actfit.ztau1,Inf,actfit.ztdriftswitch,actfit.zz,t);
    actdrift1p=driftdermod(actfit.zA1,0,actfit.zm1,0,actfit.ztau1,Inf,actfit.ztdriftswitch,t);
    actdrift10=actdrift1;
    
    
   if sum(abs(actdrift1-fu(:,1)))<sum(abs(drift1-fu(:,1))) || sum(abs(actdrift1-fu(:,1)))<sum(abs(driftsol-fu(:,1)))
       
       if sum(abs(drift1-fu(:,1)))<sum(abs(driftsol-fu(:,1)))
       
       actdrift1=drift1;
       actdrift1p=drift1p;
   
        zA10=driftfit.A1;
        zm10=driftfit.m1;
        ztau10=driftfit.tau1;
        ztdriftswitch0=driftfit.tdriftswitch;
        zz0=driftfit.z;
       else
        actdrift1=driftmod2(df.A1,0,df.m1,0,df.tau1,Inf,df.tdriftswitch,df.z,t);
        actdrift1p=driftdermod(df.A1,0,df.m1,0,df.tau1,Inf,df.tdriftswitch,t);
   
        zA10=df.A1;
        zm10=df.m1;
        ztau10=df.tau1;
        ztdriftswitch0=df.tdriftswitch;
        zz0=df.z;
           
       end
   else
    
        zA10=actfit.zA1;
        zm10=actfit.zm1;
        ztau10=actfit.ztau1;
        ztdriftswitch0=actfit.ztdriftswitch;
        zz0=actfit.zz;
   end
   
   
      mact0=actfit.mact;
    nact0=actfit.nact;
    Aon0=actfit.A;
    
%      Hillpact0=mact0*hillact+nact0*Aon0*(xact.^(nact0-1))*(dact0^nact0)./(((xact.^nact0)+dact0^nact0).^2);
%    
%     try
%       uu=Hillpact0(idact0:end)-fu(iact+idact0-1:indmax,2)/(length(t)*dt)+actdrift1p(iact+idact0-1:indmax);
%     catch e
%         getReport(e)
%     end
%     uu=abs(uu)+abs([uu(1);uu(1:end-1)])+abs([uu(2);uu(2:end)]);
%     uu=uu/3;
% %     cf(); plot(exp(-uu));
%     [~,muu]=min(uu);
     
        dact0=d0*(xact(end));
   tact=(actfit.ton:dt:tmax);
   if length(tact)==1
    tact=[actfit.ton tmax];
   end
   xact=tact-actfit.ton;
    dact0=d0*(xact(end));
    hillact=(xact.^n0)./((xact.^n0)+((dact0)^n0));
    
     rhoact0=actfit.rho;
     
     try
        A0=Aon0*hillact(end)+mact0*trapz(xact,hillact);
     catch e
      getReport(e)
      disp(e)
     end
    A00=A0;

    
    fsmooth=fu(:,1);
    
    
    [fmin2,imin2]=min(fu(indactmin:indmax,1));
    imin2=imin2+indactmin-1;
    amp=(fmax-fmin2);
    
    endcross=find(fsmooth(indmax:end)<=actdrift1(indmax:end),1)+indmax-1;
    subdrift1=true;
    if isempty(endcross)||endcross==indmax
%         endcross=length(t);
        [~,endcross]=min(fsmooth(indmax:end));
        endcross=endcross+indmax-1;
        dsmooth=fsmooth(indmax:endcross);
        
        delta=fsmooth-min(dsmooth);
        ampcross=dsmooth(1);max(delta(indmax:endcross))-min(delta(indmax:endcross));
         subdrift1=false;
    else
        dsmooth=fsmooth(indmax:endcross)-actdrift1(indmax:endcross);
        ampcross=dsmooth(1);
        delta=fsmooth-actdrift1-min(dsmooth);
    end
    
    
    
    
    
     indcrosshalf=[delta(indmax); delta(indmax:end-1)]>A0/2&delta(indmax:end)<=A0/2;
 indcrosshalf=find(indcrosshalf,1)+indmax-1;
    thalf=mean(t(indcrosshalf));
%     
%     indcrossquarter=[f(1); f(1:end-1)]>0.75*fmax&[t(1); t(1:end-1)]>tmax&f<=0.75*fmax&t>tmax&[~spike(1);~spike(1:end-1)];
%     indcrossquarter(spike)=0;
%     tquarter=mean(t(indcrossquarter));
    if isnan(thalf)
        disp('Could not detect deactivation.')
        return;
    end
    


    [~, idx] = min(abs(t-actfit.ton)); %index of closest value
    
    hillup=@(x,n,beta) (x^n)/((x^n)+(beta^n));
    hilldown=@(x,n,beta) (beta^n)/((x^n)+(beta^n));
    hillpdown=@(x,n,beta) -(n*(beta^n)*(x^(n-1)))/(((x^n)+(beta^n))^2);
    
    tdact=actfit.dact*(tmax-actfit.ton)+actfit.ton;
    rhoact0=max(3/(indmax-find(t<=tdact,1,'last')),0.5);
%     rhode0=max(-3/(indmax-find(t<=thalf,1,'last')),0.5);
    
    
    
    
    hillend=hillup((tmax-actfit.ton),actfit.nact,actfit.dact*(tmax-actfit.ton));
    mde0=max(yp(find(t==tmax)+1)/(length(t)*dt),-maxm);
%     mde0=0;
    
    dde0=(thalf-tmax)/(tend-tmax);
   
    
    [fmin3,imin3]=min(fu(indmax:end,1));
    
    imin3=imin3+indmax-1;
    ddrift=actdrift1(imin3)-actdrift1(imin2);
    Aoff0=(actfit.A*hillend+actfit.mact*(tmax-tdact)-(fmin3-fmin2)+mde0*(thalf-tmax))+ddrift;
    

%     hillend2=Aoff0*hillup(t(end)-tmax,nde0,dde0*(thalf-tmax));
%     hillpend2=Aoff0*hillpdown((t(end)-tmax),nde0,dde0*(thalf-tmax));
    
   
    % this is really excessive we need to take into accoubt drift 1 and
    % the exponential component of drift 2 otherwise the initial geuss is 
    % overcompensating
    ztau20=t(end)-tmax;
    
    

%    
%     iact=find(t>=actfit.ton,1);
%     xact=t(iact:indmax)-t(iact);
%     nact=actfit.nact;
%     hillact=(xact.^nact)./((xact.^nact)+((actfit.dact*(xact(end)))^nact));
%     dact0=actfit.dact*(xact(end));
%     idact0=find(t>actfit.ton+dact0,1)-iact;
%     idact0=max(idact0,1);
 
    
    nde0=4;
    
%        indcrosshalf=find(indcrosshalf,1)+indmax-1;
    nn=1;
    beta0=dde0*(tend-tmax);
    
    xde=t(indmax:end)-t(indmax);

    hillde=1-(xde.^nde0)./((xde.^nde0)+(beta0^nde0));
    
    
    intlinde=trapz(xde,hillde);
    %ensures linear portion is less than half
    mde0=max(-A0/(2*intlinde),mde0);
    
   [fphalf,iphalf]=min(fu(indmax:endcross,2)/(length(t)*dt));
     fpphalf=(fu(indmax:end,2)-[fu(indmax,2); fu(indmax:end-1,2)])/((length(t)*dt^2));
     fpphalf=fpphalf(indcrosshalf-indmax+1);
%    fphalf=min(fphalf,-A/4);
    Aoff0=max(0,(A0+intlinde*mde0)/((1-hillde(end))*A0));
    
    nmin0=sqrt(2);
        
    bd=(thalf-tmax);
    
    
    while nn<6

    
     if subdrift1
        nde0=-(fphalf-mde0/2-actdrift1p(indcrosshalf))*4*bd/(Aoff0*A0);
        nde2=(fpphalf-(actdrift1p(indcrosshalf)-actdrift1p(indcrosshalf-1))/dt)/(mde0/(4*bd)+Aoff0*A0/(4*bd^2));
     else
         nde0=-(fphalf-mde0/2)*4*(thalf-tmax)/(Aoff0*A0);
          nde2=(fpphalf)/(mde0/(4*bd)+Aoff0*A0/(4*bd^2));
     end
     nde2=max(nde2,1.01);
     nde2=min(nde2,maxn/sqrt(2));
     
     
     nde0=min(nde0,maxn/sqrt(2));
     nde0=max(nde0,nmin0);
     
     nde0=min(nde2,nde0);
   
    hillde=1-(xde.^nde0)./((xde.^nde0)+(beta0^nde0));    
    intlinde=trapz(xde,hillde);
    mde0=max(-A0/(2*intlinde),mde0);
    
    Aoff0=max(0,(A0+intlinde*mde0)/((1-hillde(end))*A0));
%     Aoff0=min(Aoff0,1);
     intlindehalf=trapz(xde(xde<beta0),hillde(xde<beta0));
     Delta=max(A0-(fu(indcrosshalf,1)-actdrift1(indcrosshalf)),0);
    mde0=-(Delta-Aoff0*A0*(hillde(find(xde>=beta0,1))))/intlindehalf;
    mde0=min(mde0,0);
%      A0=((fu(indcrosshalf,1)-actdrift1(indcrosshalf))-Aoff0*A0*hillde(find(xde>=beta0,1))-mde0*intlindehalf)
%         r0=r0/A0;
    nn=nn+1;
    end
    
mde0=-A0*(1-Aoff0*(1-hillde(end)))/intlinde;

   
    tf=t(end)-ztdriftswitch0;
   
    c=A0-Aoff0*A0;
    
    
     xde=t(indmax:end)-t(indmax);
         beta0=(dde0*(t(end)-tmax));
    hillde=(1-(xde.^nde0)./((xde.^nde0)+(beta0.^nde0)));

    hillendde=hillde(end);
    
    
    delta2=Aoff0*A0*hillendde...
        +mde0*trapz(xde,1-((xde.^nde0)./((xde.^nde0)+(beta0^nde0))))...
        +c;
    
    df=mean(fu(end-1:end,1))-(actdrift1(end)+delta2);
    
    driftfinal=df;
    
    
    jj=find(d(1:end-1)>[d(1);d(1:end-2)],1,'last');
    
    driftpfinal=sum(d(jj:end).*fu(end-length(d)+jj:end,2)/(length(t)*dt))/sum(d(jj:end));
    
    
    driftpfinal=mean(fu(inddrift:end,2)/(length(t)*dt));
    
    
  
    
    Hillpde0=mde0*hillde-nde0*Aoff0*A0*(xde.^(nde0-1))*(beta0^nde0)./(((xde.^nde0)+beta0^nde0).^2);
    
    driftpfinal=driftpfinal-actdrift1p(end)-Hillpde0(end);
    %%also need to subtract drift contributions to derivative
    uu=Hillpde0(1:indcrosshalf-indmax+1)-fu(indmax:indcrosshalf,2)/(length(t)*dt)+actdrift1p(indmax:indcrosshalf);
    
    uu=abs(uu)+abs([uu(1);uu(1:end-1)])+abs([uu(2);uu(2:end)]);
    uu=uu/3;
%     cf(); plot(exp(-uu));
    [~,muu]=min(uu);
    rhode0=muu/(find(t==tend,1)-indmax+1);
    
     

     jj=max(find(t>ztdriftswitch0,1));
      if isempty(jj)|| jj>=indmax
        
      else
        jj=length(f)-1;
        
      end
   
        [~,idriftpswitch]=min(abs(fu(indmax:jj+1,2)/(length(t)*dt)));
          idriftpswitch=idriftpswitch+indmax-1   ;
         driftpswitch0=fu(idriftpswitch,2)/(length(t)*dt);
         
    jj=max(find(t>ztdriftswitch0,1));
    
    if jj>=indmax
        xswitch=xde(jj-indmax+1);
        driftpswitch=driftpswitch0-actdrift1p(jj)-mde0*hillde(jj-indmax+1)+nde0*Aoff0*A0*(xswitch^(nde0-1))*(beta0^nde0)/((xswitch^nde0)+beta0^nde0)^2;
    elseif t(jj)>actfit.ton
        xswitch=xact(find(tact>=ztdriftswitch0,1));
        try
         driftpswitch=driftpswitch0-actdrift1p(jj)-mact0*hillact(jj-iact+1)-nact0*Aon0*(xswitch^(nact0-1))*(dact0^nde0)/((xswitch^nde0)+dact0^nde0)^2;
        catch ex
            getReport(ex)
            disp(e)
        end
     else
         driftpswitch=driftpswitch0-actdrift1p(jj);

    end
    
    
    zA20=(driftfinal-tf*driftpfinal)/(1-exp(-tf/ztau20)-tf*exp(-tf/ztau20)/ztau20); 
    zA21=(driftfinal-tf*driftpswitch)/(1-exp(-tf/ztau20)-tf/ztau20); 
    

    if abs(zA21)<abs(zA20)
        zA20=zA21;
        zm20=driftpswitch-zA20/ztau20;
    else
        zm20=driftpfinal-zA20*exp(-tf/ztau20)/ztau20;
    end
    
%       zm20=driftpfinal-zA20*exp(-tf/ztau20)/ztau20;

    deltaF=max(f(~spike))-min(f(~spike));
    
%     if zA20>deltaF
%         zA20=deltaF;
%     end
    

%     if abs(zm20)>maxm
%         zm20=sign(zm20)*maxm/sqrt(2);
%         zA20=(driftfinal-zm20*tf)/(1-exp(-tf/ztau20));
%     end
   
    
   ld(1:2)=-max(abs(zA20)*sqrt(2),maxA);
   ud(1:2)=max(abs(zA20)*sqrt(2),maxA);
   
%    if actfit.zA1>0
%         ld(1)=0;
%    else
%        ud(1)=0;
%    end

 
   
   ld(3:4)=-max(abs(zm20)*sqrt(2),maxm);
   ud(3:4)=max(abs(zm20)*sqrt(2),maxm);
   
    lresp=[0,0,0,0,0,n0/2,1.01,0,0,0,0,ld];
    uresp=[Inf,uact(1),1-mind,1-mind,maxm,maxn,maxn,1,2,t(end),t(ipmax2),ud];
   
%     Aoff0=min(Aoff0,1)
% 
    sresp=[Aoff0,actfit.A,actfit.dact,dde0,actfit.mact,actfit.nact,nde0,rhoact0,rhode0,tmax-actfit.ton,actfit.ton,zA10,zA20,zm10,zm20,ztau10,ztau20,ztdriftswitch0,zz0];
%      sresp=[Aoff0,actfit.A,actfit.dact,dde0,actfit.mact,mde0,actfit.nact,nde0,rhact0,rhode0,tmax-actfit.ton,actfit.ton,actfit.zA1,0,actfit.zm1,0,actfit.ztau1,t(end),actfit.ztdriftswitch,actfit.zz];
   
%     sresp=[Aoff0,actfit.A,actfit.dact,dde0,actfit.mact,mde0,actfit.nact,4,rhact0,rhode0,tmax-actfit.ton,actfit.ton,actfit.zA1,actfit.zA2,actfit.zm1,actfit.zm2,actfit.ztau1,actfit.ztau2,actfit.ztdriftswitch,actfit.zz];
    fitting=1;
    try 

        
        
        
    w=ones(size(f(~spike)));
    w(1:indactminclean)=d(1:indactminclean);
    w(1:indactminclean)=w(1:indactminclean)./(sum(w(1:indactminclean))/double(indactminclean));
    w(end-(length(d)-indactminclean-1):end)=d(indactminclean+1:end);
    w(end-(length(d)-indactminclean-1):end)=w(length(w)-(length(d)-indactminclean-1):end)./(sum(w(end-(length(d)-indactminclean-1):end))/double((length(d)-indactminclean)));
    
    
    wyp=w;
        

     dadj=sum(w)*d/sum(d);

    xresp=[t(~spike);t(~spike);tdrift;tdrift;tdrift;t];
    yresp=[f(~spike);fufine(~spike,2)/(length(t)*dt);fdrift;ypdrift;zeros(size(tdrift));zeros(size(t))];
   
    
    indmax=find(t>=tmax,1);
    gg=snr^-2*max(abs(yp(min(inddrift+1,length(yp)):end)))/max(yp(indactmin:indmax));
    try
    wresp=[(snr^2)*w;(snr)*wyp;gg*dadj;gg*dadj;snr*wint*dadj;snr*(wint)*mean(w)*ones(size(t))];
    catch e
        disp(e)
    end
    if abs(actfit.zz)>sigma/2 || abs(actfit.zA1)<3*sigma
       
        if actfit.zz>0
            ld(end)=0;
            sresp(end)=max(sresp(end),0);
        else
            ud(end)=0;
            sresp(end)=min(sresp(end),0);
        end
        
        lresp=[0,0,0,0,0,0,1.01,0,0,0,0,ld];
        uresp=[Inf,max(f)-min(f),1-mind,1-mind,maxm,maxn,maxn,1,2,t(end),t(ipmax2),ud];
        
        [respfit,gofresp]=fit(xresp,yresp,respfittype,'Lower',lresp,'Upper',uresp,'StartPoint',sresp,'Weight',wresp,'Display','off','Robust','Bisquare','Algorithm','Trust-Region');
    else
        respfittype=fittype('responseModelDrift11(Aoff,Aon,dact,dde,mact,nact,nde,rhoact,rhode,tdecay,ton,zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz,x)','problem','zz');
        [respfit,gofresp]=fit(xresp,yresp,respfittype,'Lower',lresp(1:end-1),'Upper',uresp(1:end-1),'StartPoint',sresp(1:end-1),'Weight',wresp,'problem',0,'Robust','Bisquare','Display','off','Algorithm','Trust-Region');
%             [respfit,gofresp]=fit(xresp,yresp,respfittype,'Lower',lresp(1:end-1),'Upper',uresp(1:end-1),'StartPoint',sresp(1:end-1),'Weight',wresp,'problem',0,'Display','off');
    end
    
    catch e
        disp('problem with resp fit')
        disp(getReport(e));
    end
%     disp(lastwarn)

% disp(lastwarn)
% regexp(lastwarn,'starting')
%     if ~isempty(regexp(lastwarn,'starting'))
%         disp(lastwarn)
%     end
    
    fitting=0;
    tclean=t(~spike);
    fitting=1;
    r2=respfit([tclean;tclean;tdrift;tdrift;tdrift;tclean]);
    fitting=0;
    
    fr2=r2(1:length(tclean));
    fpr2=r2(length(tclean)+1:2*length(tclean));
   dr2=r2(2*length(tclean)+2:2*length(tclean)+1+length(tdrift));
    fpdr2=r2(2*length(t)+length(tdrift)+1:end);
    
%     cf(44);plot(t,fr2,'-r',t(~spike),f(~spike),'-b',tdrift,dr2,'-r',tdrift,fdrift,'ob',tdrift,fpdr2,'-r',tdrift,ypdrift,'.b');
%     cf(44);plot(tclean,fr2,'-r',t(~spike),f(~spike),'-b',tdrift,dr2,'-r',tdrift,fdrift,'ob',tclean,fpr2,'-r',t,fufine(:,2)/(length(t)*dt),'.b');
%      cf(144);plot(tclean,fpr2,'-r',t,fufine(:,2)/(length(t)*dt),'.b',t(~spike),wyp);
%    
    try
    out.ton=find(fufine(indactmin:length(tact),1)-fufine(indactmin,1)>=sigma,1)+(indactmin-1);
    catch e
        disp(e)
        out.ton=indactmin;
    end

%     end

    
%     warning ('off','all');
%         zm=respfit.zm;
%         respfit.zm=0;
%         zA=respfit.zA;
%         respfit.zA=0;
%         zz=respfit.zz;
%         respfit.zz=0;
%      warning ('on','all');

    respest=respfit(t);
%     

    indon=out.ton;
    [minfit, iminfit]=min(respest(indon:end));
    iminfit=iminfit+(indon-1);
    %if we have something monotonically increasing or decreasing
    if iminfit==length(t) || iminfit==indon
        intend=length(t);
    else% otherwise the minimum in the signal is local
        intend=iminfit;
    end
%     
 
    
    tfine=linspace(t(1),t(end),length(t)*10)';
    

    
    
    
        
    if isempty(out.ton)
        out.ton=t(indactmin);
    else
        max(out.ton,find(t<=actfit.ton,1,'last'));
        try
        while respest(out.ton+1)<respest(out.ton) && out.ton<length(t)
            out.ton=out.ton+1;
        end
        catch e
            disp(getReport(e))
        disp('poop')    
        end
        out.ton=t(max(out.ton,find(t<=actfit.ton,1,'last')));
    end
    
    
    respdrift=driftmod2(respfit.zA1,respfit.zA2,respfit.zm1,respfit.zm2,respfit.ztau1,respfit.ztau2,respfit.ztdriftswitch,respfit.zz,t);
       
    
    zA1=respfit.zA1;
    zA2=respfit.zA2;
    zm1=respfit.zm1;
    zm2=respfit.zm2;
    ztau1=respfit.ztau1;
    ztau2=respfit.ztau2;
    ztdriftswitch=respfit.ztdriftswitch;
    zz=respfit.zz;
    
    
    warning ('off','all');
    
    respfit.zA1=0;
    respfit.zA2=0;
    respfit.zm1=0;
    respfit.zm2=0;
    respfit.ztau1=Inf;
    respfit.ztau2=Inf;
    respfit.ztdriftswitch=0;
    respfit.zz=0;
    
    warning ('on','all');
    
    tend=t(end);
    
    detrend=respfit(t);
    
    [respmax,indmax]=max(detrend);
    tmax=t(indmax);
    
    idwindle=find((detrend(indmax:end)-detrend(end))<sigma/4,1)+indmax-1;
    if ~isempty(idwindle)
        intend=idwindle;
        tdwindle=t(idwindle);
        
    end
    
    detrendfine=respfit(tfine);
    
    t10=find(tfine>respfit.ton&tfine<tmax&detrendfine>=0.1*respmax,1);
    t50up=find(tfine>respfit.ton&tfine<tmax&detrendfine>=0.5*respmax,1);
    t50down=find(tfine>tmax&(detrendfine-detrend(end))<=0.5*(respmax-detrend(end)),1);
    t90=find(tfine>respfit.ton&detrendfine>=0.9*detrend(indmax),1);
    snr=((max(f)-min(f))/noise);
    
%     tint=find();
    ipass=find(respest<respdrift,1);
    isig=find(detrend(indmax)<respmax/snr,1)+indmax-1; 
    if ~isempty(isig)
        if ~isempty(ipass)
            ilow=max(ipass,isig);
        else
            ilow=isig;
        end
    else
        ilow=ipass;
    end

    if ~isempty(ilow)
       intend=min(intend,ilow);
    end
    

    ion=find(t>=respfit.ton,1);
    out.AUC=trapz(t(ion:intend),detrend(ion:intend));
    
    
    respfit.zA1=zA1;
    respfit.zA2=zA2;
    respfit.zm1=zm1;
    respfit.zm2=zm2;
    respfit.ztau1=ztau1;
    respfit.ztau2=ztau2;
    respfit.ztdriftswitch=ztdriftswitch;
    respfit.zz=zz;
    
    %this is abuse of the variable name
        respfine=respfit(tfine);
    
    out.tact=tfine(t90)-tfine(t10);
    
    if ~isempty(t50down) && ~isempty(t50up)
        out.FWHM=tfine(t50down)-tfine(t50up);
    else
        return;
    end
    
    
    

    beta=respfit.dde*(tend-respfit.ton-respfit.tdecay);
    n=respfit.nde;
%     if n>1
%         tconc = beta*((n-1)/(n+1))^(1/n);
%     else
        tconc=0;
%     end

    ffit=respfit(t);
    
    idecay0=find(t>=(respfit.ton+respfit.tdecay+tconc),1);
    ffit=ffit(idecay0:end);
    
    iconv = find(((detrend(idecay0+1:end-1)-detrend(idecay0:end-2))-(detrend(idecay0+2:end)-detrend(idecay0+1:end-1)))<0,1)-1+idecay0;
    if ~isempty(iconv)
        decay = t(iconv);


        %needs to be able to handle the case where dde is beyond the timeframe
        %of the experiment...also that just shouldnt happen
        fdecay=f(t>decay);
        tdecay=t(t>decay)-decay;

       Amax=2*(1+2*sigma/((max(f)-min(f))))*(max(f)-min(f));

        out=characterizeDecay(tdecay,fdecay,decay,Amax,out);
    end
      tactivation=out.ton;
    
      
        i50up=find(t>=tfine(t50up),1);
        i50down=find(t>=tfine(t50down),1);
        resp5050=respfine(t50up:t50down);
        mr5050=(respfine(t50down)-respfine(t50up))/(tfine(t50down)-tfine(t50up));
        br5050=respfine(t50up)-mr5050*tfine(t50up);
        lr5050=mr5050*t(i50up:i50down)+br5050;
    
    if plt
        figure(1);
        
        delete(findall(gcf,'type','textarrow'));
        xlabel('')
        pos=get(gca,'Position');
        set(gca,'Position',[pos(1) pos(2)+pos(4)/2 pos(3) pos(4)/2])
        set(gca,'fontsize',12);
%         annotateActivation();
        set(gca,'XTickLabel','')
        
        hold on
        subplot(2,1,2)
        set(gca,'fontsize',12);

% cf(2)
        

        patch([t(ion:intend);t(intend:-1:ion)],[respest(ion:intend); respdrift(intend:-1:ion)]+offset,'red','FaceColor',[1 0.8 0.8],'LineStyle',':');
             hold on
        h=plot(t(~spike),f(~spike)+offset,'.b');% if ~isempty(maxi) && length(t)==length(signal)

        xlim([0,t(end)]);
        
        ylimf=ylim;
        hold on

        
        d5050=respdrift(i50up:i50down);
        md5050=(respdrift(i50down)-respdrift(i50up))/(t(i50down)-t(i50up));
        bd5050=respdrift(i50up)-md5050*t(i50up);
        ld5050=md5050*t(i50up:i50down)+bd5050;
        
        l5050=lr5050-(ld5050-respdrift(i50up:i50down));
    
        warning ('off','all');
    
        h=plot(t,respfit(t)+offset,'-r',t(i50up:i50down),l5050+offset,'-g',t,driftmod2(respfit.zA1,respfit.zA2,respfit.zm1,respfit.zm2,respfit.ztau1,respfit.ztau2,respfit.ztdriftswitch,respfit.zz,t)+offset);
        set(h(end),'Color','k','LineStyle','-.')
        hlabel=[hlabel; h(1); h(end)];
        hleg=legend(hlabel,'Data','Estimated Drift','Estimated Activation','TV Data Estimate','Initial Drift Fitting','Activation Fitting','Response Fitting','Final Baseline Estimate','Location','NorthWest');
        
        warning ('on','all');
        
        hh=findall(gca,'DisplayName','Initial Drift Fitting');
        if ~isempty(hh)
            set(hh,'LineStyle',':')
        end
        
        cornerLegend(hleg)
        
        
        ylim(ylimf);
        set(h,'linewidth',2)
        set(gca,'Box','on');
        
        
        t5050=mean(tfine([t50up,t50down]));
        i5050=find(t>=t5050,1);
        
        r5050=l5050(i5050-i50up+1);
        x=x_to_norm(t5050,t5050);
        y=ylim;
        y=y_to_norm(y(1)+offset,r5050+offset);

        tdec=tmax;
        rmax=detrend(indmax)+respdrift(indmax);

        
        h=annotation('textbox','String',{['FWHM = '  num2str(out.FWHM) 's'];['t_{decline} = '  num2str(tdec) 's'];['AUC = '  num2str(out.AUC) ' A.U. \cdot s']},'Color','black','HorizontalAlignment','left','BackgroundColor','w','LineStyle','none');drawnow;
        
        pos=get(h,'position');
        set(h,'position',[x(1)-pos(3)/2  0 pos(3) pos(4)])
        

        
        pos=get(h,'position');
        set(h,'position',[x(1)-pos(3)/2  0 pos(3) pos(4)])
        annotation('line',x,[pos(4),y(2)],'Color','green','LineWidth',2);
        
        
        x=x_to_norm(tdec,tdec);
        y=y_to_norm(rmax+offset,l5050(indmax-i50up+1)+offset);
        annotation('line',x,y,'LineStyle',':','Color','green','LineWidth',2);

        %           x=x_to_norm(respfit.ton,respfit.ton);
        %         ton=respfit.ton;
        %         y=y_to_norm(0.5*fmax,respfit(ton));
        %         annotation('textarrow',x,y,'String',['t_{on} =  ' num2str(respfit.ton) ' s'])  
        %         


        x=x_to_norm(0,out.ton);
        y=ylim();
        y=y_to_norm(mean(y),f(find(t>=out.ton,1))+offset);
%         annotation('textarrow',x,y,'String',{['t_{onset} =  ' num2str(out.ton) ' s'];['Amplitude = ' num2str(out.amplitude) ' A.U.']},'HorizontalAlignment','left','TextBackgroundColor','w','LineWidth',2)  

         
              
        h=annotation('textbox','String',{['t_{onset} =  ' num2str(out.ton) ' s'];['Amplitude = ' num2str(out.amplitude) ' A.U.'];['t_{10%-90%} =  ' num2str(out.tact) ' s'];},'Color','black','HorizontalAlignment','right','BackgroundColor','w','LineStyle','none');drawnow;
          
        pos=get(h,'position');
        set(h,'position',[0 y(1)-pos(4)/2 pos(3) pos(4)])

%         y(1)=y(1);
        x(1)=pos(3);    
        
        annotation('arrow',x,y,'LineWidth',2);
        
            m=(x(2)-x(1))/(y(2)-y(1));
            b=x(1)-(m*y(1));
            
        xt=x;
        yt=y;
        
        x=x_to_norm(out.ton,tdec);
        y=y_to_norm(rmax+offset,rmax+offset);
        
        
        x(1)=xt(1);
        y(1)=yt(1);
            
        annotation('arrow',x,y,'LineWidth',2) 
        
        
        %         
        if isfield(out,'decaytime')

           idecay0=idecay0-1;
            idecay=find(t>=decay,1);
            x=x_to_norm(t(end),decay);
            y=ylim();
            
            y=y_to_norm(mean(y),ffit(idecay-idecay0)+offset);
            
            
            
              
            h=annotation('textbox','String',{['t_{0,decay} = ' num2str(decay) ' s'] ;['  \tau_{decay} =  ' num2str(out.decaytime) ' s']},'Color','black','HorizontalAlignment','left','BackgroundColor','w','LineStyle','none');drawnow;
          
            pos=get(h,'position');
            set(h,'position',[1-pos(3) y(1)-pos(4)/2 pos(3) pos(4)])
            
            y(1)=y(1);
            
            annotation('arrow',[1-pos(3) x(2)],y,'LineWidth',2);
            
            mt=min(decay+out.decaytime,t(end));
            
            idecay2=find(t>=(mt),1);         
           
%             y=ylim;
            xt=x_to_norm(mt,mt);
            m=(y(2)-y(1))/(x(2)-x(1));
            b=y(1)-(m*x(1));
            yt=m*(xt(1))+b;
            
            y=y_to_norm(ffit(idecay2-idecay0)+offset,ffit(idecay2-idecay0)+offset);
            y(1)=yt;
            x=xt;
            
            annotation('textarrow',x,y,'String','','LineWidth',2);
            
%             x=x_to_norm(decay,mt);
%             y=ylim;
%             y=y_to_norm(ffit(idecay-idecay0)+0.75*(y(2)-ffit(idecay-idecay0)),ffit(idecay-idecay0)+0.75*(y(2)-ffit(idecay-idecay0)));
%              annotation('line',x,y,'LineWidth',2);
        end
        
            hold off
    end
    
%     
%     Nosc=2;
%     argstr='Aoff,Aon,dact,dde,mact,nact,nde,rhoact,rhode,tdecay,ton,zA1,zA2,zm1,zm2,ztau1,ztau2,ztdriftswitch,zz';
%     prob={};
%     for i=1:Nosc
%         ii=int2str(i);
%         argstr=[argstr ',t0' ii  ',fp0' ii ',t1' ii ',f1' ii ',t2p' ii ',f2p' ii ];
%         prob={prob{:}, ['t1' ii]};
%     end
%     argstr=[argstr ',x'];
%     resposcfittype=fittype(['responseModelDrift11Osc(' argstr  ')'],'problem',prob)
%     
%     
    
    
    
 
    
    
   
%     out=characterizeDecay(tdecay,fdecay,decay,out);
       

%     out.decayrate=-((n*d^(-1+2*n))/(d^n+d^n)^2)+(n*d^(-1+n))/(d^n+d^n)+m;
    
    ind=t>=tfine(t50up);

    tcurr=t(ind);
    if length(tcurr)>=10
        dt=mean(tcurr(2:end)-tcurr(1:end-1));
        Tmax=35;
        n=0;




    Tmin=5;

  respflat=respest;
  respflat(i50up:i50down)=lr5050;
% respflat=respdrift;
detrend=fufine(ind,1)-respflat(ind);
  detrend(find(detrend(1:i50down-i50up+1)<0))=0;
  detrend(1)=0;
  detrend(end)=0;

    [crest,tcrest,wcrest,pcrest]=findpeaks(detrend,tcurr,'MinPeakProminence',3*sigma,'MinPeakDistance',Tmin,'WidthReference','HalfProm');
    [trough,ttrough,wtrough,ptrough]=findpeaks(-detrend,tcurr,'MinPeakProminence',3*sigma,'MinPeakDistance',Tmin,'WidthReference','HalfProm');

out.peakdeviations=length(crest);

   if length(crest)+length(trough)>1 % 
       
   iwcrest= find(tcrest(2:end)<[tcrest(1:end-1)+wcrest(1:end-1)]);
   wcrest(1+iwcrest)=min([wcrest(iwcrest),wcrest(1+iwcrest)],[],2);
   
   
   iwtrough= find(ttrough(2:end)<[ttrough(1:end-1)+wtrough(1:end-1)]);
   wtrough(1+iwtrough)=min([wtrough(iwtrough),wtrough(1+iwtrough)],[],2);
    
    
     [T,cohcrest,cohtrough,Nosc] = periodCluster4(abs(pcrest),tcrest,wcrest,abs(ptrough),ttrough,wtrough,dt,sigma,tcurr(1));

   
   
     
     

        goodcrest=find(cohcrest);
        goodtrough=find(cohtrough);
        
        magcrest=zeros(size(goodcrest));
        dutycrest=zeros(size(goodcrest));
        Tcrest=zeros(size(goodcrest));
        crestcounted=false(size(crest));
        
        Nosc=0;
        
        if length(goodcrest)>1
            
            j=find(ttrough<tcrest(goodcrest(1)),1,'last');
            k=find(ttrough>tcrest(goodcrest(1))&ttrough<tcrest(goodcrest(2)),1);
            jk=[j;k];
            if ~isempty(jk)
                magcrest(1)=crest(goodcrest(1))+mean(trough(jk));
            else
                magcrest(1)=pcrest(goodcrest(1));
            end
            
            for i=2:length(goodcrest)-1
                ii=goodcrest(i);
                j=find(ttrough<tcrest(ii)&ttrough>tcrest(goodcrest(i-1)),1,'last');
                k=find(ttrough>tcrest(ii)&ttrough<tcrest(goodcrest(i+1)),1);
                jk=[j;k];
                if ~isempty(jk)
                    magcrest(i)=crest(goodcrest(i))+mean(trough(jk));
                else
                    magcrest(i)=pcrest(goodtrough(i));
                end
                
                if goodcrest(i-1)==goodcrest(i)-1
                    crestcounted(goodcrest(i))=true;
                    Nosc=Nosc+1;
                    Tcrest(i)=diff(tcrest(ii-1:ii));
                    dutycrest(i)=wcrest(ii)/Tcrest(i);
                end
                
            end
            
            i=length(goodcrest);
            ii=goodcrest(i);
            j=find(ttrough<tcrest(goodcrest(i))&ttrough>tcrest(goodcrest(i-1)),1,'last');
            k=find(ttrough>tcrest(goodcrest(i)),1);
            jk=[j;k];
            if ~isempty(jk)
                magcrest(i)=crest(goodcrest(i))+mean(trough(jk));
            else
                magcrest(i)=pcrest(goodtrough(i));
            end

            if goodcrest(i-1)==goodcrest(i)-1
                crestcounted(goodcrest(i))=true;
                Nosc=Nosc+1;
                Tcrest(i)=diff(tcrest(ii-1:ii));
                dutycrest(i)=wcrest(ii)/Tcrest(i);
            end
            

        else
        end
            
        magtrough=zeros(size(goodtrough));
        Ttrough=zeros(size(goodtrough));
        
          if length(goodtrough)>1
            
            j=find(tcrest<ttrough(goodtrough(1)),1,'last');
            k=find(tcrest>ttrough(goodtrough(1))&tcrest<ttrough(goodtrough(2)),1);
            jk=[j;k];
            if ~isempty(jk)
                magtrough(1)=abs(trough(goodtrough(1))+mean(crest(jk)));
            else
                magtrough(1)=ptrough(goodtrough(1));
            end
            
            for i=2:length(goodtrough)-1
                ii=goodtrough(i);
                j=find(tcrest<ttrough(ii)&tcrest>ttrough(goodtrough(i-1)),1,'last');
                k=find(tcrest>ttrough(ii)&tcrest<ttrough(goodtrough(i+1)),1);
                jk=[j;k];
                if ~isempty(jk)
                    magtrough(i)=abs(trough(ii)+mean(crest(jk)));
                else
                    magtrough(i)=ptrough(ii);
                end
                
                if goodtrough(i-1)==ii-1
                    if ~crestcounted(j)
                        Nosc=Nosc+1;
                    end
                    jj=find(goodcrest==j);
                    if  ~isempty(jj) && (Tcrest(jj)==0 || abs(Tcrest(jj)-diff(ttrough(ii-1:ii)))/Tcrest(jj)<0.5)
                        Ttrough(i)=diff(ttrough(ii-1:ii));
                    end

                end
                
            end
            
            i=length(goodtrough);
            ii=goodtrough(i);
            j=find(tcrest<ttrough(ii)&tcrest>ttrough(goodtrough(i-1)),1,'last');
            k=find(tcrest>ttrough(ii),1);
            jk=[j;k];
            if ~isempty(jk)
                magtrough(i)=abs(trough(ii)+mean(crest(jk)));
            else
                magtrough(i)=ptrough(ii);
            end

                if goodtrough(i-1)==ii-1
                    if ~crestcounted(j)
                        Nosc=Nosc+1;

                    end
                    jj=find(goodcrest==j);
                    if ~isempty(jj) && (Tcrest(jj)==0 || abs(Tcrest(jj)-diff(ttrough(ii-1:ii)))/Tcrest(jj)<0.5)
                     Ttrough(i)=diff(ttrough(ii-1:ii));
                    end
                   

                end
            

        else
            
          end
          
          
           dutycycle=mean(wcrest(goodcrest))/mean(Tcrest(Tcrest~=0));
          
          if ~isempty(dutycycle) && isfinite(dutycycle) && (length(goodcrest)>2||length(goodtrough)>2)||((length(goodcrest)>1||length(goodtrough)>1)&&dutycycle>0.1) 
             
              %%%%
              %%EXPORT MEASURED PARAMETERS
              %%%
              
              Nosc=Nosc+1;
              out.Nosc=Nosc;
              out.periods=mean([Tcrest(Tcrest~=0);Ttrough(Ttrough~=0)]);
              out.Oscmag=mean([magcrest(magcrest~=0); magtrough(magtrough~=0)]);

              out.FWHMosc=mean(wcrest(cohcrest));
               out.dutycycle=dutycycle;
               
               out.oscillatoryPersistence=tcrest(goodcrest(end))-tcrest(goodcrest(1));
                
               
                if plt
                figure(1);

                for i=1:length(goodcrest)
                    
                    ii=goodcrest(i);
                    x=x_to_norm(tcrest(ii),tcrest(ii));
                    fi=fufine(find(t==tcrest(ii),1),1)+offset;
                    y=y_to_norm(fi+0.1*fmax,fi);
                    annotation('arrow',x,y) 

                end



                for i=1:length(goodtrough)
                    
                    ii=goodtrough(i);
                    x=x_to_norm(ttrough(ii),ttrough(ii));
                    fi=fufine(find(t==ttrough(ii),1),1)+offset;
                    y=y_to_norm(fi-0.1*fmax,fi);
                    annotation('arrow',x,y) 

                end

                hold on;
%                 plot(t,fufine(:,1));
                hold off;

                end
    
        
        
                   
               end
          end
        

        %there must be atleast three coherent periods to declare
        %someting unambigously periodic
        
 


    end %end enough data for MPR
    else
        
        disp('could not detect any activation!')
    end %end enough drift prior to activation



end %end sig diff from drift


 end %end function

