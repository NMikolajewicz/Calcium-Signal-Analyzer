function [T,cohcrest,cohtrough,Nosc] = periodCluster4(crest,tcrest,wcrest,trough,ttrough,wtrough,ts,sigma,t0)
global tend


ncrest=length(tcrest);
ntrough=length(ttrough);
ndt=ncrest+ntrough-2;
dt0=[diff(tcrest);diff(ttrough)];
tfor0=[tcrest(2:end);ttrough(2:end)];
tback0=[tcrest(1:end-1);ttrough(1:end-1)];
[tfor,st]=sort(tfor0);
tback=tback0(st);
dt=dt0(st);

if length(dt)<=1
    T=dt;
    cohcrest=true(size(crest));
    cohtrough=true(size(trough));
    Nosc=1;
    return;
end

% dt0=[diff(tcrest);diff(ttrough);];

wback0=[wcrest(1:end-1);wtrough(1:end-1)];
wback=wback0(st);

wfor0=[wcrest(2:end);wtrough(2:end);];
wfor=wfor0(st);

Afor0=[crest(2:end);trough(2:end)];
Afor=Afor0(st);

Aback0=[crest(1:end-1);trough(1:end-1)];
Aback=Aback0(st);

pc=find(st<ncrest&[st(2:end);st(end)]>=ncrest);
for i=1:length(pc)
    ii=pc(i);
    Afor(ii:ii+1)=mean(Afor(ii:ii+1));
    Aback(ii:ii+1)=mean(Aback(ii:ii+1));
end


nmax=max(ncrest,ntrough);
options = statset('MaxIter',2^(ncrest+ntrough+1),'Display','off'); 

    function [copt,kopt]=optClusterLocalExtrema(dt)
        
        if length(dt)>2
            if length(dt)>3
                X=[dt,wfor,wback];
                kmax=ceil(nmax/2)+1;
            else
                X=[dt,(wfor+wback)/2];
                kmax=2;
            end

    %         kmax=ndt-1;

            Rmax=ceil((nmax+1)^2);
            clust=zeros(length(dt),kmax);

            gmfit={};
            cl={};
            for k=1:kmax

                try
                    gmfit{k} = fitgmdist(X,k,'CovarianceType','diagonal','SharedCovariance',false,'Start','plus','Replicate',Rmax,'Regularization',ts/5,'Options',options);
                    cl{k}=cluster(gmfit{k},X);
                    if sum(~ismember(1:k,cl{k}))==0
                        clust(:,k)=cl{k};
                    else
                        kmax=k-1;
                        break;
                    end
                catch e
                    getReport(e)
                    kmax=k-1;
                    break;
                end
            end

            clustfun=@(x,k) clust(:,k);

            try
                eva=evalclusters(dt,@(x,k) clustfun(x,k),'gap','Klist',[1:kmax]);
            catch ex
                getReport(ex)
            end

            kopt=eva.OptimalK;
            copt=clust(:,kopt);        
            fopt=gmfit{kopt};
        else
            copt=[1;1];
            kopt=1;
        end
        
%         cf(25);
%         subplot(2,1,1);
%         h=gscatter(dt,wfor,cl{kopt});
%           subplot(2,1,2);
%           wave=(wfor+wback)/2;
%         h=gscatter(tfor,dt,cl{kopt});
%         drawnow;
    end



[c,kopt]=optClusterLocalExtrema(dt);
stagnant = false;

if kopt==1
    stagnant=true;
end

dist=sort(sum(repmat(c,1,kopt)==repmat(1:kopt,size(c))),'descend');

cold=c;
distold=dist;


coh=false(size(dt));
incoh=false(size(dt));
% %sequential maxima in the same cluster are deemed coherent
% coh(1:ncrest-1)=[c(1:ncrest-2)==c(2:ncrest-1);c(ncrest-1)==c(ncrest-2)];
% 
% %sequential minima in the same cluster are deemed coherent
% coh(ncrest:ndt)=[c(ncrest:end-1)==c(ncrest+1:end);c(end)==c(end-2)];


rho=(sigma/max([crest;trough]))^(ts/mean(dt));
rho=min(rho,1);

% f=std(dt)/mean(dt);+(1-rho)*0.25;
f=rho*max(std(dt)/mean(dt),0.5)+(1-rho)*0.25;

%first pass coherence check;
[a,b]=hist(c,unique(c));
alone=false(size(c));
for i=1:length(b)
    if a(i)==1
        alone(c==b(i))=true;
    end
end
hanging=((tback-wback)<t0)|((tfor+wfor)>=tend);
cohw=(((wfor>f/2*dt&wfor<(3*f/2)*dt)&(wback>f/2*dt&wback<(3*f/2)*dt)));

coh1=([c(1:end-1)==c(2:end);c(end)==c(end-1)]|alone)&(cohw|hanging);



uncertdt=std(dt)/mean(dt);


if sum(coh1)~=0
    Atime=optimfit1(tfor(coh1),Afor(coh1),false);
    uncert=periodicUncertainty(sigma./Atime(tfor),dt/ts);
    rho=(sum(coh1)/length(coh1)).^(sigma./Atime(tfor));
    uncert=rho.*uncert+(1-rho).*uncertdt;
    
    cohw2=(((wfor>(f/2-6*uncert).*dt&wfor<(3*f/2+6*uncert).*dt)&(wback>(f/2-6*uncert).*dt&wback<(3*f/2+6*uncert).*dt)));
else
    cohw2=coh1;
end


coh2=([c(1:end-1)==c(2:end);c(end)==c(end-1)]|alone)&(cohw2|hanging);

% cf(57);
% subplot(2,1,1);

fplane=optimfit1(wfor(coh2&~hanging),dt(coh2&~hanging),false);


% subplot(2,1,2);
fplane2=optimfit1(wback(coh2&~hanging),dt(coh2&~hanging),false);


%extend the coherence based on a trend in the width-dt plane
%0.3 is an arbitrary parameter ~ 3x timing uncertainty
%of oscillations with 50% white noise with adequate sampling

begins=st==1|st==ncrest;
ends=st==ncrest-1|st==(ncrest+ntrough-1);
central=~begins&~ends;

% cohw1=(wfor>f*dt&wfor<(1-f)*dt)&(wback>f*dt&wback<(1-f)*dt);
% cohw2=true(size(cohw1));
% for i=1:length(cohw1)
%     if i>1 && cohw1(i-1) && cohw1(i)
%         cohw2(i)=cohw1(i)
%     elseif i<length(cohw1)
%         cohw2=(wfor(i)>f*dt&wfor(i)<(1-f)*dt)
%     end
% end
% cohwfor=(wfor>f*dt&wfor<(1-f)*dt);
% cohwback=   

if sum(coh1&coh2&~hanging)~=0
uncert=periodicUncertainty(sigma./Atime(tfor),dt/ts);
rho=(sum(coh1)/length(coh1)).^(sigma./Atime(tfor));
uncert=rho.*uncert+(1-rho).*uncertdt;


cohw3=abs((fplane(wfor)-dt)./dt)<6*uncert&cohw2;

cohw4=abs((fplane2(wback)-dt)./dt)<6*uncert&cohw2;

else
    
    cohw3=coh1;
    cohw4=coh1;

end


if sum(cohw3|cohw4)~=0
    Atime2=optimfit1(tfor(cohw3|cohw4),Afor(cohw3|cohw4),false);
    cohA2=(Atime2(tfor)-Afor)/std((Atime2(tfor)-Afor))<2;
    Atime3=optimfit1(tback(cohw3|cohw4),Aback(cohw3|cohw4),false);
    cohA3=(Atime3(tback)-Aback)/std((Atime3(tback)-Aback))<2;
else
    cohA2=cohw3;
    cohA3=cohw3;
end





Atime4=optimfit1(tfor(cohA2|cohA3),Afor(cohA2|cohA3),false);

if sum(cohw3|cohw4)~=0

    uncert2=periodicUncertainty(sigma./Atime2(tfor),dt/ts);
    rho=(sum(cohw3|cohw4)/length(coh)).^(sigma./Atime2(tfor));
    uncert2=rho.*uncert2+(1-rho).*uncertdt;
    % cf(58);
    fper=optimfit1(tfor(cohw3|cohw4),dt(cohw3|cohw4),false);
    fper2=optimfit1(tback(cohw3|cohw4),dt(cohw3|cohw4),false);

    cohh1=abs((fper(tfor)-dt)./dt)<3*uncert2&hanging;

    cohh2=abs((fper2(tback)-dt)./dt)<3*uncert2&hanging;

else
    cohh1=cohw3;

    cohh2=cohw3;
end

coh=(cohw3|cohw4|cohh1|cohh2)&(cohA2|cohA3);


% cf(57);
% fper=optimfit2(tfor(coh|alone),dt(coh&~incoh),true);




incoh=st(~coh);
crestfor=true(size(tcrest));
crestback=true(size(tcrest));

troughfor=true(size(ttrough));
troughback=true(size(ttrough));
%cohst(cohst>ncrest-1)-(ncrest-1)
for i=1:length(incoh)
    ii=incoh(i);
    if ii<ncrest
        crestback(ii)=false;
        crestfor(ii+1)=false;
        if ii==1
            crestfor(ii)=false;
        elseif ii==ncrest-1
            crestback(ii+1)=false;
        end
    else
        ii=ii-ncrest+1;
        troughback(ii)=false;
        troughfor(ii+1)=false;
        if ii==1
            troughfor(ii)=false;
        elseif ii==ntrough-1
            troughback(ii+1)=false;
        end

    end
end

cohcrest=true(size(tcrest));
cohcrest(~crestfor&~crestback)=false;

% st<ncrest-1

cohtrough=true(size(ttrough));
cohtrough(~troughfor&~troughback)=false;

T=mean(dt(coh));

dtcoh=diff(tfor(coh));
Nosc=0;

goodcrest=find(cohcrest);
goodtrough=find(cohtrough);

Nosc=length(goodcrest);
tcrest(cohcrest);

for i=1:length(goodtrough)-1
    if isempty(tcrest(cohcrest)>ttrough(goodtrough(i))&tcrest(cohcrest)<ttrough(goodtrough(i+1)))
        Nosc=Nosc+1;
    end
end

% if ~cohcrest(end)
%     disp('gone down the poop hole!!')
% end

end

