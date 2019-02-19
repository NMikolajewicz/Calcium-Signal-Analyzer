function out = characterizeDecay(t,f,decay,Amax,in)
global plt;
out=in;

if length(t)>=4
    decaymax=max(f);
    decaymin=min(f);
    inde=find(f-decaymin<(decaymax-decaymin)/exp(1),1);
    ts=t(inde);
    if isempty(ts)
        ts=t(end)/2;
    end
    w=ones(size(f));
    w(f>f(1))=0;
    [defit,gof]=fit(t,f,@(a,b,c,m,x) a*exp(-x/b)+c+m*(x+decay),'Lower',[0,0,-Inf,-Inf],'Upper',[Amax,2*t(end),Inf,Inf],'StartPoint',[(decaymax-decaymin),ts,f(end),(f(end)-f(1))/t(end)],'Weight',w);
   [defit2,gof2]=fit(t,f,@(a,b,c,x) a*exp(-x/b)+c,'Lower',[0,0,-Inf],'Upper',[Amax,2*t(end),Inf],'StartPoint',[(decaymax-decaymin),ts,f(end)],'Weight',w);

   
    if gof2.sse<gof.sse
        w=ones(size(f));
        w(abs(defit2(t)-f)>gof2.rmse)=0;
    else
        w=ones(size(f));
        w(abs(defit(t)-f)>gof.rmse)=0;
    end
    
    [defit,gof]=fit(t,f,@(a,b,c,m,x) a*exp(-x/b)+c+m*(x+decay),'Lower',[0,0,-Inf,-Inf],'Upper',[Amax,2*t(end),Inf,Inf],'StartPoint',[(decaymax-decaymin),ts,f(end),(f(end)-f(1))/t(end)],'Weight',w);
    [defit2,gof2]=fit(t,f,@(a,b,c,x) a*exp(-x/b)+c,'Lower',[0,0,-Inf],'Upper',[Amax,2*t(end),Inf],'StartPoint',[(decaymax-decaymin),ts,f(end)],'Weight',w);

    
%     if plt
%         cf(3);
%         plot(t,f,'.b',t,defit(t),'-r',t,defit2(t),'-g');
%         legend('Data','Exp + Lin','Exp')
% %         plot(t,f,'.b',t,defit(t),'-r');
% %         legend('Data','Exp + Lin')
%     end
    
   AIC1=2*4+2*gof.sse*length(t);
   AIC1c=AIC1-2*4*(4+1)/(length(t)-4-1);
   AIC2=2*3+2*gof2.sse*length(t);
   AIC2c=AIC2-2*3*(3+1)/(length(t)-3-1);
   
   if AIC2c<AIC1c
        fitting = defit2;
   else
       fitting=defit;
   end
   
    c=min(fitting(t));
%     c =     defit(t(imin));
    A = fitting(t(1))-c;
    idecay = find((fitting(t)-c)<A/exp(1),1);
    
    if AIC2c>AIC1c
        
    if defit.m<0

            out.decaytime = max(t(idecay),defit.b);

    else
        out.decaytime = min(t(idecay),defit.b);
    end
    
    else
        if defit2.b<t(end)
            out.decaytime = max(t(idecay),defit2.b);
        else
            out.decaytime=max(t(idecay),defit.b);
        end
    end
    
end
end