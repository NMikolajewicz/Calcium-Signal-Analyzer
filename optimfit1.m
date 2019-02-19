function [f,gof] = optimfit(x,y,plotflag)

if nargin==2
    plotflag=false;
end
if length(x)<2
    f=@(x) repmat(y,size(x));
    gof=[];
    return;
end

[flin,goflin]=fit(x,y,'poly1','Lower',[-Inf,0],'Upper',[Inf,Inf]);
AIClquad=4+2*goflin.sse;

if plotflag
    plot(x,y,'ob',x,flin(x),'-r');
    legend('Data','Linear Fit');
end

if length(x)>=3
AIClquad=AIClquad+12/(length(x)-3);

if mean(diff(y))>0
    b0=2/(x(end)-x(1));
    a0=(y(end)-y(1))/(exp(2)-1);
    c0=y(1);
else
    b0=-2/(x(end)-x(1));
    a0=(y(end)-y(1))/(exp(-2)-1);
    c0=y(1);
    
end
x0=x(1);
[fexp,gofexp]=fit(x,y,@(a,b,c,x) a*(exp(b*(x-x0))-1)+c,'Upper',[Inf,Inf,Inf],'Lower',[-Inf,-Inf,-Inf],'Start',[a0,b0,c0]);
AICexp=6+2*gofexp.sse;
if length(x)>=4
AICexp=AICexp+24/(length(x)-4);
end

if plotflag
    hold on;
    plot(x,fexp(x),'-g');
    legend('Data','Linear Fit','Exponential Fit');
    hold off;
end

if AICexp<AIClquad
    f=fexp;
else
    f=flin;
end

else
    f=flin;
end


end

