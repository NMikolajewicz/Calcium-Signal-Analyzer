function expmnts = splitSheet(f,sheet)

Data=xlsread(f,sheet);
count=0;
i=1;
expmnts=cell(0,0);
begin=1;
blocks=false;
while i<size(Data,2)
    if sum(isnan(Data(:,i)))~=0
        blocks=false;
        break;
    end
    i=i+1;
end

i=1;

if blocks

    while i<size(Data,1)
        if sum(isnan(Data(i,:)))~=0
            disp(Data(i,:))
            break;
        end
         disp(Data(i,:))
        begin=begin+1;
        i=i+1;
    end
end

begin=i;
i=1;
Data=Data(begin:end,:);
while i<=size(Data,2)
    if Data(1,i)==0%new experiment
        i0=i;
        count=count+1;
        l=0;
        i=i+1;
        while i<=size(Data,2) && sum(~isnan(Data(:,i)))~=0
            i=i+1;
            l=l+1;
        end
        ind=~isnan(Data(:,i0));
        for j=i0+1:i0+l
            Data(ind,j)=inpaintn(Data(ind,j));
        end
        expmnts{count}=struct('time',Data(ind,i0),'fluo',Data(ind,i0+1:i0+l));
    end
    i=i+1;
end

end

