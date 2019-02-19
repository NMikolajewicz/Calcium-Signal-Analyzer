function  writeSheet(f,s,expmts)



javaaddpath('poi_library/poi-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('poi_library/xmlbeans-2.3.0.jar');
javaaddpath('poi_library/dom4j-1.6.1.jar');
javaaddpath('poi_library/stax-api-1.0.1.jar');




map=struct();
map.amplitude='Amplitude (A.U.)';
map.ton='Time of Onset (s)';
map.tact='Activation Time (s)';
map.FWHM='FWHM (s)';
map.AUC='Area Under Curve';
map.decaytime='Decay Time (s)';


map.Nosc='Number of Oscillations';
map.Oscmag='Magnitude of Oscillations';
map.periods='Period';
map.periodstd='Period Standard Deviation';
map.FWHMosc='FWHM of Oscillation';
map.oscillatoryPersistence='Persistence of Oscillation';
map.dutycycle='Duty Cycle';


map.NoscCOH='Number of Oscillations (Coherent)';
map.OscmagCOH='Magnitude of Oscillations (Coherent)';
map.periodsCOH='Period (Coherent)';
map.periodstdCOH='Period Standard Deviation (Coherent)';
map.FWHMoscCOH='FWHM of Oscillation (Coherent)';
map.oscillatoryPersistenceCOH='Persistence of Oscillation (Coherent)';
map.dutycycleCOH='Duty Cycle (Coherent)';



map.peakdeviations='Number of deviation peaks';


keys=fieldnames(map);
maxLength=zeros(size(keys));
count=0;
for i=1:length(expmts)
    for j=1:length(expmts{i}.resps)
        count=count+1;
        for k=1:length(keys)
            if ~isempty(expmts{i}.resps{j}) && isfield(expmts{i}.resps{j},keys{k})
                    maxLength(k)=max(maxLength(k),length(expmts{i}.resps{j}.(keys{k})));
            end
        end
    end
end

hdr={{'Recording'}};
for i=1:length(keys)
    for j=1:maxLength(i)
        
        h=map.(keys{i});
        
        if maxLength(i)>1
            h=strsplit(map.(keys{i}),'(');
            h{end-1}=[h{end-1} ' ' int2str(j) ' '];
            h{end}=strtrim(h{end});
            h=strjoin(h,'(');
        end

        hdr{length(hdr)+1}={h};
    end
end


prnt=cell(count+1,length(hdr));
prnt(1,:)=hdr;
count=1;
for i=1:length(expmts)
    expmnt=expmts{i};
    for j=1:length(expmnt.resps)
        count=count+1;
        resp=expmnt.resps{j};
        tmp=cell(1,length(hdr));
        tmp(1)={['Exp:' int2str(i) ';ROI:' int2str(j),]};
        if isempty(resp)
            tmp(2:length(hdr))={'-'};
%             tmp={{['Exp:' int2str(i) ';ROI:' int2str(j),]},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'},{'-'}};
        else
            pos=2;
            for k=1:length(keys)
                if isfield(resp,keys{k})
                    lr=length(resp.(keys{k}));
                    tmp(pos:pos+lr-1)=num2cell(resp.(keys{k}));
                else
                    lr=0;
                end
                if lr<maxLength(k)
                    tmp(pos+lr+1:pos+maxLength(k)-1)={'-'};
                end
                pos=pos+maxLength(k);
            end
        end
        prnt(count,:)=tmp;
    end

end



    xlwrite(f,prnt,s);

end

