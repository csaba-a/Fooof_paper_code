function dataES=getES(datamucl,datasucl,datamnorm,datasnorm,flag)


%flag ==0 effect size returned as cohens d if possible, then z if not, then nan
%otherwise

%flag==1 Effect size returned as z score of pat to norm

nroi=size(datamucl,1);
dataES=nan(nroi,1);
if flag==0
for i=1:nroi
    if datasucl(i)~=0 && datasnorm(i)~=0
        dataES(i)=(datamucl(i)-datamnorm(i))/mean([datasucl(i),datasnorm(i)]);
    elseif datasucl(i)==0 && datasnorm(i)==0
        dataES(i)=nan;
    elseif datasucl(i)==0 && datasnorm(i)~=0
        dataES(i)=datamucl(i)-datamnorm(i)/datasnorm(i);
    elseif datasucl(i)~=0 && datasnorm(i)==0
        dataES(i)=datamucl(i)-datamnorm(i)/datasucl(i);
    end
end
elseif flag==1
    dataES=(datamucl-datamnorm)./datasnorm;
end


end