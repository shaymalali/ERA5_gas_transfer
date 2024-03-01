function [kmean]=calcClimateMean(k,lat,lon,dt)

[y1,m1]=ymd(dt); 

%Get the groups 
G=findgroups(m1);

for i=1:max(G)

   idx=find(G == i);

   for j=1:length(idx)
       s(:,:,j)=k(:,:,idx(j));
   end
   
   %kmean(:,:,i)=mean(s,3,"omitnan");
   kmean(:,:,i)=nanmean(s,3);
   clear s

end 
