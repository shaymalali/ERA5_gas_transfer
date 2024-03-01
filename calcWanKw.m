function [kw_w]=calcWanKw(era,sc)

kw_w=NaN(size(era.wspd));
for i=1:size(era.wspd,3)
    % kw_w=0.251.*(era.u10(:,:,i)).^2.*((sc(:,:,i)./660).^(-1/2));
    kw_w(:,:,i)=(era.wspd(:,:,i)).^2.*((sc(:,:,i)./660).^(-1/2));   
end
