function [d14c]=load_results(path)
d14c=[]; 
for i=1:272 %change length to number of files in the model output directory
        fname=sprintf('model_C12_C13_C14_%i.mat',i);
        s=[path '/' fname];
        
        load(s)

        % solve for Delta C14
        R14oxa=1.176*1e-12;
        R13pdb=1.12372*1e-2;
        R13 = @(c13,c12) c13./(c13+c12);
        R14 = @(c14,c13,c12) c14./( c13 + c12 );
        delta13c = @(c13,c12) ( R13(c13,c12) ./ ( (1 - R13(c13,c12) ) * R13pdb ) - 1 ) * 1e3;
        Delta14c = @(c14,c13,c12) (( R14(c14,c13,c12) ./ R14oxa ) .* ( 0.975 ./ ( 1 + delta13c(c13,c12) / 1e3 ) ).^2 -1 ) * 1e3;

        nwet=200160; %number of wet points in ocean
        
        DIC=Xout.C12(1:nwet,:);
        DI13C=Xout.C13(1:nwet,:);
        DI14C=Xout.C14(1:nwet,:);

        tmp=Delta14c(DI14C,DI13C,DIC);
        d14c=[d14c; tmp(:)];   

end




