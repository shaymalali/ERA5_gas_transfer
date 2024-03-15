function [alpha, min_x]=solveparabola(x,y,method)
%inputs: x: the coefficients
%        y: the rms errors
%        method: 0 if A\B ; 1 if using the long langrage interpolation method


if method == 0
  A=[x.^2 x repelem(1,5,1)];
  alpha=A\y;

end

if method == 1

  for i=1:length(x)
    p=1;
    for j=1:length(x)
        if j~=i
            c = poly(x(j))/(x(i)-x(j));
            p = conv(p,c);
        end
    end
    term = p*y(i);
    sum= sum + term;
   end
   alpha=sum;
end

%find the minimum by diffrentiating with respect to x
% 2ax +b =0
% x = -b/2a

min_x= -alpha(2)/2*alpha(1);