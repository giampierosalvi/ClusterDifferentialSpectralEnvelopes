function d=bhattacharyya(mu1, mu2, C1, C2)
% BHATTACHARYYA  Bhattacharyya distance between two Gaussian classes
%
% 2010, Giampiero Salvi <giampisalvi@gmail.com>

%Check inputs and output
error(nargchk(4,4,nargin));
error(nargoutchk(0,1,nargout));

% check dimensions
m = size(C1,1);
assert(size(C1,2)==m, 'C1 must be square.')
assert(size(C2,1)==size(C2,2), 'C2 must be square.')
assert(size(mu1,2)==size(mu2,2), 'Dimensions of M1 and M2 mismatch.');
assert(size(mu1,2)==size(C1,2), 'Dimensions of M1 and C1 mismatch.');
assert(size(mu2,2)==size(C2,2), 'Dimensions of M2 and C2 mismatch.');

C=(C1+C2)/2;
dmu=(mu1-mu2)/chol(C);
try
    d=0.125*(dmu*dmu')+0.5*log(det(C/chol(C1*C2)));
catch
    d=0.125*(dmu*dmu')+0.5*log(abs(det(C/sqrtm(C1*C2))));
    warning('MATLAB:divideByZero','Data are almost linear dependent. The results may not be accurate.');
end
% d=0.125*dmu*dmu'+0.25*log(det((C1+C2)/2)^2/(det(C1)*det(C2)));
