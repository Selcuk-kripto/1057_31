function b=trc(i,P)
% rempwr2(i,j,N): remainder of i*2^j divided by N
n=size(P,2);
b=0;
N=2^n-1;
for j=0:n-1
    b=xor(b,P(rempwr2(i,j,N)+1,:));
end