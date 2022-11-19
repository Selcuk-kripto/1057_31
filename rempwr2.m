function b=rempwr2(b,j,N)

for i=1:j
    b=rem(b*2,N);
end