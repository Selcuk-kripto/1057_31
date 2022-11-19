clear all

n=15;
N=2^n-1;
%Enter kr = 1, 3, or 5 for the PW type 1-RSBFs, 3-RSBFs, or 5-RSBFs, resp.
kr=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%
% p represents the polynomial x^15+x+1
p(1:n)=0;
p(n-1:n)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%
% P stores the vectorial representation of the field elements
% dobin(a,n): returns n-bit representation of the integer a
P(1:N,1:n)=0;
P(1,:)=dobin(1,n);
Pd(1)=1;
x(1:n)=0;
x(n)=1;
t=x;
for i=1:2^n
    if x(1)==1
        y=[x(2:n) 0];
        y=xor(y,p);
    else
        y=[x(2:n) 0];
    end
    if sum(y==t)==n
        i
        break
    end
    x=y;
    P(i+1,:)=y;
    Pd(i+1)=todec(y);
end

% TR is the output string of the trace function Tr_1^n(\zeta), where n=15,
% i=0,1,...,2^15-2, and \zeta is a root of x^15+x+1
% trc(i-1,P): returns Tr_1^15(\zeta^(i-1))
TRb(1:N,1:n)=0;
for i=1:2^n-1
    TRb(i,:)=trc(i-1,P);
end

TR=TRb(:,n)';

% Adk: (1057,31)-interlaved sequence form of TR
tp=5;
Ntp=2^tp-1;
Mtp=N/Ntp;
k=1;
for i=1:Ntp
    Adk(i,1:Mtp)=TR(k:k+Mtp-1);
    Adkt(i,1:Mtp)=Pd(k:k+Mtp-1);
    k=k+Mtp;
end

% Z: Column numbers of 'all zero' columns of Adk
% Kzr: Number of 'all zero' columns
k=1;
for i=1:Mtp
    if sum(Adk(:,i)==0)==Ntp
        Z(k)=i-1;
        k=k+1;
    end
end
Kzr=k-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%
% E(k): the smallest integer in the kth equivalence class
% SZ(k): the size of the kth equivalence class
% EC(k,:): the elements of the kth equivalence class
Nr=1057;

B(1:Nr)=0:Nr-1;
k=0;
SZa(1:n)=0;
ec=1;
for i=1:Nr
    if B(i)~=-1
        k=k+1;
        j=0;
        jt=0;
        while (jt>=0)
            if B(rempwr2(i-1,jt,Nr)+1)==-1;
                break;
            else
                B(rempwr2(i-1,jt,Nr)+1)=-1;    
            end 
            EC(k,j+1)=rempwr2(i-1,jt,Nr);
            j=j+1;
            jt=jt+kr;
        end
        SZa(j)=SZa(j)+1;
        SZ(k)=j;
        if j==3
            EC3(ec)=k;
            ec=ec+1;
        end        
    end
end

E=EC(:,1)';
for i=1:n
    if SZa(i)~=0
        fprintf('\nNumber of equivalence class(es) with size %2d: %d',i,SZa(i));
    end
end
fprintf('\nThe representatives of these equivalence classes are:\n');
for i=1:k
    fprintf('%d,',E(i));
end
fprintf('%d\n',E(k));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 4 %%%%%%%%%%%%%%%%%%%%%%%%
% E(L(r)) is the representative of the equivalence class containing j
Kec=size(EC,1);

for i=1:Kec
    for j=1:SZ(i)
        L(EC(i,j)+1)=i;
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 5 %%%%%%%%%%%%%%%%%%%%%%%%

for I=1:Kec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 6 %%%%%%%%%%%%%%%%%%%%%%%%

    C(1:Kec)=0;
    K(1:Kzr)=0;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 7 %%%%%%%%%%%%%%%%%%%%%%%%    

    for i=1:Kzr
        K(i)=mod(Z(i)-E(I),Mtp);
        m=mod(K(i),Nr);
        C(L(m+1))=C(L(m+1))+1;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 8 %%%%%%%%%%%%%%%%%%%%%%%%

    CA(I,:)=C;
end
% The coefficients of the system of linear inequalities are saved in the
% file 'CFs.txt'
fid = fopen('CFs.txt', 'w');
for i=1:size(CA,1)
    for j=1:size(CA,2)
        fprintf(fid,'%d ',CA(i,j));
    end
end
fclose(fid);
return

