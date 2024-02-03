function [V_y,Q_y]=eigsort(R1, R2)




[L,L1]=size(R1);

if (nargin == 2)
    [V_y,Q_y] = eig(R1, R2);
else
    [V_y,Q_y] = eig(R1);
end

for nnn=1:L
    for nn=1:L-1
        if Q_y(nn,nn)<Q_y(nn+1,nn+1)
            temp_1=Q_y(nn,nn);
            Q_y(nn,nn)=Q_y(nn+1,nn+1);
            Q_y(nn+1,nn+1)=temp_1;
            V_temp=V_y(:,nn);
            V_y(:,nn)=V_y(:,nn+1);
            V_y(:,nn+1)=V_temp;
        end
    end
end