function [point]=kmedoid_algorithm(Y,n)

%Band Normalization
Y1=Y;
[Bands,num] = size(Y);
for i=1:1:Bands
    R=Y1(i,:);
    meanr=mean(R);
    stdr=std(R);
    R1=(R-meanr)/stdr;
    Y(i,:)=R1;
end

%Band selection
R = (Y*Y')/num;
u = mean(Y,2);
K = R-u*u';
K1 = K.*eye(Bands);
K = K - K1;
imagesc(K);
[q1,q2]=max(K);
[q3,q4]=max(q1);
n_y1=q2(q4);
n_x1=q4;

band_x1=Y(n_x1,:)';
band_y1=Y(n_y1,:)';

%Convex set for the selected bands
p11=[band_x1,band_y1];
for i=0:0.001:1
    conv_points = boundary(p11,i);
    [m1,n1]=size(conv_points);
    if(m1>n)
        break
    end
end

cov_ED=conv_points;

z=(size(cov_ED)-1);
for i=1:z(1)
    [t1]=p11(cov_ED(i),:);
    [t2]=p11(cov_ED(i+1),:);
    dist(i,1)=sqrt(sum((t1-t2).^2));
end
e=[p11(cov_ED),p11(cov_ED,2)];
e(end,:)=[];
cov_ED(end,:)=[];

%Extra points removal using Kmedoids
vec=kmedoids(e,n,'Distance','euclidean');
point=[];
for k=1:1:n
    for l=1:1:size(vec)
        if (vec(l)==k)
            point=[point;cov_ED(l)];
            break;
        end
    end
end

end