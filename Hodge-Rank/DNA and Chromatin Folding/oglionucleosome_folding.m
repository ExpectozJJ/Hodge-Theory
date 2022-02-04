clear

threshold = 10;
dat=pdb2mat('330d-P.pdb');
%coordinates of protein frames
data(:,1)=dat.X;
data(:,2)=dat.Y;
data(:,3)=dat.Z; 

natom=length(data);

no_vertex = natom;
W = zeros(no_vertex,no_vertex);
Y = zeros(no_vertex,no_vertex);

for i=1:(no_vertex-1)
    for j=(i+1):no_vertex
        amino1 = data(i,:);
        amino2 = data(j,:);
        dis = sqrt(sum((amino1-amino2).^2));
        if dis < threshold
            W(i,j) = 1;
            W(j,i) = 1;
            Y(i,j) = dis;
            Y(j,i) = -dis;
        end
    end
end

Yvec = uppertri(Y,W);

%Edge (s0)
s0 = [];
for i=1:(no_vertex-1)
    for j=(i+1):no_vertex
        if W(i,j)~=0
            s0 = [s0;[i,j]];
        end
    end
end
s0 = s0';
no_edge = size(s0,2);
 
%Triangle (s1)
s1 = [];
for i=1:(no_vertex-2)
    for j=(i+1):(no_vertex-1)
        for k=(j+1):no_vertex
            if ((W(i,j)~=0)&&(W(j,k)~=0)&&(W(k,i)~=0)) %if there is such a triangle
                s1 = [s1;[i,j,k]];
            end
        end
    end
end
s1 = s1';
no_triangle = size(s1,2);

%Coboundary matrix (d0) edge-vertex
d0 = zeros(no_edge,no_vertex);
for i=1:no_edge
    d0(i,s0(1,i)) = -1;
    d0(i,s0(2,i)) = 1;
end

%Laplacian matrix (L0) vertex-vertex
L0 = (d0')*(d0);

%Coboundary matrix (d1) triangle-edge
d1 = zeros(no_triangle,no_edge);
for i=1:no_triangle
    tri1 = s1(1,i);
    tri2 = s1(2,i);
    tri3 = s1(3,i);
    d1(i,find(s0(1,:) == tri1 & s0(2,:) == tri2)) = 1;
    d1(i,find(s0(1,:) == tri1 & s0(2,:) == tri3)) = -1;
    d1(i,find(s0(1,:) == tri2 & s0(2,:) == tri3)) = 1;
end

%Laplacian matrix (L1) triangle-triangle
L1 = (d1)*(d1');

%Divergence
div = (-d0')*(Yvec);

%Optimization (to find the gradient flow)
s = lsqr(L0,-div,1e-10,50000);

%Gradient Flow
Yg = d0*s;

%Triangular Curl
curl = d1*Yvec;

%Optimization (to find the curl flow)
z = lsqr(L1,curl,1e-10,50000);

%Curl Flow
Yc = (d1')*(z);

%Harmonic Flow
Yh = Yvec-Yg-Yc;

%Inconsistencies Rate
%local_in2 = sum(abs(Yc))/no_edge;
%global_in2 = sum(abs(Yh))/no_edge;

local_in2 = sqrt(sum((Yc./Yvec).^2));
local_in1 = sum(abs((Yc+Yh)./Yvec))
local_in1avg = sum(abs((Yc+Yh)./Yvec))/no_vertex