clc
clear

TADindexValues = [108 33 114 20 42 92]; % TAD indices
indexLength = size(TADindexValues,2);

chrno = 'chr10';

%Domain Extraction
[chr,do_start,do_end] = textread('rep2.txt','%s %f %f');
chr_edit = [];
for i=1:length(chr)
     x=string(chr{i});
     chr_edit = [chr_edit,x];
end
do_start = round(do_start/40000); %40kb bp is the resolution of the data
do_end = round(do_end/40000);
chr_index = find(chr_edit == chrno);
domain_data = [do_start(chr_index),do_end(chr_index)];

%Data Extraction
data = load('uij.chr10');

%oneDistanceALIs = zeros(6,1); %for constant edge flow model
normalDistanceALIs = zeros(6,1);
threshold = 0.4;
alpha = 0.25;

for u = 1:indexLength
domain = domain_data(TADindexValues(u),:);

final_data = data(domain(1):domain(2),domain(1):domain(2));
throw_ind = find(sum(final_data) == 0);
final_data(throw_ind,:) = [];
final_data(:,throw_ind) = [];

M = final_data;

%Matrix Representation of Graph (W) and Pairwise Comparison Matrix (Y)
no_vertex = size(final_data,1);

W = zeros(no_vertex);
Y = zeros(no_vertex);

%Transformation of Hi-C data to alpha-distance matrix
for i=1:(no_vertex-1)
    for j=(i+1):no_vertex
        if M(i,j) > 0
            Y(i,j) = 1/(M(i,j)^alpha);
            %Y(i,j) = 1; % for constant edge flow model
        elseif M(i,j) == 0
            Y(i,j) = Inf;
        end
    end
end

%Implementing cutoff distance (threshold)
for i=1:(no_vertex-1)
    for j=(i+1):no_vertex
        if Y(i,j) < threshold
            W(i,j) = 1;
            W(j,i) = 1;
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
            if ((W(i,j)~=0)&&(W(j,k)~=0)&&(W(k,i)~=0))
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

%required for minimization of curl flow
L1 = (d1)*(d1');

%Divergence
div = (-d0')*(Yvec);

%Optimization (to find the gradient flow)
s = lsqr(L0,-div,1e-15,50000);

%Gradient Flow
Yg = d0*s;

%Triangular Curl
curl = d1*Yvec;

%Optimization (to find the curl flow)
z = lsqr(L1,curl,1e-15,50000);

%Curl Flow
Yc = (d1')*(z);

%Harmonic Flow
Yh = Yvec-Yg-Yc;

% %Inconsistencies 
% local_in1 = sum(abs(Yc))/sum(abs(Yvec))
% global_in1 = sum(abs(Yh))/sum(abs(Yvec));
% 
% local_in2 = sum(Yc.^2)/sum(Yvec.^2)
% global_in2 = sum(Yh.^2)/sum(Yvec.^2);

local_in1 = sum(abs((Yc+Yh)./Yvec))/no_vertex;
                                      
%oneDistanceALIs(u) = local_in1; % for constant edge flow model
normalDistanceALIs(u) = local_in1;
end
