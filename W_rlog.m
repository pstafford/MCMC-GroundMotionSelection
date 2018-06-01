function [z_wO_phi3] = W_rlog(x,n,y,m,ep,dim,Weight_input);
% this function is used to obtain potential energy with applying logarithmic distance
% the part of this function calculating the Euclidean distance between two different matrix uses 
% the code of Pairwise distance between 2 group of vectors made by Trung Duong 2011.

% Input 
% x              : selected records
% n              : number of selected records
% y              : samples of target distribution
% m              : number of samples
% ep             : Maximum PDF value
% dim            : dimentions of target intensity measure vector
% Weight_input   : Normalized weight factor used for intensity measures

% Output
% z_wO_phi3      : value of energy

modify=1/(m*ep*dim);
% diff matrix
sW = sqrt(Weight_input);
x1 = sW(ones(1,n),:).*x; 
y1 = sW(ones(1,m),:).*y; 
xx1= sum(x1.*x1,2); 
yy1= sum(y1.*y1,2)'; 
dis_c = sqrt(xx1(:,ones(1,m)) + yy1(ones(1,n),:) - 2*x1*y1');

% Same matrix
x_reduced=x;
dis_a=[];

for i=1:length(x(:,1))-1
 x_reduced(1,:)=[]; 
 x1 = sW.*x(i,:);
 y1 = sW(ones(1,length(x_reduced(:,1))),:).*x_reduced;
 xx1= sum(x1.*x1,2); 
 yy1= sum(y1.*y1,2)'; 
 dis_a_1 = sqrt(xx1(:,ones(1,length(x_reduced(:,1)))) + yy1(1,:) - 2*x1*y1');

 dis_a = [dis_a,dis_a_1];
end
% Logarithm distance
dis_a=-log(dis_a+modify);
dis_c=-log(dis_c+modify);
% Energy 
z_wO_phi3 = (1/(n*(n)))*sum(dis_a)-(1/(n*m))*sum(dis_c(:));

