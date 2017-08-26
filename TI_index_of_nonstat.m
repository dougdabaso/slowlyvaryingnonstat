
function [TI_lb] = TI_index_of_nonstat(Eg,MInC,flag_theta,vector_theta)


if nargin<3
% default is not use those angles
flag_theta = 0;    
end

% Which is the index of the IMF that corresponds to the trend?
if nargin<2    
MInC=Eg(end);
end



% vector containing the angle of the scalar products between IMFs and
% residue
if nargin<4
vector_theta=zeros(1,length(Eg)-1);
end

% defaults
lower_bound = 0.00001;
upper_bound = 0.99999;
Ex = sum(Eg);
number_imfs = length(Eg);
location_factor =  0.5;
A=zeros(1,number_imfs);



% Ex is the energy of the whole signal, which is distributed among the imfs
lb = Ex.*(lower_bound).*ones(1,number_imfs); 
ub = Ex.*(upper_bound).*ones(1,number_imfs); 

% The sum of the energies of the IMFs should return EX
Aeq = ones(1,number_imfs);
beq = [Ex];

% Ex is the energy of the whole signal, which is distributed among the imfs
% Starting guess at the solution
% The energies of the imfs need take values somewhere between zero and Ex.
% with the location factor, I set where at which point at this scale the 
% starting guess will start. Default = 50%

x0 = Ex.*(location_factor).*ones(1,number_imfs);     

% Obtaining the matrix of initial conditions A
% will return a negative sign if Eg(i) > Eg(i+j)
% if Eg(i) > Eg(i+j) => A(i,:) = [A(i,1),...,A(i,j),...,A(i,length(Eg))] = [-1,...,1,...,0];
% if Eg(i) < Eg(i+j) => A(i,:) = [A(i,1),...,A(i,j),...,A(i,length(Eg))] = [1,...,-1,...,0];

current_Eg_shift = 1;
for i = 1:length(Eg)
for j = 1:length(Eg)-i
A(current_Eg_shift,i) = sign(Eg(i+j)-Eg(i));
A(current_Eg_shift,i+j) = -sign(Eg(i+j)-Eg(i));
current_Eg_shift = current_Eg_shift + 1;
end
end
b = zeros(size(A,1),1);
% [x,fval] = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub)
% [x,fval] = fmincon(@(x,MInC) 1+(sum(x(MInC:length(x)-1))/sum(x(1:MInC-1))),x0,A,b,Aeq,beq,lb,ub)

options = optimset('Algorithm','interior-point','Display','Off'); % run medium-scale algorithm

if flag_theta == 0
% do not use the angles
[x,fval,exitflag,output] = fmincon(@(x)myfun(x,MInC),x0,A,b,Aeq,beq,lb,ub,[],options);
TI_lb = fval;
else

[x,fval,exitflag,output] = fmincon(@(x)myfun(x,MInC),x0,A,b,Aeq,beq,lb,ub,[],options);    
[x_theta,fval_theta,exitflag_theta,output_theta] = fmincon(@(x)myfun_theta(x,MInC,vector_theta),x0,A,b,Aeq,beq,lb,ub,[],options);
TI_lb = [fval,fval_theta];
end



%chosen_algorithm=output.algorithm;


