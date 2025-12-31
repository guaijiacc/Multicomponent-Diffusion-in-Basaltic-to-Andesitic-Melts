%% this program uses BFGS method to obtain eigenvectors and eigenvalues,
%% interface positions, and initial conditions

close all;
clear all;
start_time = clock();

%% ----------------------------define a structure----------------------------
data.name = "BS41&42A";
data.type= 1;% 1:diffusion couple(DC); 2:mineral dissolution(MD)
data.label= "B1260";
data.np = 100; % number of points for each exp
data.time = 1000;% exp duration
data.boundary = [];%boundary condition;size:[8,2]
data.w_mean = [];
data.delta_w = [];
data.s_mean = []; 
data.delta_s = [];
data.x0 = 0;
data.x = [];% x-coordinate;[np,1]
data.w = [];% concentration;[np,8]
data.w_calc = []; % predicted concentration;[np,8];
data.Z = []; % Z = Inv(Q)*w
data.Z_calc = []; % Z_calc = Inv(Q)*w_calc
data.x_P = [];
data.w_P = []; % predicted y;[500,8];
data.Z_P = [];
data.error=[]; % error of concentration meassured by EMPA
data.RCS = 1; % reduced chi-square for each exp

%% ------------------------------initialization---------------------------
T = [1260,1350,1500];% temp

% exp duration at different temp
%% round2 correction: assuming an average activation energy of 181 kJ for andesite
t1 = [1826.0, 1809.6, 1325.5, 1269.2, 950.8, 579.5, 740.0, 577.4, 599.9,...
    1005.53, 1065.47, 1032.73]; % exp duration at 1260
t2 = [1492.4, 1243.2, 899.0, 974.4, 563.9, 393.6, 522.0, 335.7, 346.3, 441.87,...
    949.87, 711.63, 690.89, 772.29, 402.63, 744.42, 536.93]; % at 1350
t3 = [276.1, 223.6, 247.8, 213.7, 158.5, 157.1, 154.8, 184.9, 221.4,...
    290.96, 365.99, 390.93]; % at 1500

t=[t1,t2,t3];

% num of points of each exp
Np1 = [112,121,128,140,97,146,123,128,120,...
    220, 176, 195]; % number of points of exp at 1260
Np2 = [116,84,100,109,105,101,64,122,123,141,...
    82,100,119,249,207, 211, 252];% at 1350
Np3 = [92,124,116,118,120,141,131,141,166,...
    224,243,197]; % at 1500
Np=[Np1,Np2,Np3];

Ndm = [12,17,12]; % number of exp at 1260, 1350 and 1500
threshold = [0,12,29]; % no actual meaning, just for reading data

% x0 and L of each exp
Nexp = 41; % total number of exp

name1 = ["BS1&2C";"BS3&4C";"BS5&6C";"BS7&8C";"BS9&10C";"BS11&12C";"BS13&14C";"BS17&18C";"BS19&20C";...
    "AND7&8C2";"AND9&10C1";"AND11&12C1"];
name2 = ["BS1&2A";"BS3&4A";"BS5&6A";"BS7&8A";"BS9&10A";"BS11&12A";"BS13&14A";"BS17&18A";"BS19&20A";...
    "BS43&44A";"AND1&2A1";"AND3&4A1";"AND5&6A1";"AND7&8A3";"AND9&10A1";"AND11&12A1"; "BS14&AND3A1"];
name3 = ["BS1&2B";"BS3&4B";"BS5&6B";"BS7&8B";"BS9&10B";"BS11&12B";"BS13&14B";"BS17&18B";"BS19&20B";...
    "AND7&8B1"; "AND9&10B1";"AND11&12B1"];
name = [name1;name2;name3];


Label1 = ["B1260";"B1260";"B1260";"B1260";"B1260";"B1260";"B1260";"B1260";"B1260";...
    "A1250";"A1250";"A1250"];
Label2 = ["B1350";"B1350";"B1350";"B1350";"B1350";"B1350";"B1350";"B1350";"B1350";...
    "B1350";"A1350";"A1350";"A1350";"A1350";"A1350";"A1350";"B&A1350"];
Label3 = ["B1500";"B1500";"B1500";"B1500";"B1500";"B1500";"B1500";"B1500";"B1500";...
    "A1500";"A1500";"A1500"];
Label = [Label1; Label2; Label3];

filename = "Initial.xlsx";
% filename = "fitted_BA_AND.xlsx";
realboundary = xlsread(filename,"boundary");
realboundary = realboundary(1:7,:);
realinterface = xlsread(filename,"interface");
realinterface = realinterface(:,1);
% read initial beta
read_b = xlsread(filename,"beta");
read_b = read_b(1:14,1:7);

%read data
for i=1:3
    T_str = string(T(i));
    realdata = xlsread(strcat(T_str,"_RealData_8Comp.xlsx"));
    realerror = xlsread(strcat(T_str,"_RealError_8Comp.xlsx"));
    
    for j=1:Ndm(i)
        j1 = j+threshold(i);
        data(j1).name = name(j1);
        data(j1).type = 1;
        data(j1).label = Label(j1);
        data(j1).np= Np(j1);
        data(j1).time= t(j1);
        data(j1).w_mean = realboundary(:,4*j1-3);
        data(j1).delta_w = realboundary(:,4*j1-2);
        data(j1).s_mean = realboundary(:,4*j1-1);
        data(j1).delta_s = realboundary(:,4*j1);
        data(j1).x0 = realinterface(j1);
        data(j1).x = realdata(1:Np(j1),12*j-11);
        data(j1).x_P = linspace(min(min(data(j1).x),-1500),max(max(data(j1).x),1500),500)';
        data(j1).w = realdata(1:Np(j1),(12*j-10):(12*j-3));
        data(j1).error=realerror(:,j);
    end
end

% read penalty factor
penalty_factor = xlsread("Penalty_factor.xlsx");

%% ---------------------Choose data to fit--------------------------------
data_d = data;  % BA + AND + BS14&AND3A1
factor = penalty_factor;

[~,m]=size(data_d); % number of exp involved in fitting

Ni = 7; %7 independent components; 

% initialize para
para = [];
for i = 1:m
    para = [para; data_d(i).x0; data_d(i).w_mean; data_d(i).delta_w; data_d(i).s_mean; data_d(i).delta_s];
end

%% ----------------- claculate degree of freedom-------------------------
Npoint = 0;
for i = 1:m
    Npoint = Npoint + data_d(i).np;
end


N_eig = 7;
DF = 8*Npoint - (29*m + 42 + Ni*N_eig); % three temperature included

%% -------------------------initialize error----------------------------
error =[];
for i = 1:m
    Np_i = data_d(i).np;
    error_i = kron(ones(Np_i,1),data_d(i).error);
    error = [error;error_i];
end

Q = read_b(N_eig + 1:end,:);% matrix of eigenvector of D

%% ----------------------start of BFGS method-------------------------------
% BFGS method was discovered and published independently by  
% Broyden, Fletcher, Goldfarb, and Shanno in 1970. And it is the most
% famous and widely used quasi-Newton method. Part of variables defined in
% the following part follows LI & FUKUSHIMA, 2001. 

% x, y1, y2 are used for plotting
x = [];
y1 = [];
y2 = [];

n_loop1= 0; % number of outer loops 
n_loop2= 0; % number of inner loops

broadcast_factor = kron(ones(7,1),factor')*5E9; % penalty factor
Lambda_condition = reshape([zeros(15,m);broadcast_factor;broadcast_factor/4],[29*m,1]);
Lambda_beta = zeros((N_eig+7)*7,1);
Lambda = diag([Lambda_condition;Lambda_beta]);
beta =[para;reshape(Q,Ni^2,1);reshape(read_b(1:N_eig,:)',[],1)];% parameters to be fitted;
[n,~] = size(beta);

w0=[]; % initial concentration
for i=1:m
    transform_w = reshape((data_d(i).w)',[],1);
    w0 = [w0;transform_w];
end
w = con(data_d,beta); %calculated concentration
r = (w0-w)./error; % error
S = r'*r + beta'*Lambda*beta; % penalized sum of square error
S_DF = S/DF;   % reduced chi-sqaure
J = calculate_J_par(data_d,beta); % calculate Jacobian
g = 2*J'*r + 2*Lambda*beta; %gradient of S (LI & FUKUSHIMA, 2001)
H = eye(29*m+49+7*N_eig); % H matrix (LI & FUKUSHIMA, 2001)

x = [x;n_loop1];
y1 = [y1;S_DF];
y2 = [y2;n_loop2]; 
yyaxis left
plot(x,y1,'-b');
xlabel('num of loop1');
ylabel('reduced chi-square');
yyaxis right
plot(x,y2,'-*r');
ylabel('num of loop2');
drawnow;

max1 = 3000;
% max1 = 0;
max2 = 100;
sigma = 0.5; % parameter used in Armijo-type line search (LI & FUKUSHIMA, 2001)
rho = 0.5; % parameter used in Armijo-type line search (LI & FUKUSHIMA, 2001)
n_loop1 = 1;
while n_loop1 < max1 && norm(g)/DF/2 > 1E-5
% when n_loop1 exceeds max1 or ||S_DF|| > 1E-5, iteration ends

    p = -H*g; 
    % the BFGS direction obtained by solving the linear equation: p=-H*g
    step=1; 
    n_loop2 = 0;

    % search for proper step using Armijo-type line search.
    while n_loop2 < max2
        beta1 = beta + step*p;
        w1 = con(data_d,beta1);
        r1 = (w0-w1)./error;
        S1 = r1'*r1 + beta1'*Lambda*beta1;
        l= S + sigma*step*g'*p; % parameter used in Armijo-type line search (LI & FUKUSHIMA, 2001)

        % if S1<l, the step is OK, otherwise, reduce the step by a factor of rho
        if S1 < l
            break;
        end
        step = rho*step; 
        n_loop2 = n_loop2+1;
    end

    if n_loop2==max2
        fprintf("n_loop2 = %d \n",max2);
        break;
    end

    s = beta1 - beta; % increase in beta (LI & FUKUSHIMA, 2001)
    % update beta
    beta = beta1;
    % normalize eigenvectors
    Q = reshape(beta(29*m+1:29*m+49),7,7);
    [~, maxIndices] = max(abs(Q)); % find the indices of the maximum in each column vector of Q
    for i =1:7
        Max = beta(29*m+7*i-7+maxIndices(i));
        beta(29*m+7*i-6:29*m+7*i) = beta(29*m+7*i-6:29*m+7*i) / norm(beta(29*m+7*i-6:29*m+7*i)) * ( Max/abs(Max));
    end

    % update w,r,S,S_DF,J and y
    w = w1;
    r = r1;
    S = S1;
    S_DF = S/DF;
    J = calculate_J_par(data_d,beta);
    y = 2*J'*r + 2*Lambda*beta - g; 

    % update H matrix to guarantee the positive definiteness of H
    if norm(g)>1
        if (y'*s)/(s'*s)>= 10^-6*norm(g)^0.01
            H = H + (y'*s + y'*H*y)*(s*s')/(y'*s)^2 - (H*y*s'+s*y'*H)/(y'*s);
        end
    else
        if (y'*s)/(s'*s)>= 10^-6*norm(g)^3
            H = H + (y'*s + y'*H*y)*(s*s')/(y'*s)^2 - (H*y*s'+s*y'*H)/(y'*s);
        end
    end

    % update g
    g = 2*J'*r + 2*Lambda*beta;

    x = [x;n_loop1];
    y1 = [y1;S_DF];
    y2 = [y2;n_loop2]; 
    yyaxis left
    plot(x,y1,'-b');
    xlabel('num of loop1');
    ylabel('reduced chi-square');
    yyaxis right
    plot(x,y2,'-*r');
    ylabel('num of loop2');
    drawnow;

    n_loop1 = n_loop1 +1;
end

% normalize eigenvectors
Q = reshape(beta(29*m+1:29*m+49),7,7);
[~, maxIndices] = max(abs(Q)); % find the indices of the maximum in each column vector of Q
for i =1:7
    Max = beta(29*m+7*i-7+maxIndices(i));
    beta(29*m+7*i-6:29*m+7*i) = beta(29*m+7*i-6:29*m+7*i) / norm(beta(29*m+7*i-6:29*m+7*i)) * ( Max/abs(Max));
end

%% ------------------------end of BFGS method---------------------------------

transform_w = reshape(w,8,[])';
j=0;

%% store calculated concentration
for i=1:m
    data_d(i).w_calc = transform_w(j+1:j+data_d(i).np,:);
    j = j + data_d(i).np;
end

%% calculate reduced chi-squre for each experiment
for i=1:m
    DF_i = data_d(i).np * 8 ; % degree of freedom
    
    error_i = kron(ones(data_d(i).np,1),data_d(i).error);
    w0_i = reshape((data_d(i).w)',[],1);
    w_i = reshape((data_d(i).w_calc)',[],1);
    r_i = (w0_i - w_i)./error_i;
    data_d(i).RCS = r_i'*r_i/DF_i;
end

%% calculate eigenvalues
%trans_beta is a N_eig+7 by 7 matrix. 
% The last seven rows: matrix of eigenvector

trans_beta = zeros(N_eig+7,7);  % change dimensions of beta
trans_beta(1:N_eig,:) = reshape(beta(29*m+50:29*m+49+7*N_eig) ,7,N_eig)'; % calculate lamda of each exp
trans_beta(N_eig+1:N_eig+7,:) = reshape(beta(29*m+1:29*m+49),7,7);

%% sequence eigenvalues using bubble sort
for i = 1:6
    for j = 1:7-i
        if trans_beta(1,j)>trans_beta(1,j+1)
            store = trans_beta(:,j+1);
            trans_beta(:,j+1)=trans_beta(:,j);
            trans_beta(:,j)= store;
        end
    end
end

%% calculate error of trans_beta
beta(29*m+1:29*m+49) = reshape(trans_beta(N_eig+1:N_eig+7,:),49,1);
beta(29*m+50:29*m+49+7*N_eig)= reshape(trans_beta(1:N_eig,:)',7*N_eig,1);
J7 = calculate_J7_par(data_d,beta);

Lambda_condition = reshape([zeros(15,m);broadcast_factor;broadcast_factor/4],[29*m,1]);
Lambda_beta_prime = zeros((N_eig+6)*7,1);
Lambda_prime = diag([Lambda_condition;Lambda_beta_prime]);
J_prime = [J7; sqrt(Lambda_prime)];

[U,W,V] = svd(J_prime, 'econ');
inv_A = V/(W'*W)*V';
inv_A = (diag(inv_A)).^0.5;

error_beta = zeros(29*m+49+7*N_eig,1);
for i = 0:6
    error_beta(29*m+7*i+1:29*m+7*i+6)= inv_A(29*m+6*i+1:29*m+6*i+6);
end
error_beta(29*m+50:29*m+49+7*N_eig) = inv_A(29*m+43:29*m+42+7*N_eig);

error_trans_beta = zeros(N_eig+7,7);
error_trans_beta(1:N_eig,:) = reshape( error_beta(29*m+50:29*m+49+7*N_eig) , 7,N_eig)';
error_trans_beta(N_eig+1:N_eig+7,:) = reshape( error_beta(29*m+1:29*m+49), 7,7);

para = reshape(beta(1:29*m),[29,m]);
X0 = para(1,:)';
boundary = reshape(para(2:end,:),[7,numel(para(2:end,:))/7]);

error_para = reshape(inv_A(1:29*m),[29,m]);
error_X0 = error_para(1,:)';
error_boundary = reshape(error_para(2:end,:),[7,numel(error_para(2:end,:))/7]);

%% write fitting parameters into files
writematrix([boundary; error_boundary], "fitted_BA_AND.xlsx", "Sheet", "boundary");
writematrix([X0, error_X0], "fitted_BA_AND.xlsx", "Sheet", "interface");
writematrix([trans_beta; error_trans_beta], "fitted_BA_AND.xlsx", "Sheet", "beta");

%% using obtained beta to predict concentration profiles of all exp
w_predict = prediction(data_d,beta);
transform_w = reshape(w_predict,8,[])';
j=0;
for i=1:m
    [np,~] = size(data_d(i).x_P);
    data_d(i).w_P = transform_w(j+1:j+np,:);
    j = j + np;
end

%% calculate Z,Z_P
Q = reshape(beta(29*m+1:29*m+49),7,7);
for i = 1:m
    data_d(i).Z = (inv(Q)*data_d(i).w(:,2:8)')';
    data_d(i).Z_calc = (inv(Q)*data_d(i).w_calc(:,2:8)')';
    data_d(i).Z_P = (inv(Q)*data_d(i).w_P(:,2:8)')';
end

elapsed_time = etime(clock(),start_time); 