% calculate and return the derivates of r with respect to P (matrix of eigenvector) and beta

function J=calculate_J_par(data_d,beta)
    [~,m]=size(data_d); %number of experiments included
    P = reshape(beta(29*m+1:29*m+49),7,7); % matrix of eigenvector
    N_eig = floor(size(beta,1) - 29*m - 49)/7;

    num = 0; 
    n_start = 0;
    n_end = 0;
    arr.JT = [];
    
    parfor i=1:m
        % w1 = data_d(i).boundary(2:8,1);
        % w2 = data_d(i).boundary(2:8,2);
        %for DC, w1 = w_left, w2 = w_right
        %for MD, w1 = w_interface, w2 = w_initial
        
        t = data_d(i).time;
        x0 = beta(29*i-28);
        w_mean = beta(29*i-27:29*i-21);
        w_diff = beta(29*i-20:29*i-14);
        k_mean = beta(29*i-13:29*i-7);
        k_diff = beta(29*i-6:29*i);
        
        error1 = data_d(i).error(2:8);
        error2 = data_d(i).error(1);
        
        if data_d(i).label == "B1260"
            num = 0;
            Beta = beta(29*m+50:29*m+56);
        elseif data_d(i).label == "B1350"
            num = 1;
            Beta = beta(29*m+57:29*m+63);
        elseif data_d(i).label == "B1500"
            num = 2;
            Beta = beta(29*m+64:29*m+70);
        elseif data_d(i).label == "A1250"
            num = 3;
            Beta = beta(29*m+71:29*m+77);
        elseif data_d(i).label == "A1350"
            num = 4;
            Beta = beta(29*m+78:29*m+84);
        elseif data_d(i).label == "A1500"
            num = 5;
            Beta = beta(29*m+85:29*m+91); 
        else
            num = 6;
            Beta = beta(29*m+92:29*m+98); 
        end
        
        n_start = num*7+1;
        n_end = (num+1)*7;
        
        arr(i).JT = zeros(size(beta,1), 8*data_d(i).np);
        for j=1:data_d(i).np
            x= data_d(i).x(j);
            Y= (x-x0)*exp(-Beta/2)/(sqrt(4*t)); % Eq.(6) in "derivation1.pdf"
            E = diag( erf(Y) ); %  Eq.(7) in "derivation1.pdf"
            F = sqrt(t/pi)*diag(exp(-Y.*Y+Beta/2));
            dE1 = -diag(exp(-Y.*Y - Beta/2))/sqrt(pi*t);
            dE2 = -diag(exp(-Y.*Y - Beta/2))*(x-x0)/sqrt(4*pi*t);
            dF1 = -dE2;
            dF2 = diag(exp(-Y.*Y + Beta/2)).*(diag(Y.*exp(-Beta/2)*(x-x0)/sqrt(4*t)) + eye(7)/2)*sqrt(t/pi);
            
            wx0 = (-k_mean + P*dE1*P^-1*(w_diff/2+k_diff/2*(x-x0)) - P*E*P^-1*k_diff/2 + P*dF1*P^-1*k_diff)'; % Eq.(18) in "derivation3.pdf"
            % wp is the partial derivative of concentration of 
            % independent components with respect to x0

            ww_mean = eye(7);

            ww_diff = (P*E*P^-1)'/2;

            wk_mean = (x-x0)*eye(7);

            wk_diff = (P*E*P^-1*(x-x0)/2 + P*F*P^-1)';

            wp = kron(E*P^-1*(w_diff/2+k_diff/2*(x-x0)),eye(7))...
                - kron(P^-1*(w_diff/2+k_diff/2*(x-x0)),(P*E*P^-1)')...
                + kron(F*P^-1*k_diff, eye(7))...
                - kron(P^-1*k_diff, (P*F*P^-1)');
            % Eq.(19) in "derivation3.pdf"
            % wp is the partial derivative of concentration of 
            % independent components with respect to P
            
            wBeta = ( kron((P^-1*(w_diff/2+k_diff/2*(x-x0)))',P)* diag(reshape(dE2,7^2,1))...
                + kron((P^-1*k_diff)',P)*diag(reshape(dF2,7^2,1)) )';

            % Eq.(20) in "derivation3.pdf"
            % wBeta is the partial derivative of concentration of
            % independent components with respect to beta
        
            rx01 =  -wx0/diag(error1);% Eq.(24) in "derivation3.pdf"
            % partial derivative of residue error of independent
            % components with respect with to x0

            rw_mean1 =  -ww_mean/diag(error1);

            rw_diff1 =  -ww_diff/diag(error1);

            rk_mean1 =  -wk_mean/diag(error1);

            rk_diff1 =  -wk_diff/diag(error1);

            rp1 = -wp/diag(error1); % Eq.(25) in "derivation3.pdf"
            % partial derivative of residue error of independent
            % components with respect with to P
            
            rBeta1= -wBeta/diag(error1); % Eq.(26) in "derivation3.pdf"
            % partial derivative of residue error of independent
            % components with respect with to beta

            
            rx02 = sum(wx0,2)/error2; % Eq.(27) in "derivation3.pdf"
            % partial derivative of residue error of the dependent
            % components with respect with to x0

            rw_mean2 = sum(ww_mean,2)/error2;

            rw_diff2 = sum(ww_diff,2)/error2;

            rk_mean2 = sum(wk_mean,2)/error2;

            rk_diff2 = sum(wk_diff,2)/error2;

            rp2 = sum(wp,2)/error2; % Eq.(28) in "derivation3.pdf"
            % partial derivative of residue error of the dependent
            % components with respect with to P
            
            rBeta2 = sum(wBeta,2)/error2; % Eq.(29) in "derivation3.pdf"
            % partial derivative of residue error of the dependent
            % components with respect with to beta 
            
            rx0 = [rx02,rx01];
            % partial derivative of residue error of all
            % components with respect with to x0

            rw_mean = [rw_mean2,rw_mean1];

            rw_diff = [rw_diff2,rw_diff1];

            rk_mean = [rk_mean2,rk_mean1];

            rk_diff = [rk_diff2,rk_diff1];

            rp = [rp2,rp1];
            % partial derivative of residue error of all
            % components with respect with to P
        
            rBeta3=[rBeta2,rBeta1];% 7^2 * (7+1)
            % partial derivative of residue error of all
            % components with respect with to beta

            r_exp = zeros(29*m,8);
            r_exp(29*i-28:29*i,:) = [rx0;rw_mean;rw_diff;rk_mean;rk_diff];
            
            rBeta=zeros(7*N_eig,8);
            for k=1:8
                rBeta(n_start:n_end,k)= diag(reshape(rBeta3(:,k),7,7));
            end
            arr(i).JT(:,8*j-7:8*j) =[r_exp;rp;rBeta];
        end
            
    end

    JT = [];
    for i = 1:m
        JT = [JT, arr(i).JT];
    end
    J = JT';
end