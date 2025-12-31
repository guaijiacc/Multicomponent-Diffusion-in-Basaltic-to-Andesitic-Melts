function ww = con(data_d,beta)
    [~,m]=size(data_d);
    P = reshape(beta(29*m+1:29*m+49),7,7);
    ww = [];
    
    for i=1:m
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
        
        if data_d(i).label == "B1260"
            Beta = beta(29*m+50:29*m+56);
        elseif data_d(i).label == "B1350"
            Beta = beta(29*m+57:29*m+63);
        elseif data_d(i).label == "B1500"
            Beta = beta(29*m+64:29*m+70);
        elseif data_d(i).label == "A1250"
            Beta = beta(29*m+71:29*m+77);
        elseif data_d(i).label == "A1350"
            Beta = beta(29*m+78:29*m+84);
        elseif data_d(i).label == "A1500"
            Beta = beta(29*m+85:29*m+91);
        else
            Beta = beta(29*m+92:29*m+98);
        end

        for j=1:data_d(i).np
            x = data_d(i).x(j);
            Y = (x-x0)*exp(-Beta/2)/(sqrt(4*t));
            E = diag( erf(Y) );
            F = sqrt(t/pi)*diag(exp(-Y.*Y+Beta/2));
            ww1 = w_mean + k_mean*(x-x0) + P*E*P^-1*(w_diff/2 + k_diff/2*(x-x0)) + P*F*P^-1*k_diff; %concentration of independent components
            ww2 = 100-sum(ww1,"all"); %concentration of the dependent component,i.e.,SiO2
            ww = [ww;ww2;ww1];
        end

    end
end