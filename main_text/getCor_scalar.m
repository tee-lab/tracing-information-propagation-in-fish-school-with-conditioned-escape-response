function [cor_ut,cor_u,delay,sig_ci]=getCor_scalar(u1,u2,max_delay)
    %%%% get delay time based on correlation in speed,
    %%%% delay time is u1 relative to u2, i.e., delay>0, u1 turn first,
    %%%% verse verse
    %%%%%%%%%%% start the calculation  

    cor_ut = -max_delay:1:max_delay; % delays to calculate
    cor_u = nan(1,length(cor_ut)); % storing correlations at given lag
    T = length(u1); % length of speed data
    
    % u1 = u1 - mean(u1, 'omitmissing');
    % u2 = u2 - mean(u2, 'omitmissing');

    for k=1:length(cor_ut)
        
        cross_cor_temp = nan(1,(T-abs(cor_ut(k)))); 

        if cor_ut(k) >= 0

            for t = 1:(T-cor_ut(k))

                cross_cor_temp(t) = u1(t)*u2(t+cor_ut(k));

            end

            cor_u(k) = sum(cross_cor_temp, 'omitmissing');

        elseif cor_ut(k) <0

            for t = abs(cor_ut(k))+1:T

                cross_cor_temp(t - abs(cor_ut(k))) = u1(t)*u2(t+cor_ut(k));

            end
            
            cor_u(k) = sum(cross_cor_temp, 'omitmissing');

        end

    end

    norm_coru = sqrt(sum(u1.^2, 'omitmissing')*sum(u2.^2, 'omitmissing'));
    cor_u = cor_u/norm_coru;
    % cor_u_abs = abs(cor_u);
    % delay=cor_ut((cor_u_abs==max(cor_u_abs)));
    delay = cor_ut(cor_u == max(cor_u));
    sig_ci = 1.96/sqrt(sum(~isnan(u1)));

end