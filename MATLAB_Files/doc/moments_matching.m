function [xi1,xi2] = moments_matching(m1,m2) % 02/02/2020 19:17
%     disp(['m1 = ',num2str(m1),' and m2 = ',num2str(m2),'.']);
    mu   = m1;
    sig2 = m2 - m1^2;
    xi1  = - ((mu+1)*(mu^2+sig2-1)) / (2*sig2);
    xi2  =   ((mu-1)*(mu^2+sig2-1)) / (2*sig2);
    
end