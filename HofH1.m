function [h] = HofH1(H1)

C1 = 1.5501;
C2 = 0.6778;
C3 = -3.064;
C4 = 3.3;
C5 = 0.8234;
C6 = 1.1;
C7 = -1.287;

if(H1<-C6)
    if(H1<-1.6)
        h = C1*(-H1-C2)^C3+C4;
    else
        h = C5*(-H1-C6)^C7+C4;
    end
elseif(H1<=C4)
    h = 0;
elseif(H1<5.3)
    h = ((H1-C4)/C1)^(1/C3)+C2;
else
    h = ((H1-C4)/C5)^(1/C7)+C6;
end

end % endfunction
