function y=map64q(x)
    tab(1:8)=[-7 -5 -1 -3 +7 +5 +1 +3];
    x_floor = floor(length(x)/6);
    temp = x(1:x_floor*6);
    temp=temp(3:3:end)*4+temp(2:3:end)*2+temp(1:3:end)+1;
    change = tab(temp);
    y = change(1:2:end)+j*change(2:2:end);
  
    
