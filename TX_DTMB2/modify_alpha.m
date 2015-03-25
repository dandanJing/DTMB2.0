function alpha_out = modify_alpha(alpha_in, direction)
if direction > 0
    temp = alpha_in + 0.01;
elseif direction < 0
    temp = alpha_in - 0.01;
else
    temp = alpha_in
end
if temp > 0.05
    temp = 0.05;
elseif temp < 0.01
    temp = 0.01;
end
alpha_out = temp;