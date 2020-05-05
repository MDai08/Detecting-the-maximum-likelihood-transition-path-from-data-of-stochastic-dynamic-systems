[a,b,c,d,e,f,g,h] = textread('par_data.txt',' %.6f %.6f %.6f %.6f %.6f %.6f  %.6f %.6f');
[m, n] = size(a);
step_len = 10;
mean_data = zeros(m/step_len,8);
for i = 10:step_len:m
    mean_data(i/step_len,1) = mean(a(i-step_len+1:i));
    mean_data(i/step_len,2) = mean(b(i-step_len+1:i));
    mean_data(i/step_len,3) = mean(c(i-step_len+1:i));
    mean_data(i/step_len,4) = mean(d(i-step_len+1:i));
    mean_data(i/step_len,5) = mean(e(i-step_len+1:i));
    mean_data(i/step_len,6) = mean(f(i-step_len+1:i));
    mean_data(i/step_len,7) = mean(g(i-step_len+1:i));
    mean_data(i/step_len,8) = mean(h(i-step_len+1:i));
end
mean_data
save mean_data.mat
