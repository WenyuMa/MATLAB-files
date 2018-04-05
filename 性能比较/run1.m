clc;
clear;
error_time=0;
for times=1:10
    a=main();
    if a==1
        error_time=error_time+1;
    end
    disp(times);
   clf;
end
disp('error_times=');
disp(error_time);