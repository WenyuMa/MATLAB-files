clc;
clear;
error_time=0;
sum_holes=0;
for times=1:1000
    [err,sum]=main();
    error_time=error_time+err;
    sum_holes=sum_holes+sum;
    disp(times);
    clf;
end
disp(['Error_times=' num2str(error_time)]);
disp(['Total holes=' num2str(sum_holes)]);
disp(['Error rate = ' num2str(error_time) '/' num2str(sum_holes) ' out of ' num2str(times) ' different network topology']);
