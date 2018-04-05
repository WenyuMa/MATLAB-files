function b=isAinB(A,B)
len_a=length(A);
for idx=1:len_a
    x=find(B==A(idx), 1);
    if isempty(x)
        b=0;
        return;
    end
end
b=1;
