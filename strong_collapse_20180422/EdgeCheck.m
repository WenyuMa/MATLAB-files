% This function is used to check whether an edge can be dominated
function [index]=EdgeCheck(i,j)

global node edge
index1=-1;
index2=-1;
Nset1=node(edge.u(i)).neighbor;
Nset2=node(edge.v(i)).neighbor;

Nset3=union(node(edge.u(j)).neighbor,node(edge.v(j)).neighbor);
diff1=setdiff(Nset1,Nset3);
diff2=setdiff(Nset2,Nset3);

if isempty(diff1) && isempty(diff2)
    index=0;
else
    if diff1   % neighbors of u do not belong to neighbor of w\x
        for i=1:length(diff1)
            if node(diff1(i)).fence_flag==1
                diff1(i)=0;
            end
        end
        diff1=diff1(diff1~=0);
        if isempty(diff1)
            index1=0;
        elseif length(diff1)==1
            S1=node(diff1(1)).neighbor;
            if CoNeighbNum(S1)<4
                index1=1;
            end
        else
            if CoNeighbNum(diff1)==length(diff1)
                index1=2;
            else
                index1=1;
            end
        end
    end
    if diff2   % neighbors of u do not belong to neighbor of w\x
        for i=1:length(diff2)
            if node(diff2(i)).fence_flag==1
                diff2(i)=0;
            end
        end
        diff2=diff2(diff2~=0);
        if isempty(diff2)
            index2=0;
        elseif length(diff2)==1
            S2=node(diff2(1)).neighbor;
            if CoNeighbNum(S2)<4
                index2=1;
            end
        else
            if CoNeighbNum(diff2)==length(diff2)
                index2=2;
            else
                index2=1;
            end
        end
    end
    
    index=index1+index2;
end