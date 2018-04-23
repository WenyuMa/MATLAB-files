function[num]=CoNeighbNum(node_set)

global node

if length(node_set)<2
    num=0;
else
    co_neighb=node(node_set(1)).neighbor;
    for i=2:length(node_set)
        co_neighb=interset(co_neighb,node(node_set(i)).neighbor);
    end
    num=length(co_neighb);
end
