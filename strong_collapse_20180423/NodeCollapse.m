%% collapse nodes

function [flag]=NodeCollapse()

global node node_coor 

sleep_node=[];
flag=0;

%% test------------set new c_flag each round
for i=1:length(node_coor)-28
    if node(i).status>0
        node(i).c_flag=1;
    end
end
%-----------------------------------------------

for i=1:length(node_coor)   % 节点i的邻居w是否能够被collapse
    if node(i).status==1 && node(i).c_flag==1
       sleep_temp=0;
       sleep_weight=0;
       for j=1:length(node(i).neighbor)
          w=node(i).neighbor(j);
          if w~=i && node(w).status==1 
              co_neighbor=intersect(node(i).neighbor,node(w).neighbor);
              diff_node=setdiff(node(i).neighbor,co_neighbor);
              if isempty(diff_node)
                  if (length(node(w).neighbor)>sleep_weight)
                      sleep_temp=w;
                      sleep_weight=length(node(w).neighbor);
                  end
              end               
          end
       end
       if sleep_temp~=0
           disp(['node ' num2str(i) ' is dominated by node ' num2str(sleep_temp)]);
           for n=1:length(node(i).neighbor)
               if node(i).neighbor(n)==i
                   continue;
               else
                   v=node(i).neighbor(n);
                   node(v).neighbor=setdiff(node(v).neighbor,i);
               end
           end
           node(i).status=0;
           node(sleep_temp).c_flag=0;
           sleep_node=[sleep_node,i];
       end
    end
end

if ~isempty(sleep_node)
    flag=1;
end