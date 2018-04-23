%% collapse nodes

function [flag]=NodeCollapse()

global node node_coor

sleep_node=[];
flag=0;

for i=1:length(node_coor)   % 节点i的邻居w是否能够被collapse
    if node(i).status==1
       sleep_temp=[];
       for j=1:length(node(i).neighbor)
          w=node(i).neighbor(j);
          if w~=i && node(w).status==1 && node(w).c_flag==1
              co_neighbor=intersect(node(i).neighbor,node(w).neighbor);
              diff_node=setdiff(node(w).neighbor,co_neighbor);
              if isempty(diff_node)
                  node(w).status=0;
                  node(i).c_flag=0;
                  disp(['node ' num2str(w) ' is dominated by node ' num2str(i)]);
                  sleep_node=[sleep_node,w];
                  sleep_temp=[sleep_temp,w];
              end               
          end
       end
       for m=1:length(sleep_temp)
           for n=1:length(node(sleep_temp(m)).neighbor)
               node_temp=node(sleep_temp(m)).neighbor(n);
               if node_temp==sleep_temp(m)
                   continue;
               else
                   node(node_temp).neighbor=setdiff(node(node_temp).neighbor,sleep_temp(m));
               end
           end
       end
    end
end

if ~isempty(sleep_node)
    flag=1;
end