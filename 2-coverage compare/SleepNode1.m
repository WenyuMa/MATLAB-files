% This function is used for deciding whether a node should sleep
function [sleep]=SleepNode1(w)

global node1 

sleep=0;
flag_stop=0;

if node1(w).status==1
   no_neighbor = length(node1(w).neighbour);
   if no_neighbor
       for i=1:no_neighbor-1    %choose u
         if flag_stop==0  %
           u=node1(w).neighbour(i);
           if node1(u).status==1 
               for j=i+1:no_neighbor    %choose v
                 if flag_stop==0 %
                   v=node1(w).neighbour(j);
                   if node1(v).status==1
                       index1=find(node1(u).neighbour==v);  %u,v are two neighbors
                       
                       if index1              %u and v are neighbors
                           S2=[]; 
                           cover_type=0; 
                            uv=node1(u).dist(index1); %uv 
                            uw=node1(w).dist(i); %uw
                            wv=node1(w).dist(j); %wv
                            wuv=acos((uv^2+uw^2-wv^2)/(2*uv*uw));
            
                            angle1=node1(u).angle(index1);%uv
                            angle2=node1(w).angle(i);%uw  
                            
                            if angle2>=angle1+2*wuv
                               cover_type=2;   %w覆盖了uv的2个交点  
                            elseif wuv > abs(angle2 - angle1)/2 && wuv <= (angle2 + angle1)/2
                               cover_type=1;   %w覆盖了uv的1个交点
                            end
                            
                            if cover_type==0   
                                continue;        % find next v
                            else
                                flag_cover=0;

                                for m=1:length(node1(u).neighbour)   
                                    index2=find(node1(v).neighbour==node1(u).neighbour(m)); % find all common neighbors of u and v
                                    
                                    if index2
                                        x=node1(v).neighbour(index2);  
                                        if x==w
                                            continue;
                                        elseif node1(x).status==1
                                            ux=node1(u).dist(m); %ux
                                            vx=node1(v).dist(index2); %vx
                                            xuv=acos((uv^2+ux^2-vx^2)/(2*uv*ux)); %xuv                                            
                                            angle3=node1(u).angle(m); %ux
                                            
                                            if cover_type==2
                                                if angle3>=angle1+2*xuv
                                                    flag_cover=1;   % find a node can replace u
                                                    break;
                                                elseif xuv > abs(angle3 - angle1)/2 && xuv <= (angle3 + angle1)/2
                                                    S2=[S2,x];
                                                else
                                                    continue;   %find next x
                                                end
                                            elseif cover_type==1
                                                if angle3>=angle1+2*xuv    % find a node x covers two junctions
                                                   flag_cover=1;
                                                   break;
                                                elseif xuv > abs(angle3 - angle1)/2 && xuv <= (angle3 + angle1)/2  %x covers a junction
                                                    index4=find(node1(w).neighbour==x);
                                                    
                                                    if index4        %x and w must be neighbors
                                                        wx=node1(w).dist(index4);
                                                        wux=acos((uw^2+ux^2-wx^2)/(2*uw*ux));
                                                        
                                                        if  abs(wux - xuv - wuv) < 1e-5 || abs(2*pi-wux - xuv - wuv) < 1e-5
                                                            continue;       % find next x
                                                        else 
                                                            flag_cover=1;   % find a node can replace u
                                                            break;
                                                        end 
                                                    end
                                                else
                                                    continue;
                                                end
                                            end
                                        end             
                                    end 
                                end
                                if cover_type==2 && flag_cover==0   %whether exists two nodes in S2 that cover two different junctions
                                    l=length(S2);
                                    if l<2
                                       flag_cover=0;
                                    else 
                                        flag_pair=0;
                                         for f=1:l-1
                                              if flag_pair==0
                                                  for g=f+1:l
                                                      if flag_pair==1
                                                         break;
                                                      else 
                                                          neibr=find(node1(S2(f)).neighbour==S2(g));
                                                          if neibr
                                                              index5=find(node1(u).neighbour==S2(f));%uf
                                                              index6=find(node1(u).neighbour==S2(g));%ug
                                                              index7=find(node1(S2(f)).neighbour==v);%vf
                                                              index8=find(node1(S2(g)).neighbour==v);%vg
                                                              index9=neibr;  %fg
                                 
                                                              uf=node1(u).dist(index5); %uf 
                                                              ug=node1(u).dist(index6); %ug
                                                              vf=node1(S2(f)).dist(index7); %vf 
                                                              vg=node1(S2(g)).dist(index8); %vg
                                                              fg=node1(S2(f)).dist(index9);%wx
                                 
                                                              fug=acos((uf^2+ug^2-fg^2)/(2*uf*ug));                      
                                                              vug=acos((uv^2+ug^2-vg^2)/(2*uv*ug));
                                                              fuv=acos((uf^2+uv^2-vf^2)/(2*uf*uv));
                                 
                                                              if abs(fug - vug - fuv) < 1e-5 || abs(2*pi-fug - vug - fuv) < 1e-5
                                                                  flag_pair=1;
                                                                  flag_cover=1;
                                                              end  
                                                          else  
                                                              flag_pair=1;
                                                              flag_cover=1;
                                                          end 
                                                      end 
                                                  end
                                              else
                                                  break;
                                              end
                                         end
                                    end 
                                end
                               
                                if flag_cover==1
                                    flag_stop=0; % test next junction
                                else
                                    flag_stop=1;  %find a junction that is only covered by node w
                                end
                            end
                       else
                           continue;
                       end
                   end
                 else
                     break;
                 end %
               end
           end
         else
             break;
         end %
       end
       if flag_stop==0
           sleep=1;   % node w can sleep
       else
           sleep=0;   % node w can not sleep
       end
   end
end