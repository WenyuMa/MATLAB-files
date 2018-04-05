% This function is used for deciding whether a node should sleep
function [sleep]=SleepNode2(w)

global node2

sleep=0;
flag_stop=0;

if node2(w).status==1
   no_neighbor = length(node2(w).neighbour);
   if no_neighbor
       for i=1:no_neighbor-1    %choose u
         if flag_stop==0  %
           u=node2(w).neighbour(i);
           if node2(u).status==1 
               for j=i+1:no_neighbor    %choose v
                 if flag_stop==0 %
                   v=node2(w).neighbour(j);
                   if node2(v).status==1
                       index1=find(node2(u).neighbour==v);  %u,v are two neighbors
                       
                       if index1              %u and v are neighbors
                           S2=[]; 
                           cover_type=0; 
                            uv=node2(u).dist0(index1); %uv 
                            uw=node2(w).dist0(i); %uw
                            wv=node2(w).dist0(j); %wv
                            wuv=acos((uv^2+uw^2-wv^2)/(2*uv*uw));
            
                            angle1=node2(u).angle0(index1);%uv
                            angle2=node2(w).angle0(i);%uw  
                            
                            if angle2>=angle1+2*wuv
                               cover_type=2;   %w覆盖了uv的2个交点  
                            elseif wuv > abs(angle2 - angle1)/2 && wuv <= (angle2 + angle1)/2
                               cover_type=1;   %w覆盖了uv的1个交点
                            end
                            
                            if cover_type==0   
                                continue;        % find next v
                            else
                                flag_cover=0;

                                for m=1:length(node2(u).neighbour)   
                                    index2=find(node2(v).neighbour==node2(u).neighbour(m)); % find all common neighbors of u and v
                                    
                                    if index2
                                        x=node2(v).neighbour(index2);  
                                        if x==w
                                            continue;
                                        elseif node2(x).status==1
                                            ux=node2(u).dist0(m); %ux
                                            vx=node2(v).dist0(index2); %vx
                                            xuv=acos((uv^2+ux^2-vx^2)/(2*uv*ux)); %xuv                                            
                                            angle3=node2(u).angle0(m); %ux
                                            
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
                                                    index4=find(node2(w).neighbour==x);
                                                    
                                                    if index4        %x and w must be neighbors
                                                        wx=node2(w).dist0(index4);
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
                                                          neibr=find(node2(S2(f)).neighbour==S2(g));
                                                          if neibr
                                                              index5=find(node2(u).neighbour==S2(f));%uf
                                                              index6=find(node2(u).neighbour==S2(g));%ug
                                                              index7=find(node2(S2(f)).neighbour==v);%vf
                                                              index8=find(node2(S2(g)).neighbour==v);%vg
                                                              index9=neibr;  %fg
                                 
                                                              uf=node2(u).dist0(index5); %uf 
                                                              ug=node2(u).dist0(index6); %ug
                                                              vf=node2(S2(f)).dist0(index7); %vf 
                                                              vg=node2(S2(g)).dist0(index8); %vg
                                                              fg=node2(S2(f)).dist0(index9);%wx
                                 
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