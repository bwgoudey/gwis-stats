
function [s1]=addStructs(s1, s2)
%bwgoudey
%function [s1]=addStructs(s1, s2)
%Add 2 of the same structures together
%not much error checking or corner cases handled here

    %If structure 1 is empty return structure 2
    if(isempty(s1))
        s1=s2;
        return
    %If structure 2 is empty return structure 1
    elseif(isempty(s2))
        return
    end 
    %Otherwise, bash the results together and hope for the best!
    fns=fieldnames(s2);
    fns1=fieldnames(s1);
    for i=1:length(fns)

        fn=char(fns(i));
        
        if(~ismember(fns{i}, fns1))
            continue;
        end          
        %Skip if we have a single descriptive string/structure
        if(length(s1.(fn))==1 && length(s2.(fn))==1)
            continue;
        end
      
        s1.(fn)=[s1.(fn); s2.(fn)];
        clear s2.(fn);
    end
end