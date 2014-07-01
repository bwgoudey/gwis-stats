

function pair_inds = convertToPairInd(start_pair,nPair, N)
%This function should be converted to use combinatorial number system but 
%am just brute forcing it for the moment.

%%
    pair_inds=zeros(nPair,2);
    a=start_pair(1);
    b=start_pair(2);
    
    for i = 1:nPair
        b = b + 1;
        
        %If first iterator goes over total number of elements, reset
        if(b>N)            
            a=a+1;
            b=a+1;
            
            %If second iterator goes over total number of elements, reset
            if(a==N)
                pair_inds = pair_inds(1:(i-1), :);
                assert(all(pair_inds(:,1)<pair_inds(:,2)));
                return;
            end
        end
        pair_inds(i,:)=[a,b];
    end

    assert(all(pair_inds(:,1)<pair_inds(:,2)));


end

    