function [nc nvc vc] = Split_Edges(nv,A)
    nc   = 0;
    [k]  = FindNextComp(nv,A);
    while k ~= 0
        nc             = nc+1;
        nvc(nc)        = 1;
        vc(nc,nvc(nc)) = k;
        av             = 1;
        while av <= nvc(nc)
            k = vc(nc,av);
            for i = 1:nv
                if A(k,i) ~= 0
                    [value] = InComp(nvc(nc),i,vc(nc,:));
                    if value == 0
                        nvc(nc)        = nvc(nc)+1;
                        vc(nc,nvc(nc)) = i;
                    end
                    A(k,i) = 0;
                    A(i,k) = 0;
                end
            end
            av = av+1;
        end
        [k]  = FindNextComp(nv,A);
    end
    return
end


function [k] = FindNextComp(n, A)
    for i = 1:n-1
        for j = i+1:n
            if A(i,j) ~= 0
                k = i;
                return
            end
        end
    end
    k = 0;
    return
end


function [value] = InComp(n, y, x)
    value = 0;
    for i = 1:n
        if x(i) == y
            value = 1;
            return
        end
    end
    return
end