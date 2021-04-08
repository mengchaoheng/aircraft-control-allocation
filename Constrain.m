function y=Constrain(x,x_min,x_max)
    [m,~]=size(x);
    y=zeros(m,1);
    for i=1:m
        if x(i)>x_max(i)
            y(i)=x_max(i);
        elseif x(i)<x_min(i)
            y(i)=x_min(i);
        else
            y(i)=x(i);
        end
    end
end


    
    
        