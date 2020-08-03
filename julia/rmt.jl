__precompile__() # Este comando es para que julia precompile el paquete

module rmt

using LsqFit: curve_fit

function P_Orthogonal(x)
    return x.*exp(-x.^2*pi/4)*pi/2
end

function P_Unitary(x)
    return (x.^2*32/pi^2).*exp(-x.^2*4/pi)
end

function unfolding(list)
    x=linspace(minimum(list),maximum(list),1000);
    y=staircase(list,x);
    model(x,p)=p[1]+p[2]*x+p[3]*x.^2+p[4]*x.^3+p[5]*x.^4+p[6]*x.^5+p[7]*x.^6+p[8]*x.^7+p[9]*x.^8+p[10]*x.^9;
    fit = curve_fit(model, x, y,ones(10));
    p=fit.param;
    return model(list,p)
end

function staircase(list,x)
    y=zeros(length(x));
    for i in range(1,length(x))
        y[i]=count(j->(j<=x[i]),list)
    end
    return y
end

end