using Gadfly

#iter is the number of iterations, passed as an argument to the function
function createHaarUnitaryMatrices() #iter::Int64)

    #for i = 1:1:iter
    Z = randn(2,2) + im * randn(2,2)

    (Q, R) = qr(Z)
    Q = Q * diagm(sign(diag(R)))
    R = inv(Q) * Z

    a = R[1,1] / norm(R[1,1])
    b = R[2,2] / norm(R[2,2])
    #c = R[3,3] / norm(R[3,3])
    #D = [a 0 0; 0 b 0; 0 0 c]
    D = [a 0; 0 b]
    
    H = Q * D
    #d = (det(H)) ^ (-1/3)
    #e[:,:,i] = d * H
    #eph[:, i] = angle(eigvals(e[:,:,i]))

    #Gadfly.plot(x=eph[:,i], y=(1/(2*pi)*ones(3,1)))
    
    #end 
end
