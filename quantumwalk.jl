using Gadfly
#using Reel


# We use this data structure called "Settings" to set the parameters for the calculations
type Settings
    iter
    num_steps
    lower_lim
    upper_lim
    initcondition   #choose between -1 (spin-DOWN), 0 (equal superposition), and +1 (spin-UP)
    isdemo::Bool
end


function createHaarUnitaryMatrices() #iter::Int64)
    H = zeros(2,2) #,iter)
    Z = randn(2,2) + im * randn(2,2)
    (Q, R) = qr(Z)
    Q = Q * diagm(sign(diag(R)))
    R = inv(Q) * Z

    a = R[1,1] / norm(R[1,1])
    b = R[2,2] / norm(R[2,2])
    D = [a 0; 0 b]
    H = Q * D
    return H
end



function calcprobdensity!(a_up, a_down, mysettings)

    upper_lim = mysettings.upper_lim
    lower_lim = mysettings.lower_lim
    num_steps = mysettings.num_steps
    iter = mysettings.iter
    a0 = 1.0/sqrt(2.0)  #We need this for 'equal superposition' case.

    for k = 1:1:iter

        if mysettings.initcondition == 0
            a_up[num_steps+1, 1, k] = complex(a0, 0)
            a_down[num_steps+1, 1, k] = complex(0, a0)
        elseif mysettings.initcondition == 1
            a_up[num_steps+1, 1, k] = complex(1, 0)
            a_down[num_steps+1, 1, k] = complex(0, 0)
        elseif mysettings.initcondition == -1
            a_up[num_steps+1, 1, k] = complex(0, 0)
            a_down[num_steps+1, 1, k] = complex(0, 1)
        else
            error("Please set initcondition to -1 (spin-DOWN) or 0 (equal superposition) or +1 (spin-UP)")
        end

        for t = 1:1:(num_steps-1)

            #We define 'r' to be the noise inducing parameter. Ideal case: r = a0
            if lower_lim == upper_lim
                r = a0    #Noiseless case
            else
                #Add noise in the range (lower_lim, upper_lim) using a uniform random number
                r = lower_lim + (upper_lim - lower_lim) * rand()
                #H = [r sqrt(1 - (r^2)); sqrt(1 - (r^2)) -r]
                #H = createHaarUnitaryMatrices()
            end
            
            for x = 1:1:(2*num_steps)
                if x == 1
                    #the matrices are not defined for indices less than 1, can't get a[x-1,t]
                else
                    #a_up[x,t+1,k] =  (H[1,1] * a_up[x-1,t,k]) + (H[1,2] * a_down[x-1,t,k] )
                    #a_down[x,t+1,k] =(H[2,1] * a_up[x+1,t,k] ) + (H[2,2] * a_down[x+1,t,k])
                    a_up[x,t+1,k] =  (r * a_up[x-1,t,k]) + (sqrt(1 - (r^2)) * a_down[x-1,t,k] )
                    a_down[x,t+1,k] =(sqrt(1 - (r^2)) * a_up[x+1,t,k] ) - (r * a_down[x+1,t,k])
                end                
            end
        end
    end
end


function spatialprob!(a_up, a_down, prob,mysettings)
#Probability distribution as a function of the position and time
#We plot it at t = num_steps i.e. p[:,num_steps] to get final probabilities
    num_steps = mysettings.num_steps

    for t = 1:1:(num_steps)
        for x = 1:1:(2*num_steps)+1
            prob[x,t] = prob[x,t] + (abs(a_up[x,t])^2 ) + (abs(a_down[x,t])^2 )
        end
    end

end 
 

function calcmoments!(a, firstmom, mysettings, order)
#To find the variance in position
#order is 1 for first moment, 2 for second moment... higher moments are not needed/not supported
    if order < 3
        for t = 1:1:mysettings.num_steps
            for x = 1:1:(2*mysettings.num_steps)+1
                firstmom[t] = firstmom[t] + ( ((x - mysettings.num_steps)^order)* ((abs(a[x,t])) ^2) )
                #The walk starts at the center of lattice, i.e. x = num_steps 
                #hence we use use (x - num_steps) and not x
            end
        end
    else
        return println("Higher Moments not necessary!")
    end
end

function calcmoments_tot!(prob, firstmom_tot, mysettings, order)
#To find the variance in position
#order is 1 for first moment, 2 for second moment... higher moments are not needed/not supported
    if order < 3
        for t = 1:1:mysettings.num_steps
            for x = 1:1:(2*mysettings.num_steps)+1
                firstmom_tot[t] = firstmom_tot[t] + ( ((x - mysettings.num_steps)^order)* prob[x,t] )
                #The walk starts at the center of lattice, i.e. x = num_steps 
                #hence we use use (x - num_steps) and not x
            end
        end
    else
        return println("Higher Moments not necessary!")
    end
end


function walk(mysettings::Settings)  
#This function defines the matrices that hold the simulation data
#a_up: probability density for the particle to be 'spin-UP', a_down:  probability density for the particle to be 'spin-DOWN'  
#firstmom: First moment of position, secondmom:second moment; both are defined for UP and DOWN spin as function of (x,t)
#sigmasq: Variance ("spread") as a function of time
#prob: probability of finding the particle at position x, and time t i.e. prob[x,t]
#To find the final probability vs position plot, we just need prob at t = num_steps, i.e. prob[:,num_steps]

    num_steps = mysettings.num_steps
    iter = mysettings.iter    

    a_up = zeros(Complex{Float32},((2*num_steps)+1,num_steps, iter))
    a_down = zeros(a_up)
    probdensitymatrix = typeof(a_up)
  
    firstmom_up = zeros(Complex{Float32},(num_steps,1))
    firstmom_down = zeros(firstmom_up)

    secondmom_up = zeros(Complex{Float32},(num_steps,1))
    secondmom_down = zeros(secondmom_up)
    momentmatrix = typeof(firstmom_up)

    sigmasq_up = zeros(firstmom_up)
    sigmasq_down = zeros(firstmom_down)

    firstmom_tot = zeros(firstmom_up)
    secondmom_tot = zeros(secondmom_up)
    sigmasq_tot = zeros(sigmasq_up)
      
    ratio_up = zeros(Float64, num_steps, 1)
    ratio_down = zeros(ratio_up)
    ratio_tot = zeros(ratio_up)

    prob = zeros(Float64,((2*num_steps) + 1, num_steps))
    

    if mysettings.isdemo == true

        # We don't calculate variance in position for a single run. This single run is just to show the plot of probability vs position
        mysettings.iter = 1
        calcprobdensity!(a_up, a_down, mysettings)
        spatialprob!(a_up, a_down, prob, mysettings)
        calcmoments_tot!(prob, firstmom_tot, mysettings, 1)
        calcmoments_tot!(prob, secondmom_tot, mysettings, 2)

        sigmasq_tot =  abs((secondmom_tot) - ( ((firstmom_tot).^2) ) )

        plot_probdist = Gadfly.plot( x = collect(-(num_steps):1:num_steps), y = prob[:,num_steps], Geom.line,
                                     Theme(line_width=1pt, default_color=colorant"orange"),
                                     Guide.xlabel("Position (n)"), Guide.ylabel("Probability P(n)"),
                                     Guide.title("Probability distribution as function of position")  )
        draw(PDF("probdistdemo.pdf", 6inch, 3inch), plot_probdist)
        #return plot_probdist

        plot_variance = Gadfly.plot(
                                     layer(x = 1:num_steps, y = sigmasq_tot, Geom.line,
                                           Theme(line_width=1pt, default_color=colorant"red" ),
                                           ),
                                     Guide.xlabel("time (t)"), Guide.ylabel("Variance (sigma^2)"),
                                     Guide.title("Variance of position with time ")
                                     #Guide.manual_color_key("Legend", ["spin-UP", "spin-DOWN"], ["red", "blue"])
                    )
        draw(PDF("variance.pdf", 6inch, 3inch), plot_variance)

    else
        calcprobdensity!(a_up, a_down, mysettings)
        calcmoments!((sum(a_up,3)/iter), firstmom_up, mysettings, 1)
        calcmoments!((sum(a_down,3)/iter), firstmom_down, mysettings, 1)
        calcmoments!((sum(a_up,3)/iter), secondmom_up, mysettings, 2)
        calcmoments!((sum(a_down,3)/iter), secondmom_down, mysettings, 2)
        spatialprob!((sum(a_up,3)/iter), (sum(a_down,3)/iter), prob, mysettings)

# Variance in spatial distribution for spin UP and spin DOWN respectively
        sigmasq_up = abs((secondmom_up) - ( ((firstmom_up).^2) ) )
        sigmasq_down = abs((secondmom_down) - ( ((firstmom_down).^2) ) )

        calcmoments_tot!(prob, firstmom_tot, mysettings, 1)
        calcmoments_tot!(prob, secondmom_tot, mysettings, 2)
        sigmasq_tot = abs((secondmom_tot) - ((firstmom_tot).^2) )
        
        #The "diffusion coefficient" \\sigma^{2}/t for both spin orientations
        for i = 1:num_steps
            ratio_up[i] = (sigmasq_up[i])/i
            ratio_down[i] = (sigmasq_down[i])/i
            ratio_tot[i] = (sigmasq_tot[i])/i
        end
# At any time, we are interested in the plot of one of sigma^2 vs t or (sigma^2)/t vs t. 
#Here we use plot only sigma^2 vs t, and leave the variable (sigma^2)/t i.e. ratio_up/ratio_down
    
        plot_variance = Gadfly.plot(
                                    layer(x=1:num_steps, y = sigmasq_tot, Geom.line,
                                          Theme(line_width=1pt, default_color=colorant"cyan"),
                                          ),
                                  #=  layer(x = 1:num_steps, y = sigmasq_up, Geom.line,
                                          Theme(line_width=1pt, default_color=colorant"red" ),
                                          ),
                                    layer(x = 1:num_steps, y = sigmasq_down, Geom.line,
                                          Theme(line_width=1pt, default_color=colorant"blue" ),
                                          ),
                                  =#
                                    Guide.xlabel("time t"), Guide.ylabel("Variance (sigma^2)"),
                                    Guide.title("Variance of position with time (averaged)"),
                                    #Guide.manual_color_key("Legend", ["spin-UP", "spin-DOWN"], ["red", "blue"])
                                    )


        #Plot of the probability distribution, as a function of position on the lattice
        plot_probdist = Gadfly.plot( x = collect(-(num_steps):1:num_steps), y = prob[:,num_steps], Geom.line,
                                     Theme(line_width=1pt, default_color=colorant"blue"),
                                     Guide.xlabel("Position n"), Guide.ylabel("Probability P(n)"),
                                     Guide.title("Probability of finding the particle at a lattice site (averaged)")  )
        
        pdfplots = Gadfly.draw(PDF("randomwalkplots.pdf", 6inch, 6inch), vstack(plot_variance, plot_probdist))
        
    end
end


function quantumwalk(isdemo::Bool)
#Start Here
#The 'noise'parameter can vary between 0 and 1. These are are limits of that noise.
# Set both varibles equal for noiseless case.
    lower_lim = 0 # Keep between 0 and 1
    upper_lim = 1 #Keep between 0 and 1 nd upper_lim > lower_lim

    if isdemo == true  #Only a demo
        mysettings1 = Settings(1, 100, lower_lim, upper_lim, 0, true) 
        #e.g. '1' in the last but one parameter  = particle initially has 'spin-UP'
        walk(mysettings1)
    else
        mysettings2 = Settings(100, 100, lower_lim, upper_lim, 0, false)
        walk(mysettings2)  #average over several iterations
    end
end
