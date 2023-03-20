using JuMP, LinearAlgebra, GLPK, CSV, DataFrames, SparseArrays

"""
Solves the Auxilliary Problem (AP).
Returns the optimal solution x and the reduced cost RC.
"""
function solve_AP(pi::Vector{Float64}, total_length::Int, lengths::Vector{Int},demand)
    n = length(lengths)
    AP = Model(GLPK.Optimizer)
    @variable(AP, x[1:n] >= 0, Int)
    @constraint(AP, constraint_AP, sum(lengths[i] * x[i] for i in 1:n) <= total_length)
    @objective(AP, Max, sum(pi[i] * x[i] for i in 1:n))
    optimize!(AP)
    println(constraint_AP)
    x = value.(x)
    RC = 1 - objective_value(AP)
    println("obj = ", objective_value(AP))
    return x, RC
end

"""
Solves the Restricted Master Problem (RMP).
Returns the optimal dual variables pi and lambda.
"""
function solve_RMP(patterns, demand)
    p = size(patterns, 2)
    n = size(patterns,1)
    RMP = Model(GLPK.Optimizer)
    @variable(RMP, lambda[1:p] >= 0)
    @constraint(RMP, constraint_RMP[j=1:n], sum(patterns[j,i] * lambda[i] for i in 1:p) .>= demand[j])
    @objective(RMP, Min, sum(lambda[j] for j in 1:p))
    optimize!(RMP)
    println(constraint_RMP)
    pi = dual.(constraint_RMP)
    lambda = value.(lambda)
    return pi, lambda
end

"""
Print the result of the problem.
"""
function print_Result(lengths,lambda,patterns,n)
    error = 0.0001
    cpt = 0
    println("RESULTS :")
    lambda = value.(lambda)
    p = size(patterns, 2)
    println("----------------------------------------------------")
    println("Cut pattern : ",lengths)
    for j in 1:p
        res = ceil(Int, lambda[j]-error)
        if res != 0
            cpt += res
            print("We use ", res, " bar(s) with that pattern : {")
            for i in 1:(n-1)
                print(patterns[i, j], "--")
            end
            println(patterns[n, j], "}")
        end
    end
    println("----------------------------------------------------")
    println("You will need a minimum of ", cpt," raw steel bars")
    println("----------------------------------------------------")
end 

"""
Solves the given cutting stock problem.
"""
function cutting_stock()

    # Load the data from CSV file
    data = CSV.read("input.csv",DataFrame)

    # Extract the data into separate arrays
    # l_i is the desired length of steel bars.
    lengths = data[:, "lengths"]
    # n_i is the desired quantity of bars with length l_i.
    demand = data[:, "demand"]

    # total_length is the length of (raw) steel bars.
    total_length = 5600

    n = length(lengths)
    ncols = length(lengths)
    # Initialize the matrix patterns with the cutting patterns
    patterns = SparseArrays.spzeros(UInt16, n, ncols)
    for i in 1:n
        patterns[i, i] = min(floor(Int, total_length / lengths[i]), round(Int, demand[i]))
    end

    iter = 0
    RC = -1 # Objective function
    lambda = Array{Float64}
    while RC < 0
        iter += 1
        println("Loop number ", iter)

        # Solve the Restricted Master Problem (RMP)
        pi, lambda = solve_RMP(patterns, demand)
        println("lambda = ", lambda)
        println("pi = ", pi)

        # Solve the Auxilliary Problem (AP)
        x, RC = solve_AP(pi, total_length, lengths,demand)
        # Print the solution for this iteration
        println("x = ", x)
        println("RC = ",RC)

        # Update patterns matrix if necessary
        if RC < 0
            patterns = [patterns x]
        end
    end

    
    # Print the final solution
    print_Result(lengths,lambda,patterns,n)

end


cutting_stock()
    
