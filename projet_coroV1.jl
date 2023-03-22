using JuMP, LinearAlgebra, GLPK, CSV, DataFrames, SparseArrays

# path to the data file.
FILE_NAME = "./input.csv"
# total_length is the length of (raw) steel bars.
total_length = 5600

"""
Solves the Auxilliary Problem (AP).
Returns the optimal solution x and the reduced cost RC.
"""
function solve_AP(pi, total_length, lengths,demand)
    n = length(pi)
    AP = Model(GLPK.Optimizer)
    @variable(AP, x[1:n] >= 0, Int)
    @constraint(AP, constraint_AP, sum(lengths[i] * x[i] for i in 1:n) <= total_length)
    @objective(AP, Max, sum(pi[i] * x[i] for i in 1:n))
    optimize!(AP)
    x = value.(x)
    RC = 1 - objective_value(AP)
    return x, RC
end

"""
Solves the Restricted Master Problem (RMP).
Returns the optimal dual variables pi and lambda.
"""
function solve_RMP(patterns, demand, iter)
    p = size(patterns, 2)
    n = size(patterns,1)
    RMP = Model(GLPK.Optimizer)
    @variable(RMP, lambda[1:p] >= 0)
    @constraint(RMP, constraint_RMP[j=1:n], sum(patterns[j,i] * lambda[i] for i in 1:p) >= demand[j])
    @objective(RMP, Min, sum(lambda[j] for j in 1:p))
    optimize!(RMP)
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
                print(patterns[i, j], " - ")
            end
            println(patterns[n, j], "}")
        end
    end
    println("----------------------------------------------------")
    println("You will need a minimum of ", cpt," raw steel bars")
    println("----------------------------------------------------")
end 

"""
Print the number of different part and the real demand.
"""
function print_verification(length,lambda,patterns,n,ncols, demand)
    quantity = zeros(1,n)
    for j in 1:ncols
        error = 0.0001
        res = ceil(Int, lambda[j]-error)
        if res != 0
            for i in 1:n
                quantity[i] += patterns[i,j]*res
            end
        end
    end
    println("Final quantity with that cut : ", quantity)
    println("The demand was : ", demand)
    println("----------------------------------------------------")
end

"""
Solves the given cutting stock problem.
"""
function cutting_stock(file_path,total_length)
    if isfile(file_path)
        # Load the data from CSV file
        data = CSV.read("input.csv",DataFrame)

        # Extract the data into separate arrays
        # l_i is the desired length of steel bars.
        lengths = data[:, "lengths"]
        # n_i is the desired quantity of bars with length l_i.
        demand = data[:, "demand"]
    else
        println("No data file")
        exit(1)
    end

    n = length(lengths)
    ncols = length(lengths)
    # Initialize the matrix patterns with the cutting patterns
    patterns = SparseArrays.spzeros(UInt16, n, ncols)
    for i in 1:n
        patterns[i, i] = min(floor(Int, total_length / lengths[i]), round(Int, demand[i]))
    end

    iter = 0
    iterMax = 500
    RC = -1 # Objective function
    lambda = Array{Float64}
    while RC < 0 && iter < iterMax
        iter += 1

        # Solve the Restricted Master Problem (RMP)
        pi, lambda = solve_RMP(patterns, demand,iter)

        # Solve the Auxilliary Problem (AP)
        x, RC = solve_AP(pi, total_length, lengths,demand)

        # Update patterns matrix if necessary
        if RC < 0
            patterns = hcat(patterns, x)
        end
    end

    # Impose the master variables to be integer and solve.
    p = size(patterns, 2)
    n = size(patterns,1)
    RMP = Model(GLPK.Optimizer)
    @variable(RMP, lambda[1:p] >= 0,Int)
    @constraint(RMP, constraint_RMP[j=1:n], sum(patterns[j,i] * lambda[i] for i in 1:p) >= demand[j])
    @objective(RMP, Min, sum(lambda[j] for j in 1:p))
    optimize!(RMP)
    lambda = value.(lambda)

    # Print the final solution and the test
    print_Result(lengths,lambda,patterns,n)
    print_verification(length,lambda,patterns,n,p, demand)

end


cutting_stock(FILE_NAME,total_length)
