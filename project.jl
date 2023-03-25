using JuMP, LinearAlgebra, GLPK, SparseArrays, BenchmarkTools

# Path of data
FILE_NAME = "./data.txt"

"""
Import data function eg:
    i li num
    1 22 45
    2 42 38
    3 52 25
    4 53 11
    5 78 12
    100 (length of bar)
"""
function get_data(file_name::String)
    if isfile(file_name)
        f = open(file_name)
        data = readlines(f)
        line_num = length(data) - 1
        lengths = Array{Int64}(zeros(line_num))
        demands = Array{Int64}(zeros(line_num))
        for i in 1:line_num
            line = collect(split(data[i], " "))
            lengths[i] = parse(Int64, String(line[2]))
            demands[i] = parse(Int64, String(line[3]))
        end
        total_length = parse(Float64, data[length(data)])
        close(f)
        return lengths, demands, total_length
    end
    println("No data file")
    exit(1)
end

"""
Main fonction which solves the cutting stock problem
"""
function ex_cutting_stock(FILE_NAME, debug)
    lengths, demand, total_length = get_data(FILE_NAME)
    # Maximum adding generated column 
    max_gen_cols = 100
    # c = fill(1.0, length(lengths))  we consider that every c_i = 1
    n = length(lengths)
    ncols = length(lengths)

    # initial basic feasible solution
    A = SparseArrays.spzeros(UInt16, n, ncols)
    for i in 1:n
        A[i, i] = min(floor(Int, total_length / lengths[i]), round(Int, demand[i]))
    end

    # Define the RMP
    RMP = Model(GLPK.Optimizer)
    set_silent(RMP)
    @variable(RMP, x[1:ncols] >= 0)
    @objective(RMP, Min, sum(x[j] for j in 1:ncols))
    @constraint(RMP, stdemande[i=1:n], sum(A[i, j] * x[j] for j in 1:ncols) >= demand[i])
    optimize!(RMP)

    # Define the AP
    pi = dual.(stdemande)
    AP = Model(GLPK.Optimizer)
    @variable(AP, y[1:n] >= 0, Int)
    @constraint(AP, sum(y .* lengths) <= total_length)
    @objective(AP, Max, sum(y .* pi))


    while ncols - n <= max_gen_cols
        # If there is no dual problem, we cannot find the optimal solution
        if !has_duals(RMP)
            break
        end

        optimize!(AP)
        # Reduced cost : 1 - (AP), if we cannot find positive RC, it means we have already the optimal solution
        if 1.0 - sum(value.(y) .* pi) .>= 0.0
            break
        end
        newA = round.(Int, value.(y))

        # Add generated column to the RMP
        ncols += 1
        A = hcat(A, newA)
        push!(x, @variable(RMP, base_name = "x[$(ncols)]", lower_bound = 0))
        set_objective_coefficient(
            RMP,
            x[ncols],
            1 # we consider that every c_i = 1
        )
        for j in 1:n
            if newA[j] > 0
                set_normalized_coefficient(
                    stdemande[j],
                    x[ncols],
                    newA[j],
                )
            end
        end

        optimize!(RMP)

        pi = dual.(stdemande)
        # Regenerate AP
        for i in 1:n
            set_objective_coefficient(AP,x[i],pi[i])
        end

    end

    # Strain the constraint to integer value
    set_integer.(x)
    optimize!(RMP)

    # Print Result if not in benchmark mode
    if (debug)
        println("Final solution:")
        total = 0
        for i in 1:ncols
            if value(x[i]) > 0.5
                total += round(Int, value(x[i]))
                print("$(round(Int, value(x[i]))) unit(s) of pattern $(i) \t: ")
                for j in 1:n
                    curr = A[j, i]
                    if curr > 0
                        print("$(curr) of $(lengths[j]), ")
                    end
                end
                println()
            end
        end
        println("Total used : $(total) bars")
    end
    return
end


if length(ARGS) != 0
    @btime(ex_cutting_stock(FILE_NAME, false))
else
    ex_cutting_stock(FILE_NAME, true)
end
