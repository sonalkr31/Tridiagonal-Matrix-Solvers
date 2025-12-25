# ==========================================
# MASTER'S PROJECT: TRIDIAGONAL SOLVERS
# Methods: Thomas, LU, Cholesky
# ==========================================

using LinearAlgebra
using Plots
using BenchmarkTools
using Printf

# --- 1. Thomas Algorithm (TDMA) ---
function thomas_solve(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64})
    n = length(d)
    bc = copy(b)
    dc = copy(d)
    
    # Forward Elimination
    for i in 2:n
        m = a[i-1] / bc[i-1]
        bc[i] = bc[i] - m * c[i-1]
        dc[i] = dc[i] - m * dc[i-1]
    end
    
    # Backward Substitution
    x = zeros(n)
    x[n] = dc[n] / bc[n]
    for i in n-1:-1:1
        x[i] = (dc[i] - c[i] * x[i+1]) / bc[i]
    end
    return x
end

# --- 2. LU Decomposition ---
function lu_tridiagonal_solve(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64})
    n = length(d)
    l = zeros(n-1)
    u = zeros(n)
    
    u[1] = b[1]
    for i in 1:n-1
        l[i] = a[i] / u[i]
        u[i+1] = b[i+1] - l[i] * c[i]
    end
    
    y = zeros(n)
    y[1] = d[1]
    for i in 1:n-1
        y[i+1] = d[i+1] - l[i] * y[i]
    end
    
    x = zeros(n)
    x[n] = y[n] / u[n]
    for i in n-1:-1:1
        x[i] = (y[i] - c[i] * x[i+1]) / u[i]
    end
    return x
end

# --- 3. Cholesky Decomposition ---
function cholesky_tridiagonal_solve(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64}, d::Vector{Float64})
    n = length(d)
    p = zeros(n)
    g = zeros(n-1)
    
    p[1] = sqrt(b[1])
    for i in 1:n-1
        g[i] = a[i] / p[i]
        val = b[i+1] - g[i]^2
        if val <= 0; error("Matrix not positive definite"); end
        p[i+1] = sqrt(val)
    end
    
    y = zeros(n)
    y[1] = d[1] / p[1]
    for i in 1:n-1
        y[i+1] = (d[i+1] - g[i] * y[i]) / p[i+1]
    end
    
    x = zeros(n)
    x[n] = y[n] / p[n]
    for i in n-1:-1:1
        x[i] = (y[i] - g[i] * x[i+1]) / p[i]
    end
    return x
end

# --- Data Generation (Symmetric Positive Definite) ---
function generate_system(n)
    off_diag = rand(n-1)
    # Ensure Diagonal Dominance for Positive Definiteness
    main_diag = zeros(n)
    main_diag[1] = abs(off_diag[1]) + 2.0
    main_diag[n] = abs(off_diag[end]) + 2.0
    for i in 2:n-1
        main_diag[i] = abs(off_diag[i-1]) + abs(off_diag[i]) + 2.0
    end
    rhs = rand(n)
    return off_diag, main_diag, off_diag, rhs
end

# --- MAIN EXECUTION ---
function main()
    dimensions = [100, 1000, 5000, 10000, 20000]
    times_thomas = Float64[]
    times_lu = Float64[]
    times_cholesky = Float64[]

    println("Starting Benchmark (Please wait)...")
    println("-"^60)
    @printf "%-10s | %-12s | %-12s | %-12s\n" "N" "Thomas(s)" "LU(s)" "Cholesky(s)"
    println("-"^60)

    for n in dimensions
        a, b, c, d = generate_system(n)
        
        # Warmup and Measure Thomas
        thomas_solve(a, b, c, d) 
        t1 = @elapsed thomas_solve(a, b, c, d)
        push!(times_thomas, t1)
        
        # Measure LU
        lu_tridiagonal_solve(a, b, c, d)
        t2 = @elapsed lu_tridiagonal_solve(a, b, c, d)
        push!(times_lu, t2)
        
        # Measure Cholesky
        cholesky_tridiagonal_solve(a, b, c, d)
        t3 = @elapsed cholesky_tridiagonal_solve(a, b, c, d)
        push!(times_cholesky, t3)
        
        @printf "%-10d | %-12.6f | %-12.6f | %-12.6f\n" n t1 t2 t3
    end

    # Plotting
    println("\nGenerating Graph...")
    p = plot(dimensions, [times_thomas, times_lu, times_cholesky],
        label=["Thomas" "LU" "Cholesky"],
        title="Comparison of Tridiagonal Solvers",
        xlabel="Matrix Size (N)", ylabel="Time (s)",
        lw=3, marker=:circle, grid=true, legend=:topleft, dpi=300
    )
    
    savefig("result_graph.png")
    println("Success! Graph saved as 'result_graph.png' in your folder.")
end

# Run the main function
main()