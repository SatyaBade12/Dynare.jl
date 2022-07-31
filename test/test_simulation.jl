using Dynare

# get the solution for example1.mod in context
context = @dynare "models/simulations/example1.mod"

m = context.models[1]
r = context.results.model_results[1]
A = Matrix{Float64}(undef, m.endogenous_nbr, m.endogenous_nbr)
B = Matrix{Float64}(undef, m.endogenous_nbr, m.exogenous_nbr)
Dynare.make_A_B!(A, B, m, r)

T = 100; # horizon of simulation
x = zeros(T + 1, m.exogenous_nbr); # deterministic shocks. Period 1 is for initialization
x[2, 1:2] .= [0.01, 0.01] # e and u equal 0.01 in first simulation period and 0 afterwards
Y = Matrix{Float64}(undef, T+1, m.endogenous_nbr) # endogenous variables y c k a h b
c = r.trends.endogenous_steady_state # steady state of endogenous variables
y0 = copy(c) # endogenous variables are initialized to the steady state
Dynare.simul_first_order!(Y, y0, x, c, A, B, T)

tdf = Dynare.TimeDataFrame(Y, Dynare.get_endogenous(context.symboltable), Dynare.UndatedDate(1))

