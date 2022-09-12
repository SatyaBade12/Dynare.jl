using Revise
using BAT, IntervalSets
using Dynare
using LinearAlgebra
using Distributions
using ValueShapes

context = @dynare "test/models/ls2003/ls2003.mod"

n = length(context.work.estimated_parameters)

datafile = "test/models/ls2003/data_ca1.csv"
first_obs=8
last_obs=86
iterations = 100
observations = get_observations(context, datafile, first_obs, last_obs)
ssws = SSWs(context, size(observations, 2), context.work.observed_variables)

prior = NamedTupleDist(
    psi1 = context.work.estimated_parameters.prior[1],
    psi2 = context.work.estimated_parameters.prior[2],
    psi3 = context.work.estimated_parameters.prior[3],
    rho_R = context.work.estimated_parameters.prior[4],
    alpha = context.work.estimated_parameters.prior[5],
    rr = context.work.estimated_parameters.prior[6],
    k = context.work.estimated_parameters.prior[7],
    tau = context.work.estimated_parameters.prior[8],
    rho_q = context.work.estimated_parameters.prior[9],
    rho_A = context.work.estimated_parameters.prior[10],
    rho_ys = context.work.estimated_parameters.prior[11],
    rho_pies = context.work.estimated_parameters.prior[12],
    e_Re_R = context.work.estimated_parameters.prior[13],
    e_qe_q = context.work.estimated_parameters.prior[14],
    e_Ae_A = context.work.estimated_parameters.prior[15],
    e_yse_ys =context.work.estimated_parameters.prior[16],
    e_piese_pies = context.work.estimated_parameters.prior[17]
)

function likelihood(x)
    params = []

    for i in 0:16
        #println(i, x[2*i+1])
        push!(params, x[i+1])
    end
    #println(params)
    LogDVal(Dynare.loglikelihood(params, context, observations, ssws))
end

posterior = PosteriorDensity(likelihood, prior)

samples = bat_sample(posterior, MCMCSampling(mcalg = MetropolisHastings(), 
                     nsteps = 10^2, nchains = 1)).result

