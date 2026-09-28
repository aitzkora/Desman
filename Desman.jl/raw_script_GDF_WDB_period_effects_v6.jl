using Desman
using DataFrames
using CSV
using Printf
using LinearAlgebra
using Statistics
using Serialization

# =========================================================
# GDF-WDB : effet de la période d'échantillonnage
# Version v6 : ajout de M0_session pour isoler structure sessionnelle vs effet période
#
# Modèles comparés pour chaque base covariable : null et WidthWetBed
#   M0  : poolé sans période, optimizeINRIA de référence, groupes siteID
#   M0_session : poolé sans effet période, mais groupes siteID × période
#   M1a : poolé/session avec effet période sur l'échelle Weibull
#   M1b : M1a + effet période sur la shape Weibull beta
#   M1c : M1b + effet période sur lambda de frailty Gamma
#   M2  : modèles séparés hiver/été, optimizeINRIA de référence
#
# M0_session et M1a/M1b/M1c utilisent des groupes siteID × période, pour respecter
# l'indépendance des sessions d'échantillonnage.
# Effets canoniques : hiver = référence, été = 1.
# =========================================================

function _include_period_effects_likelihood_v6()
    candidates = [
        joinpath(@__DIR__, "likelihood_GDF_WDB_period_effects_v6.jl"),
        joinpath(@__DIR__, "src", "likelihood_GDF_WDB_period_effects_v6.jl"),
        joinpath(pwd(), "likelihood_GDF_WDB_period_effects_v6.jl"),
        joinpath(pwd(), "src", "likelihood_GDF_WDB_period_effects_v6.jl")
    ]
    for f in candidates
        if isfile(f)
            include(f)
            return f
        end
    end
    error("Cannot find likelihood_GDF_WDB_period_effects_v6.jl. Place it next to this raw script or in src/.")
end

const PERIOD_EFFECTS_LIKELIHOOD_FILE_V6 = _include_period_effects_likelihood_v6()
using .GDF_WDB_Period_Effects_V6

if !@isdefined(datapath)
    local cands = [
        joinpath(@__DIR__, "data"),
        joinpath(@__DIR__, "..", "data"),
        joinpath(@__DIR__, "src", "..", "data"),
        joinpath(pwd(), "data")
    ]
    local existing = filter(isdir, cands)
    datapath = isempty(existing) ? joinpath(@__DIR__, "data") : first(existing)
end

# ----------------------------
# Helpers généraux
# ----------------------------
compute_AIC(f::Float64, K::Integer) = isfinite(f) ? 2.0 * f + 2.0 * Int(K) : Inf
gradient_ok(g::Vector{Float64}, threshold::Float64) = isfinite(norm(g)) && norm(g) < threshold
clamp_delta(x::Real; bound::Float64=3.0) = Float64(min(max(x, -bound), bound))

function safe_durationg!(df::DataFrame)
    df.durationg = (x -> begin
        if x isa Missing
            return Inf
        elseif x isa AbstractString
            return (x == "NA" || isempty(x)) ? Inf : parse(Float64, replace(x, "," => "."))
        elseif isfinite(Float64(x))
            return Float64(x)
        else
            return Inf
        end
    end).(df[:, :durationg])
    df.durationd = Float64.(df[:, :durationd])
    return df
end

function read_occurrence_data_with_period(datapath::AbstractString)
    df_w = CSV.read(joinpath(datapath, "survcalibfeb2024_1.csv"), DataFrame; delim=';', decimal=',')
    df_s = CSV.read(joinpath(datapath, "survcalibaug2024_1.csv"), DataFrame; delim=';', decimal=',')

    safe_durationg!(df_w)
    safe_durationg!(df_s)

    df_w.period_label .= "hivernale"
    df_s.period_label .= "estivale"
    df_w.period_estivale .= 0.0
    df_s.period_estivale .= 1.0

    df = vcat(df_w, df_s, cols=:union)
    df.siteID_period = string.(df.siteID) .* "__" .* string.(df.period_label)

    Σ = CSV.read(joinpath(datapath, "Latrine_cov.csv"), DataFrame; delim=';', decimal=',')
    df_joined = outerjoin(df, Σ; on=:siteID, matchmissing=:equal, makeunique=true)
    return df_joined, Σ
end

function initial_lambda_from_df(df::DataFrame)
    finite_idx = findall(df.durationg .< Inf)
    isempty(finite_idx) && error("No finite durationg values available for initialization.")
    λ0 = sum((df.durationg[finite_idx] .+ df.durationd[finite_idx]) ./ 2) / length(finite_idx)
    return [Float64(λ0), 1.0, 1.0]
end

function covariate_indices(cov_names::Vector{String}, wanted::Vector{String})
    idx = Int64[]
    for nm in wanted
        pos = findfirst(==(nm), cov_names)
        pos === nothing && error("Covariate '$nm' not found. Available: $(join(cov_names, ", "))")
        push!(idx, Int64(pos))
    end
    return idx
end

function fit_reference_GDF_WDB(bio, selVar::Vector{Int64}, λ_init::Vector{Float64};
                               gradnorm_threshold::Float64=0.05,
                               print_iter::Bool=false)
    f, sol, g, it, sim = optimizeINRIA(bio, selVar, λ_init; print_iter=print_iter)
    return (
        f = Float64(f),
        sol = Float64.(sol),
        g = Float64.(g),
        gradnorm = norm(g),
        gradient_ok = gradient_ok(Float64.(g), gradnorm_threshold),
        it = Int(it),
        sim = Int(sim),
        status = (isfinite(f) && all(isfinite, sol) ? "ok" : "failed")
    )
end

function optimize_M1_effects_v6(bio_occ::GDF_WDB_Period_Effects_V6.OccurrenceBiotope,
                                selVar::Vector{Int64},
                                effect_structure::Symbol,
                                x0_without_period::Vector{Float64},
                                initial_deltas::Vector{Float64};
                                period_effect_bound::Float64=3.0,
                                gradnorm_threshold::Float64=0.05,
                                print_iter::Bool=false)
    neff = GDF_WDB_Period_Effects_V6.n_period_effects(effect_structure)
    x0 = [copy(x0_without_period); initial_deltas[1:neff]]
    nbase = length(x0_without_period)

    f = GDF_WDB_Period_Effects_V6.logLikelihood(bio_occ, selVar, effect_structure)
    g! = GDF_WDB_Period_Effects_V6.getg!(bio_occ, selVar, effect_structure)

    # alpha/beta/lambda/covariables gardent les bornes de la référence.
    # Les effets période sont signés.
    lb = [fill(1e-6, nbase); fill(-period_effect_bound, neff)]
    ub = [fill(Inf, nbase); fill(period_effect_bound, neff)]

    fout, xout, it, sim = Desman.bfgsb(f, g!, x0, lb, ub; print_iter=print_iter, ϵ=5e-5)
    gout = zeros(length(xout))
    g!(gout, xout)
    return (
        f = Float64(fout),
        sol = Float64.(xout),
        g = Float64.(gout),
        gradnorm = norm(gout),
        gradient_ok = gradient_ok(gout, gradnorm_threshold),
        it = Int(it),
        sim = Int(sim),
        status = (isfinite(fout) && all(isfinite, xout) ? "ok" : "failed")
    )
end

function base_scale_values_reference(bio, selVar::Vector{Int64}, sol::Vector{Float64})
    vals = zeros(Float64, length(bio.durd))
    for j in eachindex(vals)
        sc = sol[1]
        for (k, cov_idx) in enumerate(selVar)
            sc += bio.matCov[j, cov_idx] * sol[3+k]
        end
        vals[j] = sc
    end
    return vals
end

function weibull_quantiles(scales::Vector{Float64}, betas::Vector{Float64}; probs=[0.25, 0.5, 0.75, 0.9])
    out = Dict{Float64, Vector{Float64}}()
    for p in probs
        out[p] = scales .* ((-log(1.0 - p)) .^ (1.0 ./ betas))
    end
    return out
end

function quantile_summary_rows(base_model::String, model::String, period::String, cov_label::String,
                               scales::Vector{Float64}, betas::Vector{Float64};
                               effect_structure::String="none")
    qs = weibull_quantiles(scales, betas)
    rows = NamedTuple[]
    for p in sort(collect(keys(qs)))
        vals = qs[p]
        push!(rows, (
            base_model = base_model,
            model = model,
            effect_structure = effect_structure,
            period = period,
            covariates = cov_label,
            probability = p,
            mean_time = mean(vals),
            median_time = median(vals),
            min_time = minimum(vals),
            max_time = maximum(vals)
        ))
    end
    return rows
end

function summarize_reference_fit(label::String, period::String, bio, selVar::Vector{Int64}, fit)
    scales = base_scale_values_reference(bio, selVar, fit.sol)
    betas = fill(fit.sol[2], length(scales))
    medt = weibull_quantiles(scales, betas; probs=[0.5])[0.5]
    return (
        model = label,
        period = period,
        n_occurrences = length(bio.durd),
        n_latrines = bio.N,
        negloglik = fit.f,
        K = length(fit.sol),
        AIC_total = compute_AIC(fit.f, length(fit.sol)),
        gradnorm = fit.gradnorm,
        gradient_ok = fit.gradient_ok,
        niter = fit.it,
        nsim = fit.sim,
        alpha = fit.sol[1],
        beta = fit.sol[2],
        lambda = fit.sol[3],
        frailty_variance = 1.0 / fit.sol[3],
        mean_scale = mean(scales),
        median_scale = median(scales),
        mean_conditional_median_time = mean(medt),
        median_conditional_median_time = median(medt)
    )
end

function summarize_M1_effects_fit(label::String, effect_structure::Symbol,
                                  bio_occ::GDF_WDB_Period_Effects_V6.OccurrenceBiotope,
                                  selVar::Vector{Int64}, fit)
    scales = GDF_WDB_Period_Effects_V6.period_adjusted_scale_values(bio_occ, selVar, fit.sol, effect_structure)
    betas = GDF_WDB_Period_Effects_V6.period_adjusted_beta_values(bio_occ, fit.sol, selVar, effect_structure)
    lambdas = GDF_WDB_Period_Effects_V6.group_lambda_values(bio_occ, fit.sol, selVar, effect_structure)
    medt = weibull_quantiles(scales, betas; probs=[0.5])[0.5]
    didx = GDF_WDB_Period_Effects_V6.delta_indices(selVar, effect_structure)

    δs = didx.δs == 0 ? 0.0 : fit.sol[didx.δs]
    δb = didx.δb == 0 ? 0.0 : fit.sol[didx.δb]
    δl = didx.δl == 0 ? 0.0 : fit.sol[didx.δl]

    return (
        model = label,
        period = effect_structure === :none ? "pooled_session_no_period" : "pooled_with_period",
        n_occurrences = length(bio_occ.durd),
        n_latrines = bio_occ.N,
        negloglik = fit.f,
        K = length(fit.sol),
        AIC_total = compute_AIC(fit.f, length(fit.sol)),
        gradnorm = fit.gradnorm,
        gradient_ok = fit.gradient_ok,
        niter = fit.it,
        nsim = fit.sim,
        alpha = fit.sol[1],
        beta = fit.sol[2],
        lambda = fit.sol[3],
        frailty_variance_winter = 1.0 / fit.sol[3],
        frailty_variance_summer = 1.0 / (fit.sol[3] * exp(δl)),
        delta_scale_summer = δs,
        summer_vs_winter_scale_ratio = exp(δs),
        delta_beta_summer = δb,
        summer_vs_winter_beta_ratio = exp(δb),
        delta_lambda_summer = δl,
        summer_vs_winter_lambda_ratio = exp(δl),
        mean_scale = mean(scales),
        median_scale = median(scales),
        mean_conditional_median_time = mean(medt),
        median_conditional_median_time = median(medt),
        mean_lambda_group = mean(lambdas),
        median_lambda_group = median(lambdas)
    )
end

function parameter_rows_reference(base_id::String, model::String, period::String, cov_label::String,
                                  selVar::Vector{Int64}, cov_names::Vector{String}, fit)
    rows = NamedTuple[]
    common = (
        base_model = base_id,
        model = model,
        effect_structure = "none",
        period = period,
        covariates = cov_label,
        nVar = length(selVar),
        K = length(fit.sol),
        negloglik = fit.f,
        AIC_total = compute_AIC(fit.f, length(fit.sol)),
        gradnorm = fit.gradnorm,
        gradient_ok = fit.gradient_ok
    )
    for (nm, ix) in [("alpha", 1), ("beta", 2), ("lambda", 3)]
        push!(rows, merge(common, (param_group="distribution", param_name=nm, param_index=ix, estimate=fit.sol[ix])))
    end
    push!(rows, merge(common, (param_group="derived", param_name="frailty_variance", param_index=0, estimate=1.0 / fit.sol[3])))
    for (k, cov_idx) in enumerate(selVar)
        push!(rows, merge(common, (param_group="covariate", param_name=cov_names[cov_idx], param_index=3+k, estimate=fit.sol[3+k])))
    end
    return rows
end

function parameter_rows_M1_effects(base_id::String, model::String, effect_structure::Symbol,
                                   cov_label::String, selVar::Vector{Int64}, cov_names::Vector{String}, fit)
    rows = NamedTuple[]
    common = (
        base_model = base_id,
        model = model,
        effect_structure = String(effect_structure),
        period = effect_structure === :none ? "pooled_session_no_period" : "pooled_with_period",
        covariates = cov_label,
        nVar = length(selVar),
        K = length(fit.sol),
        negloglik = fit.f,
        AIC_total = compute_AIC(fit.f, length(fit.sol)),
        gradnorm = fit.gradnorm,
        gradient_ok = fit.gradient_ok
    )
    for (nm, ix) in [("alpha", 1), ("beta", 2), ("lambda", 3)]
        push!(rows, merge(common, (param_group="distribution", param_name=nm, param_index=ix, estimate=fit.sol[ix])))
    end
    push!(rows, merge(common, (param_group="derived", param_name="frailty_variance_winter", param_index=0, estimate=1.0 / fit.sol[3])))

    for (k, cov_idx) in enumerate(selVar)
        push!(rows, merge(common, (param_group="covariate", param_name=cov_names[cov_idx], param_index=3+k, estimate=fit.sol[3+k])))
    end

    didx = GDF_WDB_Period_Effects_V6.delta_indices(selVar, effect_structure)
    if didx.δs != 0
        δs = fit.sol[didx.δs]
        push!(rows, merge(common, (param_group="period", param_name="delta_scale_summer", param_index=didx.δs, estimate=δs)))
        push!(rows, merge(common, (param_group="period", param_name="summer_vs_winter_scale_ratio", param_index=0, estimate=exp(δs))))
    end

    if didx.δb != 0
        δb = fit.sol[didx.δb]
        push!(rows, merge(common, (param_group="period", param_name="delta_beta_summer", param_index=didx.δb, estimate=δb)))
        push!(rows, merge(common, (param_group="period", param_name="summer_vs_winter_beta_ratio", param_index=0, estimate=exp(δb))))
    end
    if didx.δl != 0
        δl = fit.sol[didx.δl]
        push!(rows, merge(common, (param_group="period", param_name="delta_lambda_summer", param_index=didx.δl, estimate=δl)))
        push!(rows, merge(common, (param_group="period", param_name="summer_vs_winter_lambda_ratio", param_index=0, estimate=exp(δl))))
        push!(rows, merge(common, (param_group="derived", param_name="frailty_variance_summer", param_index=0, estimate=1.0 / (fit.sol[3] * exp(δl)))))
    end
    return rows
end

function add_delta_weights!(df::DataFrame, aic_col::Symbol=:AIC_total; suffix::String="")
    delta_col = Symbol("delta_" * String(aic_col) * suffix)
    weight_col = Symbol("weight_" * String(aic_col) * suffix)
    df[!, delta_col] = fill(Inf, nrow(df))
    df[!, weight_col] = zeros(Float64, nrow(df))
    ok = findall(isfinite.(df[!, aic_col]))
    if !isempty(ok)
        minAIC = minimum(df[ok, aic_col])
        df[ok, delta_col] .= df[ok, aic_col] .- minAIC
        ww = exp.(-0.5 .* df[ok, delta_col])
        df[ok, weight_col] .= ww ./ sum(ww)
    end
    return df
end

function run_one_period_effects_v6(; datapath::AbstractString,
                                     base_covariates::Vector{String},
                                     analysis_label::String,
                                     gradnorm_threshold::Float64=0.05,
                                     period_effect_bound::Float64=3.0,
                                     print_iter::Bool=false,
                                     outdir::AbstractString)
    df, Σ = read_occurrence_data_with_period(datapath)
    cov_names = filter(!=("siteID"), names(Σ))
    selVar = Int64.(covariate_indices(cov_names, base_covariates))
    cov_label = isempty(base_covariates) ? "null" : join(base_covariates, " + ")

    df_w = filter(:period_label => ==("hivernale"), df)
    df_s = filter(:period_label => ==("estivale"), df)

    bio_all = Biotope(df, Σ)
    bio_w = Biotope(df_w, Σ)
    bio_s = Biotope(df_s, Σ)

    λ_all = initial_lambda_from_df(df)
    λ_w = initial_lambda_from_df(df_w)
    λ_s = initial_lambda_from_df(df_s)

    fit_M0 = fit_reference_GDF_WDB(bio_all, selVar, λ_all; gradnorm_threshold=gradnorm_threshold, print_iter=print_iter)
    fit_w  = fit_reference_GDF_WDB(bio_w, selVar, λ_w; gradnorm_threshold=gradnorm_threshold, print_iter=print_iter)
    fit_s  = fit_reference_GDF_WDB(bio_s, selVar, λ_s; gradnorm_threshold=gradnorm_threshold, print_iter=print_iter)

    sum_M0 = summarize_reference_fit("M0_pooled_no_period", "pooled", bio_all, selVar, fit_M0)
    sum_w  = summarize_reference_fit("M2_separate", "hivernale", bio_w, selVar, fit_w)
    sum_s  = summarize_reference_fit("M2_separate", "estivale", bio_s, selVar, fit_s)

    # Deltas initiaux de M1 dérivés de M2, mais M1 garde le codage canonique été vs hiver.
    δs0 = clamp_delta(log(sum_s.mean_conditional_median_time / sum_w.mean_conditional_median_time))
    δb0 = clamp_delta(log(sum_s.beta / sum_w.beta))
    δl0 = clamp_delta(log(sum_s.lambda / sum_w.lambda))
    init_deltas = [δs0, δb0, δl0]

    # M0_session et M1a/M1b/M1c : groupes siteID × période, car les sessions sont indépendantes.
    bio_occ = GDF_WDB_Period_Effects_V6.OccurrenceBiotope(
        df, cov_names; group_col=:siteID_period, period_col=:period_estivale
    )

    # M0_session : même paramètres/covariables que M0, mais vraisemblance factorisée par session
    # siteID × période, sans effet fixe de période. Ce modèle isole l'effet de la structure
    # sessionnelle de l'effet fixe hiver/été.
    fit_M0_session = optimize_M1_effects_v6(bio_occ, selVar, :none, fit_M0.sol, init_deltas;
                                            period_effect_bound=period_effect_bound,
                                            gradnorm_threshold=gradnorm_threshold,
                                            print_iter=print_iter)
    sum_M0_session = summarize_M1_effects_fit("M0_session", :none, bio_occ, selVar, fit_M0_session)

    effect_structures = [:scale, :scale_beta, :scale_beta_lambda]
    effect_labels = Dict(:scale => "M1a_scale", :scale_beta => "M1b_scale_beta", :scale_beta_lambda => "M1c_scale_beta_lambda")

    fits_M1 = Dict{Symbol,Any}()
    summaries_M1 = NamedTuple[]
    for eff in effect_structures
        fit = optimize_M1_effects_v6(bio_occ, selVar, eff, fit_M0.sol, init_deltas;
                                     period_effect_bound=period_effect_bound,
                                     gradnorm_threshold=gradnorm_threshold,
                                     print_iter=print_iter)
        fits_M1[eff] = fit
        push!(summaries_M1, summarize_M1_effects_fit(effect_labels[eff], eff, bio_occ, selVar, fit))
    end

    f_M2_total = fit_w.f + fit_s.f
    K_M2_total = length(fit_w.sol) + length(fit_s.sol)

    comparison = DataFrame(
        base_model = String[],
        model = String[],
        period_structure = String[],
        effect_structure = String[],
        covariates = String[],
        negloglik = Float64[],
        K_total = Int[],
        AIC_total = Float64[],
        gradnorm = Float64[],
        gradient_ok = Bool[],
        niter = Int[],
        nsim = Int[],
        n_occurrences = Int[],
        n_latrines = Int[]
    )

    push!(comparison, (analysis_label, "M0", "pooled_no_period", "none", cov_label,
        fit_M0.f, length(fit_M0.sol), compute_AIC(fit_M0.f, length(fit_M0.sol)),
        fit_M0.gradnorm, fit_M0.gradient_ok, fit_M0.it, fit_M0.sim, length(bio_all.durd), bio_all.N))

    push!(comparison, (analysis_label, "M0_session", "pooled_session_no_period", "none_session", cov_label,
        fit_M0_session.f, length(fit_M0_session.sol), compute_AIC(fit_M0_session.f, length(fit_M0_session.sol)),
        fit_M0_session.gradnorm, fit_M0_session.gradient_ok, fit_M0_session.it, fit_M0_session.sim,
        length(bio_occ.durd), bio_occ.N))

    for eff in effect_structures
        fit = fits_M1[eff]
        push!(comparison, (analysis_label, effect_labels[eff], "pooled_with_period", String(eff), cov_label,
            fit.f, length(fit.sol), compute_AIC(fit.f, length(fit.sol)),
            fit.gradnorm, fit.gradient_ok, fit.it, fit.sim, length(bio_occ.durd), bio_occ.N))
    end

    push!(comparison, (analysis_label, "M2", "separate_winter_plus_summer", "fully_separate", cov_label,
        f_M2_total, K_M2_total, compute_AIC(f_M2_total, K_M2_total),
        max(fit_w.gradnorm, fit_s.gradnorm), fit_w.gradient_ok && fit_s.gradient_ok,
        fit_w.it + fit_s.it, fit_w.sim + fit_s.sim, length(bio_w.durd) + length(bio_s.durd), bio_w.N + bio_s.N))

    add_delta_weights!(comparison, :AIC_total)

    diagnostics_M2 = DataFrame([
        merge((base_model=analysis_label, covariates=cov_label), sum_w),
        merge((base_model=analysis_label, covariates=cov_label), sum_s)
    ])

    summaries_M1_df = DataFrame([merge((base_model=analysis_label, covariates=cov_label), s) for s in [sum_M0_session; summaries_M1]])

    # Quantiles conditionnels par période
    qrows = NamedTuple[]
    # M0 par période observée
    scales0 = base_scale_values_reference(bio_all, selVar, fit_M0.sol)
    beta0 = fill(fit_M0.sol[2], length(scales0))
    idx_w_all = findall(df.period_label .== "hivernale")
    idx_s_all = findall(df.period_label .== "estivale")
    append!(qrows, quantile_summary_rows(analysis_label, "M0", "hivernale", cov_label, scales0[idx_w_all], beta0[idx_w_all]; effect_structure="none"))
    append!(qrows, quantile_summary_rows(analysis_label, "M0", "estivale", cov_label, scales0[idx_s_all], beta0[idx_s_all]; effect_structure="none"))
    # M0_session
    sc0s = GDF_WDB_Period_Effects_V6.period_adjusted_scale_values(bio_occ, selVar, fit_M0_session.sol, :none)
    bt0s = GDF_WDB_Period_Effects_V6.period_adjusted_beta_values(bio_occ, fit_M0_session.sol, selVar, :none)
    append!(qrows, quantile_summary_rows(analysis_label, "M0_session", "hivernale", cov_label, sc0s[idx_w_all], bt0s[idx_w_all]; effect_structure="none_session"))
    append!(qrows, quantile_summary_rows(analysis_label, "M0_session", "estivale", cov_label, sc0s[idx_s_all], bt0s[idx_s_all]; effect_structure="none_session"))
    # M2
    scales_w = base_scale_values_reference(bio_w, selVar, fit_w.sol)
    betas_w = fill(fit_w.sol[2], length(scales_w))
    scales_s = base_scale_values_reference(bio_s, selVar, fit_s.sol)
    betas_s = fill(fit_s.sol[2], length(scales_s))
    append!(qrows, quantile_summary_rows(analysis_label, "M2", "hivernale", cov_label, scales_w, betas_w; effect_structure="fully_separate"))
    append!(qrows, quantile_summary_rows(analysis_label, "M2", "estivale", cov_label, scales_s, betas_s; effect_structure="fully_separate"))
    # M1 variants
    for eff in effect_structures
        fit = fits_M1[eff]
        sc = GDF_WDB_Period_Effects_V6.period_adjusted_scale_values(bio_occ, selVar, fit.sol, eff)
        bt = GDF_WDB_Period_Effects_V6.period_adjusted_beta_values(bio_occ, fit.sol, selVar, eff)
        append!(qrows, quantile_summary_rows(analysis_label, effect_labels[eff], "hivernale", cov_label, sc[idx_w_all], bt[idx_w_all]; effect_structure=String(eff)))
        append!(qrows, quantile_summary_rows(analysis_label, effect_labels[eff], "estivale", cov_label, sc[idx_s_all], bt[idx_s_all]; effect_structure=String(eff)))
    end
    quantiles = DataFrame(qrows)

    # Paramètres longs
    param_rows = NamedTuple[]
    append!(param_rows, parameter_rows_reference(analysis_label, "M0", "pooled", cov_label, selVar, cov_names, fit_M0))
    append!(param_rows, parameter_rows_M1_effects(analysis_label, "M0_session", :none, cov_label, selVar, cov_names, fit_M0_session))
    append!(param_rows, parameter_rows_reference(analysis_label, "M2", "hivernale", cov_label, selVar, cov_names, fit_w))
    append!(param_rows, parameter_rows_reference(analysis_label, "M2", "estivale", cov_label, selVar, cov_names, fit_s))
    for eff in effect_structures
        append!(param_rows, parameter_rows_M1_effects(analysis_label, effect_labels[eff], eff, cov_label, selVar, cov_names, fits_M1[eff]))
    end
    parameters_long = DataFrame(param_rows)

    mkpath(outdir)
    CSV.write(joinpath(outdir, "period_effects_v6_comparison.csv"), comparison)
    CSV.write(joinpath(outdir, "period_effects_v6_M2_diagnostics.csv"), diagnostics_M2)
    CSV.write(joinpath(outdir, "period_effects_v6_M1_summaries.csv"), summaries_M1_df)
    CSV.write(joinpath(outdir, "period_effects_v6_quantiles.csv"), quantiles)
    CSV.write(joinpath(outdir, "period_effects_v6_parameters_long.csv"), parameters_long)

    open(joinpath(outdir, "summary_v6.txt"), "w") do io
        println(io, "GDF-WDB period effects v6")
        println(io, "Base model: ", analysis_label, " / ", cov_label)
        println(io, "M0 and M2 use Desman.optimizeINRIA reference workflow; M0_session uses siteID × period groups without fixed period effect.")
        println(io, "M0_session and M1 variants use siteID × period groups; M1 variants use canonical summer indicator.")
        println(io, "Initial deltas from M2: delta_scale=$(δs0), delta_beta=$(δb0), delta_lambda=$(δl0)")
        println(io)
        println(io, "Model comparison:")
        show(io, MIME("text/plain"), comparison)
        println(io)
        println(io, "M2 diagnostics:")
        show(io, MIME("text/plain"), diagnostics_M2)
        println(io)
        println(io, "M1 summaries:")
        show(io, MIME("text/plain"), summaries_M1_df)
    end

    return (
        comparison = comparison,
        diagnostics_M2 = diagnostics_M2,
        summaries_M1 = summaries_M1_df,
        quantiles = quantiles,
        parameters_long = parameters_long,
        fits = merge(Dict(:M0 => fit_M0, :M0_session => fit_M0_session, :M2_winter => fit_w, :M2_summer => fit_s), Dict(Symbol("M1_" * String(k)) => v for (k,v) in fits_M1)),
        initial_deltas = (delta_scale=δs0, delta_beta=δb0, delta_lambda=δl0),
        outdir = outdir
    )
end

function run_period_effects_suite_GDF_WDB_v6(; datapath::AbstractString,
                                              base_covariate_sets::Vector{Vector{String}} = [String[], ["WidthWetBed"]],
                                              suite_outdir_name::AbstractString = "period_effect_GDF_WDB_effects_v6",
                                              gradnorm_threshold::Float64=0.05,
                                              period_effect_bound::Float64=3.0,
                                              print_iter::Bool=false)
    df_check, _ = read_occurrence_data_with_period(datapath)
    n_w = count(==("hivernale"), df_check.period_label)
    n_s = count(==("estivale"), df_check.period_label)
    @info "Occurrence rows by period" n_winter=n_w n_summer=n_s n_total=nrow(df_check)

    suite_outdir = joinpath(datapath, "outputs", suite_outdir_name)
    mkpath(suite_outdir)

    all_comparisons = DataFrame[]
    all_diagnostics = DataFrame[]
    all_summaries_M1 = DataFrame[]
    all_quantiles = DataFrame[]
    all_params = DataFrame[]
    results = Dict{String,Any}()

    for covset in base_covariate_sets
        label = isempty(covset) ? "null" : join(covset, "+")
        subdir = joinpath(suite_outdir, replace(label, "+" => "_"))
        res = run_one_period_effects_v6(
            datapath = datapath,
            base_covariates = covset,
            analysis_label = label,
            gradnorm_threshold = gradnorm_threshold,
            period_effect_bound = period_effect_bound,
            print_iter = print_iter,
            outdir = subdir
        )
        results[label] = res
        push!(all_comparisons, res.comparison)
        push!(all_diagnostics, res.diagnostics_M2)
        push!(all_summaries_M1, res.summaries_M1)
        push!(all_quantiles, res.quantiles)
        push!(all_params, res.parameters_long)
    end

    comparison_all = vcat(all_comparisons..., cols=:union)
    diagnostics_all = vcat(all_diagnostics..., cols=:union)
    summaries_M1_all = vcat(all_summaries_M1..., cols=:union)
    quantiles_all = vcat(all_quantiles..., cols=:union)
    parameters_all = vcat(all_params..., cols=:union)

    add_delta_weights!(comparison_all, :AIC_total; suffix="_all")

    CSV.write(joinpath(suite_outdir, "period_effects_v6_comparison_all.csv"), comparison_all)
    CSV.write(joinpath(suite_outdir, "period_effects_v6_M2_diagnostics_all.csv"), diagnostics_all)
    CSV.write(joinpath(suite_outdir, "period_effects_v6_M1_summaries_all.csv"), summaries_M1_all)
    CSV.write(joinpath(suite_outdir, "period_effects_v6_quantiles_all.csv"), quantiles_all)
    CSV.write(joinpath(suite_outdir, "period_effects_v6_parameters_long_all.csv"), parameters_all)
    serialize(joinpath(suite_outdir, "period_effects_v6_suite_bundle.jls"), results)

    open(joinpath(suite_outdir, "summary_v6_all.txt"), "w") do io
        println(io, "GDF-WDB period effects v6 suite")
        println(io, "Likelihood file: ", PERIOD_EFFECTS_LIKELIHOOD_FILE_V6)
        println(io, "M0/M2: reference optimizeINRIA ; M0_session: siteID × period groups without fixed period effect")
        println(io, "M1a/M1b/M1c: siteID × period groups, summer canonical effect")
        println(io)
        println(io, "Comparison all:")
        show(io, MIME("text/plain"), comparison_all)
        println(io)
        println(io, "M2 diagnostics all:")
        show(io, MIME("text/plain"), diagnostics_all)
        println(io)
        println(io, "M1 summaries all:")
        show(io, MIME("text/plain"), summaries_M1_all)
    end

    println("\nFiles written in: ", suite_outdir)
    println("  - period_effects_v6_comparison_all.csv")
    println("  - period_effects_v6_M2_diagnostics_all.csv")
    println("  - period_effects_v6_M1_summaries_all.csv")
    println("  - period_effects_v6_quantiles_all.csv")
    println("  - period_effects_v6_parameters_long_all.csv")
    println("  - period_effects_v6_suite_bundle.jls")
    println("  - summary_v6_all.txt")

    println("\nPeriod effects model comparison:")
    println(comparison_all)
    println("\nM2 diagnostics:")
    println(diagnostics_all)
    println("\nM1 summaries:")
    println(summaries_M1_all)

    return (
        comparison_all = comparison_all,
        diagnostics_all = diagnostics_all,
        summaries_M1_all = summaries_M1_all,
        quantiles_all = quantiles_all,
        parameters_all = parameters_all,
        results = results,
        outdir = suite_outdir
    )
end

# =========================================================
# Appel par défaut : modèles compétitifs GDF-WDB de référence
# =========================================================
res_period_effects_v6 = run_period_effects_suite_GDF_WDB_v6(
    datapath = datapath,
    base_covariate_sets = [String[], ["WidthWetBed"]]
)
