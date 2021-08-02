export V, V2, W_attr, W_rep, Wprime_attr, Wprime_rep, mobility, mobρ, mobσ

@time_independent V(r) = - r^3
@time_independent V2(r) = - (r-1)^3
@time_independent W_attr(r) = 5log1p(abs(r))
@time_independent W_rep(r) = -5log1p(abs(r))
@time_independent Wprime_attr(r) = 5 * sign(r) / (abs(r) + 1)
@time_independent Wprime_rep(r) = - 5 * sign(r) / (abs(r) + 1)
mobility(ρ) = max(1 - ρ, 0)
mobility(ρ, σ) = max(1 - ρ - σ, 0)
mobρ(ρ, σ) = max(1 - ρ - 0.5σ, 0)
mobσ(ρ, σ) = max(1 - σ - 0.5ρ, 0)
