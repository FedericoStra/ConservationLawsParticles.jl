export V, V2, W_attr, W_rep, Wprime_attr, Wprime_rep, mobility, mobρ, mobσ

V(r) = - r^3
V2(r) = - (r-1)^3
W_attr(r) = 5log(abs(r) + 1)
W_rep(r) = -5log(abs(r) + 1)
Wprime_attr(r) = 5 * sign(r) / (abs(r) + 1)
Wprime_rep(r) = - 5 * sign(r) / (abs(r) + 1)
mobility(ρ) = max(1 - ρ, 0)
mobility(ρ, σ) = max(1 - ρ - σ, 0)
mobρ(ρ, σ) = max(1 - ρ - 0.5σ, 0)
mobσ(ρ, σ) = max(1 - σ - 0.5ρ, 0)
