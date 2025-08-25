#Symbols to be replaced before compilation

#Greek
s/α/alfa_/g ; s/`alfa_/α/g 
s/β/bta_/g  ; s/`bta_/β/g  
s/Γ/Gma_/g  ; s/`Gma_/Γ/g  
s/γ/gma_/g  ; s/`gma_/γ/g  
s/Δ/Dta_/g  ; s/`Dta_/Δ/g  
s/δ/dta_/g  ; s/`dta_/δ/g  
s/ϵ/eps_/g  ; s/`eps_/ϵ/g  
s/ε/veps_/g ; s/`veps_ε/g 
s/ζ/zta_/g  ; s/`zta_/ζ/g  
s/η/eta_/g  ; s/`eta_/η/g  
s/θ/thta_/g ; s/`thta_θ/g 
s/ι/iota_/g ; s/`iota_ι/g 
s/κ/kpa_/g  ; s/`kpa_/κ/g  
s/λ/lda_/g  ; s/`lda_/λ/g  
s/μ/mu_/g   ; s/`mu_/μ/g   
s/ν/nu_/g   ; s/`nu_/ν/g   
s/ξ/xi_/g   ; s/`xi_/ξ/g   
#s/`ο/omi_/g 
s/π/pi_/g   ; s/`pi_/π/g   
s/ρ/rho _/g ; s/`rho_/ρ/g  
s/σ/sgma_/g ; s/`sgma_/σ/g 
s/τ/tau_/g  ; s/`tau_/τ/g  
s/υ/ups_/g  ; s/`ups_/υ/g  
s/φ/phi_/g  ; s/`phi_/φ/g  
s/χ/chi_/g  ; s/`chi_/χ/g  
s/ψ/psi_/g  ; s/`psi_/ψ/g  
s/ω/omga_/g ; s/`omga_/ω/g 

#Math
#s/⁺/pls_/g
#s/⁻/mns_/g
s/˙/dot_/g  ; s/`dot_/˙/g #U+02D9 dot on top
s/·/cdot_/g ; s/`cdot_/·/g  #U+00B7
s/⨯/times_/g; s/`times_/⨯/g #U+2A2F
#s/±/pm_/g
s/∂/part_/g ; s/`part_/∂/g
s/∇/nbla_/g ; s/`nbla_/∇/g
#s/ᵀ/transp_/g
s/ᵃ/adj_/g  ; s/`adj_/ᵃ/g 
s/⁻¹/inv_/g ; s/`inv_/⁻¹/g
s/²/sq_/g   ; s/`sq_/²/g  
s/³/cbc_/g  ; s/`cbc_/³/g 
#s/⁴/forth_/g;
s/√/sqrt_/g ; s/`sqrt_/√/g
s/ₓ/x_/g    ; s/`x_/ₓ/g      #U+2093
s/ᵧ/y_/g    ; s/`y_/ᵧ/g

s/★/star_/g ; s/`star_/★/g 

##Arrows
#s/←/left_/g
#s/→/rght_/g
#s/↑/up_/g
#s/↓/dn_/g
#s/⊙/frnt_/g


#Nomenclature for public entities:
#i_   - integer type
#r_   - (single-precision) real type
#d_   - double-precision real type
#c_   - (single-precision) complex type
#z_   - double-precision complex type
#s_   - character or string type
#p_   - pointer
#t_   - derived type
#o_   - optional argument for procedures
#sub_ - subroutines
#m_   - module
#if_,is_ - logical type
#oif_,ois_ - logical type optional argument