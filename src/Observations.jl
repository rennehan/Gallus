module Observations

function schutte18(log_masses::AbstractArray)
    alpha = 8.80
    beta = 1.24
    scaled_stellar_masses = 10.0.^log_masses ./ 1.0e11
    
    return alpha .+ beta .* log10.(scaled_stellar_masses)
end

function gsmf(log_stellar_masses::AbstractArray, redshift::Float64, author::String)
    phi_1 = 0.0
    phi_2 = 0.0
    alpha_1 = 0.0
    alpha_2 = 0.0
    log_M_star = 0.0
    
    if author == "baldry"
        if redshift <= 0.06
            phi_1 = 3.96e-3
            phi_2 = 7.9e-4
            alpha_1 = -0.35
            alpha_2 = -1.47
            log_M_star = 10.66
        else
            println("Redshift not in range for Baldry+'12 results!")
        end
    elseif author == "tomczak"
        if redshift > 0.2 && redshift <= 0.5
            phi_1 = 10.0^(-2.54)
            phi_2 = 10.0^(-4.29)
            alpha_1 = -0.98
            alpha_2 = -1.90 
            log_M_star = 10.78
        elseif redshift > 0.5 && redshift <= 0.75
            phi_1 = 10.0^(-2.55)
            phi_2 = 10.0^(-3.15)
            alpha_1 = -0.39
            alpha_2 = -1.53 
            log_M_star = 10.70
        elseif redshift > 0.75 && redshift <= 1.0
            phi_1 = 10.0^(-2.56)
            phi_2 = 10.0^(-3.39)
            alpha_1 = -0.37
            alpha_2 = -1.61 
            log_M_star = 10.66
        elseif redshift > 1.0 && redshift <= 1.25
            phi_1 = 10.0^(-2.72)
            phi_2 = 10.0^(-3.17)
            alpha_1 = 0.30
            alpha_2 = -3.17 
            log_M_star = 10.54
        elseif redshift > 1.25 && redshift <= 1.5
            phi_1 = 1.6595e-3
            phi_2 = 3.7154e-4
            alpha_1 = -0.12
            alpha_2 = -1.56 
            log_M_star = 10.61
        elseif redshift > 1.5 && redshift <= 2.0
            phi_1 = 10.0^(-3.05)
            phi_2 = 10.0^(-3.38)
            alpha_1 = 0.04
            alpha_2 = -1.49 
            log_M_star = 10.74
        elseif redshift > 2.0 && redshift <= 2.5
            phi_1 = 10.0^(-3.80)
            phi_2 = 10.0^(-3.26)
            alpha_1 = 1.03
            alpha_2 = -1.33 
            log_M_star = 10.69
        elseif redshift > 2.5 && redshift <= 3.0
            phi_1 = 10.0^(-4.54)
            phi_2 = 10.0^(-3.69)
            alpha_1 = 1.62
            alpha_2 = -1.57 
            log_M_star = 10.74
        else
            println("Redshift not in range for Tomczak+'14 results.")
        end
    else
        println("No valid author supplied!")
    end
            
    log_dM = log_stellar_masses .- log_M_star
    term1 = phi_1 .* 10.0.^(log_dM .* alpha_1)
    term2 = phi_2 .* 10.0.^(log_dM .* alpha_2)

    log(10.0) .* exp.(-10.0.^log_dM) .* 10.0.^log_dM .* (term1 .+ term2)
end

end