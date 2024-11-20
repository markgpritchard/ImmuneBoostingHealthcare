
function vaccinatestaff(t::Number; boostdateoffset=0)
    if 264 <= t <= 280  # 8 to 24 December 2020
        return 0.035
    elseif 285 <= t < 305  # 29 December 2020 to 18 January 2021
        return 0.055
    elseif 305 <= t <= 323  # 18 January to 5 February 2021
        return 0.025
        # values above give 89% vaccinated by 5 February 2021
    # without offset, boosting is between 1 September and 30 November 2021
    elseif 531 + boostdateoffset <= t <= 621 + boostdateoffset  
        return 0.0125  # 68% boosted over 90 days
    else 
        return 0.0
    end
end

vaccinatestaff(date::Date; kwargs...) = vaccinatestaff(datetot(date); kwargs...)
