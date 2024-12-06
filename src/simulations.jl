
function vaccinatestaff(t::Number; boostdateoffset=0)
    # each date is offset by 10 days, to allow for no change in immunity for first 10 days 
    # after vaccination
    if 274 <= t <= 290  # 8 to 24 December 2020 + 10/7
        return 0.035
    elseif 295 <= t < 315  # 29 December 2020 to 18 January 2021 + 10/7
        return 0.055
    elseif 315 <= t <= 333  # 18 January to 5 February 2021 + 10/7
        return 0.025
        # values above give 89% vaccinated by 5 February 2021
    # without offset, boosting is between 16 September and 15 December 2021
    elseif 556 + boostdateoffset <= t <= 646 + boostdateoffset  # +10/7
        return 0.0125  # 68% boosted over 90 days
    else 
        return 0.0
    end
end

vaccinatestaff(date::Date; kwargs...) = vaccinatestaff(datetot(date); kwargs...)
