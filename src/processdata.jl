
function insertproportions!(data; I=:I, Y=:Y, M=:M, N=:N)
    insertcols!(data, :PatientsProportion => getproperty(data, Y) ./ getproperty(data, M))
    insertcols!(data, :StaffProportion => getproperty(data, I) ./ getproperty(data, N))
end 
