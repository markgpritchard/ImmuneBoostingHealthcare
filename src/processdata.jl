
function insertproportions!(data; I=:I, Y=:Y, M=:M, N=:N)
    insertcols!(data, :PatientsProportion => getproperty(data, Y) ./ getproperty(data, M))
    insertcols!(data, :StaffProportion => getproperty(data, I) ./ getproperty(data, N))
end 

datetot(date::Date) = Dates.value(date - Date("2020-03-19"))
datetot(args...) = Dates.value(Date(args...) - Date("2020-03-19"))
ttodate(t) = Date("2020-03-19") + Day(t)
