
#Взятие производной сигнала

function Derivate(signal::Vector{Float64})

    range3=collect(3:length(signal)-2)
    signal_derivated=fill(0.0,length(signal))
    for i in range3
    signal_derivated[i]= (1/8)*(-signal[i-2]-2*signal[i-1]+2*signal[i+1]+signal[i+2])
    end
return signal_derivated end