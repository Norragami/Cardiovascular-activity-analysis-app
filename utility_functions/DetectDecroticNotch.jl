
#Принимает на вход сырой сигнал давления, сигнал давления после полосового фильтра, 
#а также координаты минимумов и максимумов пульсовой волны давления на сигнале после полосового фильтра
function detectApDecroticNotch(signal_raw::Vector{Float64}, signal_bandpassed::Vector{Float64}, ap_peaksX::Vector{Int64}, ap_minsX::Vector{Int64})
    notchesXCoordinates = fill(0, length(ap_peaksX))
    notchesYCoordinates = fill(0.0, length(ap_peaksX))
    temp = Derivate(signal_bandpassed)

    for i in 1:length(ap_peaksX)

        # tempX=collect(ap_peaksX[i]:ap_minsX[i])
        Min = minimum(temp[ap_peaksX[i]:ap_minsX[i]])
        firstMinimum = ap_peaksX[i] + findfirst(x -> x == Min, temp[ap_peaksX[i]:ap_minsX[i]])
        k = 0
        range = collect(firstMinimum:ap_minsX[i])
        roughNotch = maximum(temp[firstMinimum:ap_minsX[i]])
        for j in range
            if temp[j] == roughNotch
                notchesXCoordinates[i] = j
                notchesYCoordinates[i] = temp[j]
                break
            end
        end


    end
    notchesXCoordinates_updt = fill(0, length(notchesXCoordinates))
    notchesYCoordinates_updt = fill(0.0, length(notchesXCoordinates))
    for m in eachindex(notchesXCoordinates)
        Min = minimum(signal_bandpassed[notchesXCoordinates[m]-100:notchesXCoordinates[m]])
        range = collect(notchesXCoordinates[m]-100:notchesXCoordinates[m])
        for n in range
            if signal_bandpassed[n] == Min
                notchesXCoordinates_updt[m] = n
                notchesYCoordinates_updt[m] = signal_bandpassed[n]
            end
        end

    end

    notchesYCoordinates_updt_end = fill(0.0, length(notchesXCoordinates))
    notchesXCoordinates_updt_end = fill(0, length(notchesXCoordinates))

    notchesXCoordinates_updt_end_raw = notchesXCoordinates_updt .+ 12449  #TODO Прибавляем участок который исключали из сигнала ap0
    for m in eachindex(notchesXCoordinates_updt_end_raw)
        Min = minimum(signal_raw[notchesXCoordinates_updt_end_raw[m]-100:notchesXCoordinates_updt_end_raw[m]])
        range = collect(notchesXCoordinates_updt_end_raw[m]-100:notchesXCoordinates_updt_end_raw[m])
        for n in range
            if signal_raw[n] == Min
                notchesXCoordinates_updt_end[m] = n
                notchesYCoordinates_updt_end[m] = signal_raw[n]
            end
        end

    end

    return notchesXCoordinates_updt_end # notchesYCoordinates_updt_end
end