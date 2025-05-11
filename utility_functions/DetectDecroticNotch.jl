
#Принимает на вход сырой сигнал давления, сигнал давления после полосового фильтра, 
#а также координаты минимумов и максимумов пульсовой волны давления на сигнале после полосового фильтра
function detectApDecroticNotch(signal_raw::Vector{Float64}, signal_bandpassed::Vector{Float64}, ap_peaksX::Vector{Int64}, ap_minsX::Vector{Int64})
    notchesXCoordinates = fill(0, length(ap_peaksX))
    notchesYCoordinates = fill(0.0, length(ap_peaksX))
    temp = Derivate(signal_bandpassed)

    for i in 1:length(ap_peaksX)
        peak = ap_peaksX[i]
        minx = ap_minsX[i]
        if peak < minx && minx <= length(temp)
            range_indices = peak:minx
            if !isempty(range_indices)
                Min = minimum(temp[range_indices])
                rel_idx = findfirst(x -> x == Min, temp[range_indices])
                if rel_idx !== nothing
                    firstMinimum = peak + rel_idx - 1
                    if firstMinimum < minx
                        rough_range = firstMinimum:minx
                        roughNotch = maximum(temp[rough_range])
                        for j in rough_range
                            if temp[j] == roughNotch
                                notchesXCoordinates[i] = j
                                notchesYCoordinates[i] = temp[j]
                                break
                            end
                        end
                    end
                end
            end
        end
    end

    notchesXCoordinates_updt = fill(0, length(notchesXCoordinates))
    notchesYCoordinates_updt = fill(0.0, length(notchesXCoordinates))
    for m in eachindex(notchesXCoordinates)
        x = notchesXCoordinates[m]
        if x > 100
            range = (x - 100):x
            Min = minimum(signal_bandpassed[range])
            for n in range
                if signal_bandpassed[n] == Min
                    notchesXCoordinates_updt[m] = n
                    notchesYCoordinates_updt[m] = signal_bandpassed[n]
                    break
                end
            end
        end
    end

    notchesYCoordinates_updt_end = fill(0.0, length(notchesXCoordinates))
    notchesXCoordinates_updt_end = fill(0, length(notchesXCoordinates))
    notchesXCoordinates_updt_end_raw = notchesXCoordinates_updt .+ 12449

    for m in eachindex(notchesXCoordinates_updt_end_raw)
        x = notchesXCoordinates_updt_end_raw[m]
        if x > 100 && x <= length(signal_raw)
            range = (x - 100):x
            Min = minimum(signal_raw[range])
            for n in range
                if signal_raw[n] == Min
                    notchesXCoordinates_updt_end[m] = n
                    notchesYCoordinates_updt_end[m] = signal_raw[n]
                    break
                end
            end
        end
    end

    return notchesXCoordinates_updt_end
end
