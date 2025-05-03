function calculateStrokeVolume(Sa::Float64,L::Float64,signal::Vector{Float64},notchesXCoordinates::Vector{Int64},apPeaksXCoordinates::Vector{Int64},apMinsXCoordinates::Vector{Int64})
    strokeVolumes = fill(0.0, length(apMinsXCoordinates))
    strokeVolumesAPcoordinates = fill(0.0, length(apMinsXCoordinates))
    for i in 1:length(apMinsXCoordinates)
        if (apMinsXCoordinates[i]+1000 < length(signal)) && (i+1 < length(notchesXCoordinates))&&(i+1 < length(apPeaksXCoordinates))
            
            # if (notchesXCoordinates[i+1] - apMinsXCoordinates[i] > 100)&&(notchesXCoordinates[i+1] - apMinsXCoordinates[i] < 1000)
                
                P_dic_mean = sum(signal[apMinsXCoordinates[i]:notchesXCoordinates[i+1]])/length(apMinsXCoordinates[i]:notchesXCoordinates[i+1])
                temp = P_dic_mean/signal[apPeaksXCoordinates[i+1]]
                
                strokeVolumes[i] = Sa * L * temp /1000  #Результат будет в миллилитрах
                strokeVolumesAPcoordinates[i] = apPeaksXCoordinates[i]
            # end
        end

    end
    
 return strokeVolumes, strokeVolumesAPcoordinates

end
