
function lowpassPPG(originSignal::Vector{Float64})
    #Первый фильтр через уравнение (lowpass)
    range = collect(31:length(originSignal))
    originalLowpassed=fill(0.0,length(originSignal))
    for i in range
    originalLowpassed[i]= 2*originalLowpassed[i-1]-originalLowpassed[i-2]+originSignal[i]-2*originSignal[i-15]+originSignal[i-30]
    end
return originalLowpassed end

function highpassPPG(originalLowpassed::Vector{Float64})
    #Второй фильтр через уравнение (highpass)
    range2=collect(775:length(originalLowpassed))
    ppgFullFiltered=fill(0.0,length(originalLowpassed))
    for i in range2
    ppgFullFiltered[i]= ppgFullFiltered[i-1]-(1/774)*originalLowpassed[i]+originalLowpassed[i-387]-originalLowpassed[i-388]+(1/774)*originalLowpassed[i-774]
    end
return ppgFullFiltered end


function formatPPGSignal(ppgFullFiltered::Vector{Float64})
    ppgFullFiltered=ppgFullFiltered[388:end] #Для выравнивания с сигналом ЭКГ потом также откинуть первые 2000 значений
    ppgFullFiltered=ppgFullFiltered[2000:end]
    ppgFullFiltered=ppgFullFiltered.+500 #Поднятие по амплитуде
    return ppgFullFiltered
end

function convertToSSF(ppgFormatted::Vector{Float64},Awind::Int64)
    
    # Преобразуем сигнал в форму SSF с заданным окном
    SSF= fill(0.0, length(ppgFormatted))
    S=0.0
    k=0
    range3 = collect(1:length(ppgFormatted))
    for i in range3
        if k== Awind
            SSF[i]= S
            k=0
            S=0.0
        end
            
        if ppgFormatted[i] > 0
            S+=ppgFormatted[i]
        end
        k+=1
    end
    return SSF
end

function detectPpgPeaks(ssfSignal::Vector{Float64},ppgFormatted::Vector{Float64})
    
    #Находим расположение пиков и записываем их в вектор peaks_detected
    SIG = ssfSignal[1:end]
    Learning=ssfSignal[1:3000] #Обучающий участок для определения первоначального порога
    SIG_LEV=0.7*maximum(Learning)
    NOISE_LEV=0.5*mean(ssfSignal)
    THR_SIG = 0.7*maximum(Learning)
    #THR_SIG= NOISE_LEV + 0.25*(SIG_LEV - NOISE_LEV)
    k=0 #Зона нечувствительности
    ind=0 #счетчик индексов
    peaks_detected=Vector{Int64}[]     # Координаты обнаруженных пиков
    for i in collect(1:length(SIG))
        if k> 600
        
            if SIG[i]>= THR_SIG
                if i-100<0
                    O=maximum(SIG[1:i+100])
                    rangeH = collect(1:i+100)
                elseif i+100 > length(SIG)
                    O=maximum(SIG[i:length(SIG)])
                    rangeH = collect(i:length(SIG))
                else
                    O=maximum(SIG[i-100:i+100])
                    rangeH = collect(i-100:i+100)
                end
                for j in rangeH
                    if SIG[j] == O
                        if j-100<0
                            O1=maximum(SIG[1:j+100])
                            rangeH1 = collect(1:j+100)
                        elseif j+100 > length(SIG)
                            O1=maximum(SIG[j:length(SIG)])
                            rangeH1 = collect(j:length(SIG))
                        else
                            O1=maximum(SIG[j-100:j+100])
                            rangeH1 = collect(j-100:j+100)
                        end
                        for n in rangeH1
                            if SIG[n]==O1
                                insert!(peaks_detected,ind+1,[n])
                                SIG_LEV = 0.125*SIG[n]+0.575*SIG_LEV 
                                k=0
                                ind+=1
                                break
                            end
                        end 
                    end 

                end
            
            end 
        
        end 
        if (SIG[i] < THR_SIG) && (SIG[i]!=0)
            NOISE_LEV= 0.125*SIG[i]+0.575*NOISE_LEV
        end

        k+=1
    THR_SIG = 0.325*NOISE_LEV + 0.25*(SIG_LEV - NOISE_LEV)  
    end


    Peaks_x= collect(Iterators.flatten(peaks_detected)) # Переходим от "вектора векторов" к простому вектору значений

    #Уточняем положения найденных ранее пиков
        ind6=0
        Peaks_x_updt=Vector{Int64}[]
        for i in Peaks_x
            if i-200<0
                M=maximum(ppgFormatted[1:i+200])
                rangePeak = collect(1:i+200)
            elseif i+200 > length(ppgFormatted)
                M=maximum(ppgFormatted[i-200:length(ppgFormatted)])
                rangePeak = collect(i-200:length(ppgFormatted))
            else
                M=maximum(ppgFormatted[i-200:i+200])
                rangePeak = collect(i-200:i+200)
            end
            for j in rangePeak
                if ppgFormatted[j]==M
                    insert!(Peaks_x_updt,ind6+1,[j])
                    ind6+=1
                    break
                    
                end
            end
            
        end

        Peaks_x_updt=collect(Iterators.flatten(Peaks_x_updt))
        Peaks_y_updt=fill(0.0,length(Peaks_x_updt))
        index1=0
        for i in Peaks_x_updt
            Peaks_y_updt[index1+1]=ppgFormatted[i]
            index1+=1
        end
    # Peaks_x_updt и Peaks_y_updt координаты и значения пиков соответственно




    ind7=0
    Mins_x_updt=Vector{Int64}[]
        for i in Peaks_x_updt
            
            if i+700 > length(ppgFormatted)
                M=minimum(ppgFormatted[i:length(ppgFormatted)])
                rangeMin = collect(i:length(ppgFormatted))
            else
                M=minimum(ppgFormatted[i:i+700])
                rangeMin = collect(i:i+700)
            end
            for j in rangeMin
                if ppgFormatted[j]==M
                    if j+200 > length(ppgFormatted)
                        M1=minimum(ppgFormatted[j:length(ppgFormatted)])
                        rangeMin1 = collect(j:length(ppgFormatted))
                    elseif j-200 < 0 
                        M1=minimum(ppgFormatted[1:j])
                        rangeMin1 = collect(1:j)
                    else
                        M1=minimum(ppgFormatted[j-200:j+200])
                        rangeMin1 = collect(j-200:j+200)
                    end
                    for n in rangeMin1
                        if ppgFormatted[n] == M1
                            insert!(Mins_x_updt,ind7+1,[n])
                            ind7+=1
                            break
                        end
                    end
                end
            end
        end      
        Mins_x_updt=collect(Iterators.flatten(Mins_x_updt))
        Mins_y_updt=fill(0.0,length(Mins_x_updt))
        index2=0
        for i in Mins_x_updt
            Mins_y_updt[index2+1]=ppgFormatted[i]
            index2+=1
        end


    # Mins_x_updt и Mins_y_updt координаты и значения минимумов соответственно
    return Peaks_x_updt,Peaks_y_updt,Mins_x_updt,Mins_y_updt
end
