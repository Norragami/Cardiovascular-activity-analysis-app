#Задаем фильтры для ЭКГ
function designFilters(fs::Int64, fc::Int64, fc2::Int64)

    ftype = Lowpass(fc / fs)
    ftype2 = Highpass(fc2 / fs)
    df = digitalfilter(ftype, Butterworth(2)) |> DF2TFilter
    dh = digitalfilter(ftype2, Butterworth(2)) |> DF2TFilter
    return df, dh
end


#Форматирование сигналов. Убираем неинформативный период
function formatECG(ecg_original::Vector{Float64}, ecg_lowpassed::Vector{Float64}, ecg_bandpassed::Vector{Float64}, ecg_derivated::Vector{Float64}, ecg_squared::Vector{Float64}, ecg_integ::Vector{Float64})
    ecg_original_end = ecg_original[2000:end]
    ecg_lowpassed_end = ecg_lowpassed[2000:end]
    ecg_bandpassed_end = ecg_bandpassed[2000:end]
    ecg_derivated_end = ecg_derivated[2000:end]
    ecg_squared_end = ecg_squared[2000:end]
    ecg_integ_end = ecg_integ[1998:end]

    return ecg_original_end, ecg_lowpassed_end, ecg_bandpassed_end, ecg_derivated_end, ecg_squared_end, ecg_integ_end
end



#Работа с пороговыми значениями на интегрированном сигнале
function findRPeaks(ecg_integ_end::Vector{Float64})


    range_THR = collect(1:length(ecg_integ_end))
    #Обучающий интервал 2с
    SPKI = 0.25 * maximum(ecg_integ_end[1:2000])
    NPKI = 0.5 * mean(ecg_integ_end)
    THRESHOLDI1 = NPKI + 0.25 * (SPKI - NPKI)
    k = 200 #Зона нечувствительности
    ind = 0 #Счетчик индексов
    skp = 0 #Переменная для пропуска, если вдруг поймает Т-волну
    Pred_peak = 0
    peaks_ecg_detected = Vector{Int64}[]
    for i in range_THR
        if k >= 200
            if ecg_integ_end[i] >= THRESHOLDI1
                if (k < 500) && (peaks_ecg_detected != []) #Проверка, не поймало ли Т-волну
                    if i + 10 < length(ecg_integ_end)
                        tgi = abs((ecg_integ_end[i+10] - ecg_integ_end[i]) / 10)
                        tg_pred = abs((ecg_integ_end[Pred_peak+10] - ecg_integ_end[Pred_peak]) / 10)
                        if tgi < 0.5 * tg_pred
                            skp = 1
                        end
                    else
                        skp = 1
                    end
                end

                if skp == 0
                    insert!(peaks_ecg_detected, ind + 1, [i])
                    SPKI = 0.125 * ecg_integ_end[i] + 0.875 * SPKI
                    k = 0
                    ind += 1
                    Pred_peak = i
                end
            end

        end

        if ecg_integ_end[i] < THRESHOLDI1
            NPKI = 0.125 * ecg_integ_end[i] + 0.875 * NPKI
        end
        skp = 0
        k += 1
        THRESHOLDIl = NPKI + 0.25 * (SPKI - NPKI)
    end

    #Переходим от вектора векторов к вектору значений
    Peaks_ecg_x = collect(Iterators.flatten(peaks_ecg_detected))
    Peaks_ecg_y = fill(0.0, length(peaks_ecg_detected))
    index = 0
    for i in Peaks_ecg_x
        Peaks_ecg_y[index+1] = ecg_integ_end[i]
        index += 1
    end
    return Peaks_ecg_x, Peaks_ecg_y
end

#Функция для нахождения грубых положений Q и R
function FindQ_R_raw(ecg_integ_end::Vector{T} where {T}, Peaks_ecg_x::Vector{T} where {T})
    ind1 = 0
    Q_raw = Vector{Int64}[]
    for i in Peaks_ecg_x
        if i + 100 > length(ecg_integ_end)
            V = minimum(ecg_integ_end[i-200:length(ecg_integ_end)])
            range_Q = collect(i-200:length(ecg_integ_end))

        elseif i - 200 < 0
            V = minimum(ecg_integ_end[1:i+100])
            range_Q = collect(1:i+100)
        else
            V = minimum(ecg_integ_end[i-200:i+100])
            range_Q = collect(i-200:i+100)
        end

        for j in range_Q
            if ecg_integ_end[j] >= 0.3 * ecg_integ_end[i]
                #ecg_integ_end[j]>(25*V)
                insert!(Q_raw, ind1 + 1, [j])
                ind1 += 1
                break
            end

        end

    end
    Q_raw_x_end = collect(Iterators.flatten(Q_raw))
    Q_raw_y_end = fill(0.0, length(Q_raw_x_end))
    index = 0
    for i in Q_raw_x_end
        Q_raw_y_end[index+1] = ecg_integ_end[i]
        index += 1
    end

    R_raw = Vector{Int64}[]
    ind2 = 0
    for i in Peaks_ecg_x
        if i == 1
            rangeR = collect(i+1:i+100)
        elseif i + 100 > length(ecg_integ_end)
            rangeR = collect(i:length(ecg_integ_end)-1)
        else
            rangeR = (i:i+100)
        end
        for j in rangeR
            if (round(ecg_integ_end[j], sigdigits=2) == round(ecg_integ_end[j-1], sigdigits=2)) && (round(ecg_integ_end[j], sigdigits=2) == round(ecg_integ_end[j+1], sigdigits=2))
                insert!(R_raw, ind2 + 1, [j])
                ind2 += 1
                break
            end
        end
    end


    R_raw_x_end = collect(Iterators.flatten(R_raw))
    R_raw_y_end = fill(0.0, length(R_raw_x_end))
    index = 0
    for i in R_raw_x_end
        R_raw_y_end[index+1] = ecg_integ_end[i]
        index += 1
    end

    return Q_raw_x_end, Q_raw_y_end, R_raw_x_end, R_raw_y_end
end

#Функция для нахождения положений Q, R и S на отфильтрованном сигнале
function AccurateQ_R_S(ecg_bandpassed_end::Vector{Float64}, Q_raw_x_end::Vector{Int64}, R_raw_x_end::Vector{Int64})
    R_x_end = Vector{Int64}[]
    Q_x_end = Vector{Int64}[]
    S_x_end = Vector{Int64}[]
    ind3 = 0
    for i in Q_raw_x_end
        if i - 100 < 0
            M = minimum(ecg_bandpassed_end[1:i+30])
            rangeQ_end = collect(1:i+30)
        elseif i + 30 > length(ecg_bandpassed_end)
            M = minimum(ecg_bandpassed_end[i-100:length(ecg_bandpassed_end)])
            rangeQ_end = collect(i-100:length(ecg_bandpassed_end))
        else
            M = minimum(ecg_bandpassed_end[i-100:i+30])
            rangeQ_end = collect(i-100:i+30)
        end
        for j in rangeQ_end
            if ecg_bandpassed_end[j] == M
                insert!(Q_x_end, ind3 + 1, [j])
                ind3 += 1
                break
            end
        end
    end

    ind4 = 0
    for i in R_raw_x_end
        if i - 50 < 0
            L = maximum(ecg_bandpassed_end[1:i+50])
            rangeR_end = collect(1:i+50)

        elseif i + 50 > length(ecg_bandpassed_end)
            L = maximum(ecg_bandpassed_end[i-50:length(ecg_bandpassed_end)])
            rangeR_end = collect(i-50:length(ecg_bandpassed_end))
        else
            L = maximum(ecg_bandpassed_end[i-50:i+50])
            rangeR_end = collect(i-50:i+50)
        end
        for j in rangeR_end
            if ecg_bandpassed_end[j] == L
                insert!(R_x_end, ind4 + 1, [j])
                ind4 += 1
                break
            end
        end
    end
    R_x_end1 = collect(Iterators.flatten(R_x_end))
    ind5 = 0
    for i in R_x_end1
        if i + 100 > length(ecg_bandpassed_end)
            H = minimum(ecg_bandpassed_end[i:length(ecg_bandpassed_end)])
            rangeS_end = collect(i:length(ecg_bandpassed_end))
        else
            H = minimum(ecg_bandpassed_end[i:i+100])
            rangeS_end = collect(i:i+100)
        end
        for j in rangeS_end
            if ecg_bandpassed_end[j] == H
                insert!(S_x_end, ind5 + 1, [j])
                ind5 += 1
                break
            end
        end
    end
    Q_x_end1 = collect(Iterators.flatten(Q_x_end))
    R_x_end1 = collect(Iterators.flatten(R_x_end))
    S_x_end1 = collect(Iterators.flatten(S_x_end))
    return Q_x_end1, R_x_end1, S_x_end1
end


#Находим значения сигнала, для построения скаттерограммы
function findAmplitudeQRS(Q_x_end, R_x_end, S_x_end, ecg_bandpassed_end)
    
    Q_y_end = fill(0.0, length(Q_x_end))
    index = 0
    for i in Q_x_end
        Q_y_end[index+1] = ecg_bandpassed_end[i]
        index += 1
    end

    R_y_end = fill(0.0, length(R_x_end))
    index = 0
    for i in R_x_end
        R_y_end[index+1] = ecg_bandpassed_end[i]
        index += 1
    end

    S_y_end = fill(0.0, length(S_x_end))
    index = 0
    for i in S_x_end
        S_y_end[index+1] = ecg_bandpassed_end[i]
        index += 1
    end

    return Q_y_end, R_y_end, S_y_end
end