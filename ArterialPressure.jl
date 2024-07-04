


include("functions/functions.jl")
filepath = raw"signals/Мельникова_Елизавета_Дмитриевна_21-04-22_11-43-20_.hdr"
num_ch, fs, ibeg, iend, timestart, names, lsbs, units, type = readhdr(filepath)
named_channels, fs, timestart, units = readbin(filepath)

#Для получения канала просто запрашиваем его через точку:
named_channels.LR
#все каналы
keys(named_channels)
ir = named_channels.Ir ./ 1000 # ИК сигнал
red = named_channels.Red ./ 1000# красный сигнал
ap = named_channels.FPrsNorm1  # давление не делим на 1000 чтоб получить привычные значения
ecg = named_channels.LR ./ 1000
fs = 1000

ap0 = ap[3:end]
ap1 = ap[12450:end]

#Первый фильтр через уравнение (lowpass)
range = collect(31:length(ap1))
ap_lowpassed = fill(0.0, length(ap1))
for i in range
    ap_lowpassed[i] = 2 * ap_lowpassed[i-1] - ap_lowpassed[i-2] + ap1[i] - 2 * ap1[i-15] + ap1[i-30]
end


#Второй фильтр через уравнение (highpass)
range2 = collect(775:length(ap_lowpassed))
ap_bandpassed = fill(0.0, length(ap_lowpassed))
for i in range2
    ap_bandpassed[i] = ap_bandpassed[i-1] - (1 / 774) * ap_lowpassed[i] + ap_lowpassed[i-387] - ap_lowpassed[i-388] + (1 / 774) * ap_lowpassed[i-774]
end

ap_bandpassed = ap_bandpassed[388:end] #Убираем задержку фильтрации


# Преобразуем сигнал в форму SSF с окном 64 отсчета
Awind = 64
SSF = fill(0.0, length(ap_bandpassed))
S = 0.0
k = 0
range3 = collect(1:length(ap_bandpassed))
for i in range3
    if k == Awind
        SSF[i] = S
        k = 0
        S = 0.0
    end

    if ap_bandpassed[i] > 0
        S += ap_bandpassed[i]
    end
    k += 1
end


#Находим расположение пиков и записываем их в вектор peaks_detected
SIG = SSF[1:end]
Learning = SSF[1:3000] #Обучающий участок для определения первоначального порога
SIG_LEV = 0.7 * maximum(Learning)
NOISE_LEV = 0.5 * mean(SSF)
THR_SIG = 0.7 * maximum(Learning)
#THR_SIG= NOISE_LEV + 0.25*(SIG_LEV - NOISE_LEV)
k = 0 #Зона нечувствительности
ind = 0 #счетчик индексов
Pred_p = 0
ap_peaks_detected = Vector{Int64}[]     # Координаты обнаруженных пиков
for i in collect(1:length(SIG))
    if k > 500

        if SIG[i] >= THR_SIG
            if i - 100 < 0
                O = maximum(SIG[1:i+100])
                rangeH = collect(1:i+100)
            elseif i + 100 > length(SIG)
                O = maximum(SIG[i:length(SIG)])
                rangeH = collect(i:length(SIG))
            else
                O = maximum(SIG[i-100:i+100])
                rangeH = collect(i-100:i+100)
            end
            for j in rangeH
                if SIG[j] == O
                    if j - 100 < 0
                        O1 = maximum(SIG[1:j+100])
                        rangeH1 = collect(1:j+100)
                    elseif j + 100 > length(SIG)
                        O1 = maximum(SIG[j:length(SIG)])
                        rangeH1 = collect(j:length(SIG))
                    else
                        O1 = maximum(SIG[j-100:j+100])
                        rangeH1 = collect(j-100:j+100)
                    end
                    for n in rangeH1
                        if SIG[n] == O1
                            insert!(ap_peaks_detected, ind + 1, [n])
                            SIG_LEV = 0.125 * SIG[n] + 0.575 * SIG_LEV
                            k = 0
                            ind += 1
                            Pred_p = n
                            break
                        end
                    end
                end

            end

        end

    end

    k += 1
    if Pred_p != 0
        THR_SIG = 0.7 * SSF[Pred_p]
    end
end

ap_Peaks_x = collect(Iterators.flatten(ap_peaks_detected)) # Переходим от "вектора векторов" к простому вектору значений




# ap_Peaks_x= collect(Iterators.flatten(ap_peaks_detected)) # Переходим от "вектора векторов" к простому вектору значений

#Уточняем положения найденных ранее пиков
ind17 = 0
ap_Peaks_x_updt = Vector{Int64}[]
for i in ap_Peaks_x
    if i - 200 < 0
        M = maximum(ap_bandpassed[1:i+200])
        rangePeak = collect(1:i+200)
    elseif i + 200 > length(ap_bandpassed)
        M = maximum(ap_bandpassed[i-200:length(ap_bandpassed)])
        rangePeak = collect(i-200:length(ap_bandpassed))
    else
        M = maximum(ap_bandpassed[i-200:i+200])
        rangePeak = collect(i-200:i+200)
    end
    for j in rangePeak
        if ap_bandpassed[j] == M
            insert!(ap_Peaks_x_updt, ind17 + 1, [j])
            ind17 += 1
            break

        end
    end

end

ap_Peaks_x_updt = collect(Iterators.flatten(ap_Peaks_x_updt))

ap_Peaks_y_updt = fill(0.0, length(ap_Peaks_x_updt))
index = 0
for i in ap_Peaks_x_updt
    ap_Peaks_y_updt[index+1] = ap_bandpassed[i]
    index += 1
end
# Peaks_x_updt и Peaks_y_updt координаты и значения пиков соответственно



ind18 = 0
ap_Mins_x_updt = Vector{Int64}[]
for i in ap_Peaks_x_updt

    if i + 700 > length(ap_bandpassed)
        M = minimum(ap_bandpassed[i:length(ap_bandpassed)])
        rangeMin = collect(i:length(ap_bandpassed))
    else
        M = minimum(ap_bandpassed[i:i+700])
        rangeMin = collect(i:i+700)
    end
    for j in rangeMin
        if ap_bandpassed[j] == M
            if j + 200 > length(ap_bandpassed)
                M1 = minimum(ap_bandpassed[j:length(ap_bandpassed)])
                rangeMin1 = collect(j:length(ap1))
            elseif j - 200 < 0
                M1 = minimum(ap_bandpassed[1:j])
                rangeMin1 = collect(1:j)
            else
                M1 = minimum(ap_bandpassed[j-200:j+200])
                rangeMin1 = collect(j-200:j+200)
            end
            for n in rangeMin1
                if ap_bandpassed[n] == M1
                    insert!(ap_Mins_x_updt, ind18 + 1, [n])
                    ind18 += 1
                    break
                end
            end
        end
    end
end
ap_Mins_x_updt = collect(Iterators.flatten(ap_Mins_x_updt))
ap_Mins_y_updt = fill(0.0, length(ap_Mins_x_updt))
index = 0
for i in ap_Mins_x_updt
    ap_Mins_y_updt[index+1] = ap_bandpassed[i]
    index += 1
end


# Mins_x_updt и Mins_y_updt координаты и значения минимумов соответственно


ap_Peaks_x_end = ap_Peaks_x_updt .+ 12449

#Уточняем положения найденных ранее пиков для исходного сигнала
ind19 = 0
ap_Peaks_x_updt_end = Vector{Int64}[]
for i in ap_Peaks_x_end
    if i - 200 < 0
        M = maximum(ap0[1:i+200])
        rangePeak = collect(1:i+200)
    elseif i + 200 > length(ap0)
        M = maximum(ap0[i-200:length(ap0)])
        rangePeak = collect(i-200:length(ap0))
    else
        M = maximum(ap0[i-200:i+200])
        rangePeak = collect(i-200:i+200)
    end
    for j in rangePeak
        if ap0[j] == M
            insert!(ap_Peaks_x_updt_end, ind19 + 1, [j])
            ind19 += 1
            break

        end
    end

end

ap_Peaks_x_updt_end = collect(Iterators.flatten(ap_Peaks_x_updt_end))

ap_Peaks_y_updt_end = fill(0.0, length(ap_Peaks_x_updt))
index = 0
for i in ap_Peaks_x_updt_end
    ap_Peaks_y_updt_end[index+1] = ap0[i]
    index += 1
end
# Peaks_x_updt и Peaks_y_updt координаты и значения пиков соответственно



ind20 = 0
ap_Mins_x_updt_end = Vector{Int64}[]
for i in ap_Peaks_x_updt_end

    if i + 700 > length(ap0)
        M = minimum(ap0[i:length(ap0)])
        rangeMin = collect(i:length(ap0))
    else
        M = minimum(ap0[i:i+700])
        rangeMin = collect(i:i+700)
    end
    for j in rangeMin
        if ap0[j] == M
            if j + 200 > length(ap0)
                M1 = minimum(ap0[j:length(ap0)])
                rangeMin1 = collect(j:length(ap0))
            elseif j - 200 < 0
                M1 = minimum(ap0[1:j])
                rangeMin1 = collect(1:j)
            else
                M1 = minimum(ap0[j-200:j+200])
                rangeMin1 = collect(j-200:j+200)
            end
            for n in rangeMin1
                if ap0[n] == M1
                    insert!(ap_Mins_x_updt_end, ind20 + 1, [n])
                    ind20 += 1
                    break
                end
            end
            break
        end
    end
end
ap_Mins_x_updt_end = collect(Iterators.flatten(ap_Mins_x_updt_end))
ap_Mins_y_updt_end = fill(0.0, length(ap_Mins_x_updt_end))
index = 0
for i in ap_Mins_x_updt_end
    ap_Mins_y_updt_end[index+1] = ap0[i]
    index += 1
end


# SAP,DAP,PulseAP,MeanAP=variabilityAP(ap_Peaks_x_updt_end,ap_Mins_x_updt_end,ap_Peaks_y_updt_end,ap_Mins_y_updt_end)
# SAP
# DAP
# PulseAP
# MeanAP


# press=DataFrame(Peaks_coordinates= ap_Peaks_x_updt_end,Mins_coordinates=ap_Mins_x_updt_end)
# CSV.write("Мельникова_11-43_ArterialPressure-Test",press)
# ap_Peaks_x_updt_end[181:end]




plotly()

plot([ap0], title="Артериальное давление", ylabel="Давление, мм рт.ст.", xlabel="Время, с", layout=(1, 1), legend=false)


plot(SSF)
# data = collect(30:1/fs:50)
# scatter!(t_SAD, t_SAD_y)
# scatter!(Test_PPG, Test_PPG_y)
scatter!(ap_Peaks_x_updt_end, ap_Peaks_y_updt_end)
scatter!(ap_Mins_x_updt_end, ap_Mins_y_updt_end)
scatter!(notchesXCoordinates, notchesYCoordinates)
plot([y1], title="Сигнал с обозначением макс. и мин. ПВ давления", ylabel="Давление, мм рт.ст.", xlabel="Время, с", layout=(1, 1), legend=false)
index77 = 0
Test_PPG_y = fill(0.0, length(Test_PPG))
for i in Test_PPG
    Test_PPG_y[index77+1] = y1[i]
    index77 += 1
end
index78 = 0
t_SAD_y = fill(0.0, length(t_SAD))
for i in t_SAD
    t_SAD_y[index78+1] = y1[i]
    index78 += 1
end
f = reduce(vcat, peaks_detected)

Mins = ap_Mins_x_updt_end ./ fs
Peaks = ap_Peaks_x_updt_end ./ fs

# TODO не брать пики идущие непосредственно до и после "плохих" участков сигнала

# рассчитать значения сердечного выброса на каждом ударе

#Вызов функции для нахождения дикротических впадин
notchesXCoordinates_updt_end,notchesYCoordinates_updt_end = detectDecroticNotch(ap0,ap_bandpassed, ap_Peaks_x_updt, ap_Mins_x_updt)

temp = Derivate(ap_bandpassed)

plot([ap0],layout=(1,1),legend=false)
scatter!(ap_Peaks_x_updt, ap_Peaks_y_updt)
scatter!(ap_Mins_x_updt_end, ap_Mins_y_updt_end)

scatter!(notchesXCoordinates_updt_end, notchesYCoordinates_updt_end)

temp=collect(ap_Peaks_x_updt_end[1]:ap_Mins_x_updt_end[1])
plot!([ap0[ap_Peaks_x_updt_end[1]:ap_Mins_x_updt_end[1]]],layout=(1,1),legend=false)
plot(collect(1:length(temp)),fill(mean(ap0[ap_Peaks_x_updt_end[1]:ap_Mins_x_updt_end[1]]),length(temp)),layout=(1,1),legend=false)
ap_Peaks_x_updt
ap_Mins_x_updt