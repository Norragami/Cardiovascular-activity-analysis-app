include("utility_functions/functions.jl")
include("functions_for_main_processing/AP_processing/exported_AP_processing functions.jl")



function formHttpResponseAP(req::HTTP.Request)

    data = JSON.parse(String(req.body))

    outputAP, outputAP_Peaks_x, outputAP_Peaks_y, outputAP_Mins_x, outputAP_Mins_y =   getAPData(data["path"], data["startPoint"], data["endPoint"])

    dataToSend = Dict("outputAP" => outputAP, "outputAP_Peaks_x" => outputAP_Peaks_x, "outputAP_Peaks_y" => outputAP_Peaks_y, "outputAP_Mins_x" => outputAP_Mins_x, "outputAP_Mins_y" => outputAP_Mins_y)
    json_data = JSON.json(dataToSend)
    return json_data
end



function getAPData(path::String, startPoint::Int64, endPoint::Int64)
    
    # filepath = raw"signals\Мельникова_Елизавета_Дмитриевна2_21-04-22_13-02-11_.hdr"
    filepath = "signals/"*path*".hdr"
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
    rawAPSignal = ap[12450:end] #Убираем первичную накачку манжеты. Может быть другое значение!!!


    ap_lowpassed = lowpassAP(rawAPSignal)



    ap_bandpassed = highpassAP(ap_lowpassed)



    apFormatted = formatAP(ap_bandpassed)


    ssfSignal = converToSSF(apFormatted,64)


    ap_Peaks_x_updt = findPeaks(ssfSignal,ap_bandpassed)

    ap_Peaks_x_end = ap_Peaks_x_updt .+ 12449 # ВАЖНО! Переводим координаты пиков в координаты исходного сигнала


    apPeaks_x, apPeaks_y, apMins_x, apMins_y = apPeaksCorrection(ap_Peaks_x_end, ap0)
    

    #Формируем данные по необходимому участку для отправки

    outputAP = ap0[startPoint:endPoint]
    outputAP_Peaks_x = []
    outputAP_Peaks_y = []
    outputAP_Mins_x = []
    outputAP_Mins_y = []

     for i in eachindex(apPeaks_x)
        if apPeaks_x[i]>startPoint && apPeaks_x[i]<endPoint
            push!(outputAP_Peaks_x,apPeaks_x[i])
            push!(outputAP_Peaks_y,apPeaks_y[i])
        end
        if apMins_x[i]>startPoint && apMins_x[i]<endPoint
            push!(outputAP_Mins_x,apMins_x[i])
            push!(outputAP_Mins_y,apMins_y[i])
        end
    end

    return outputAP, outputAP_Peaks_x, outputAP_Peaks_y, outputAP_Mins_x, outputAP_Mins_y

end
# SAP,DAP,PulseAP,MeanAP=variabilityAP(ap_Peaks_x_updt_end,ap_Mins_x_updt_end,ap_Peaks_y_updt_end,ap_Mins_y_updt_end)
# SAP
# DAP
# PulseAP
# MeanAP


# press=DataFrame(Peaks_coordinates= ap_Peaks_x_updt_end,Mins_coordinates=ap_Mins_x_updt_end)
# CSV.write("Мельникова_11-43_ArterialPressure-Test",press)
# ap_Peaks_x_updt_end[181:end]




# plotly()

# plot([ap0], title="Артериальное давление", ylabel="Давление, мм рт.ст.", xlabel="Время, с", layout=(1, 1), legend=false)
# scatter!(apPeaks_x, apPeaks_y)
# scatter!(apMins_x, apMins_y)

#plot(SSF)
# data = collect(30:1/fs:50)
# scatter!(t_SAD, t_SAD_y)
# scatter!(Test_PPG, Test_PPG_y)
#scatter!(ap_Peaks_x_updt_end, ap_Peaks_y_updt_end)
#scatter!(ap_Mins_x_updt_end, ap_Mins_y_updt_end)
#scatter!(notchesXCoordinates, notchesYCoordinates)
#plot([y1], title="Сигнал с обозначением макс. и мин. ПВ давления", ylabel="Давление, мм рт.ст.", xlabel="Время, с", layout=(1, 1), legend=false)

# index77 = 0
# Test_PPG_y = fill(0.0, length(Test_PPG))
# for i in Test_PPG
#     Test_PPG_y[index77+1] = y1[i]
#     index77 += 1
# end
# index78 = 0
# t_SAD_y = fill(0.0, length(t_SAD))
# for i in t_SAD
#     t_SAD_y[index78+1] = y1[i]
#     index78 += 1
# end
# f = reduce(vcat, peaks_detected)



# TODO не брать пики идущие непосредственно до и после "плохих" участков сигнала!!!

# рассчитать значения сердечного выброса на каждом ударе

#Удалить первую впадину и последний минимум -> предусмотреть чтобы ошибочное определение впадин не шло в расчет
#Первый минимум соответствует второму максимуму и второй дикротической впадине соответствует первому минимуму

# Sa = 2 - 3.5 см2 -> 300 мм2
# L = 120 мм
#Вызов функции для нахождения дикротических впадин
#notchesXCoordinates_updt_end,notchesYCoordinates_updt_end = detectApDecroticNotch(ap0,ap_bandpassed, ap_Peaks_x_updt, ap_Mins_x_updt)

# plot([ap0],layout=(1,1),legend=false)
# scatter!(ap_Peaks_x_updt_end, ap_Peaks_y_updt_end)
# scatter!(ap_Mins_x_updt_end, ap_Mins_y_updt_end)

# scatter!(notchesXCoordinates_updt_end, notchesYCoordinates_updt_end)


# time = collect(80:1/fs:100)

# Peaks = ap_Peaks_x_updt_end ./ fs
# Mins = ap_Mins_x_updt_end ./ fs
# Notches=notchesXCoordinates_updt_end ./ fs

# plot([ap0],title="Артериальное давление", ylabel="Давление, мм рт.ст.", xlabel="Время, с",layout=(1,1),legend=false)
# scatter!(Peaks[55:78],ap_Peaks_y_updt_end[55:78])
# scatter!(Mins[55:78],ap_Mins_y_updt_end[55:78])
# scatter!(Notches[55:78],notchesYCoordinates_updt_end[55:78])




#907 мм2 площадь сечения аорты для девушки (34 мм диаметр), 120 мм ширина манжеты (подобрано)

# strokeVolumes = calculateStrokeVolume(907.0, 120.0, ap0, notchesXCoordinates_updt_end, ap_Peaks_x_updt_end, ap_Mins_x_updt_end)

# bar(strokeVolumes,fillcolor=:green,title="Ритмограмма значений ударного объема",xlabel="N сокращения",ylabel="Ударный объем, мл",legend=false)

