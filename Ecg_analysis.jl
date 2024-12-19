include("utility_functions/functions.jl")

include("functions_for_main_processing/ECG_processing/exported_ECG_processing functions.jl")

function formHttpResponse(req::HTTP.Request)

    data = JSON.parse(String(req.body))

    outputECG, outputQ_x, outputR_x, outputS_x, outputQ_y, outputR_y, outputS_y =   getECGData(data["path"], data["startPoint"], data["endPoint"])

    dataToSend = Dict("outputECG" => outputECG, "outputQ_x" => outputQ_x, "outputR_x" => outputR_x, "outputS_x" => outputS_x, "outputQ_y" => outputQ_y, "outputR_y" => outputR_y, "outputS_y" => outputS_y)
    json_data = JSON.json(dataToSend)
    return json_data
end

function getECGData(path::String, startPoint::Int64, endPoint::Int64)
    
    # filepath = raw"signals/Мельникова_Елизавета_Дмитриевна_21-04-22_11-43-20_.hdr"
    filepath = "signals/"*path*".hdr"
    num_ch, fs, ibeg, iend, timestart, names, lsbs, units, type = readhdr(filepath)
    named_channels, fs, timestart, units = readbin(filepath)

    #Для получения канала просто запрашиваем его через точку:
    named_channels.LR
    #все каналы
    keys(named_channels)
    ir = named_channels.Ir ./ 1000 # ИК сигнал
    red = named_channels.Red ./ 1000# красный сигнал
    ap = named_channels.FPrsNorm1 ./ 1000 # давление
    ecg = named_channels.LR ./ 1000
    fs=1000


    ecg_original=ecg[3:end]


    df,dh=designFilters(1000,25,1)
    ecg_lowpassed=filt(df,ecg_original)
    ecg_bandpassed=filt(dh,ecg_lowpassed)



    #Производная

    ecg_derivated=Derivate(ecg_bandpassed)

    #Возведение в квадрат
    ecg_squared = ecg_derivated.^2
    #Скользящее среднее
    ecg_integ=Slide_Mean(ecg_squared,0.150,1000)

    ecg_original_end, ecg_lowpassed_end, ecg_bandpassed_end, ecg_derivated_end, ecg_squared_end, ecg_integ_end = formatECG(ecg_original,ecg_lowpassed,ecg_bandpassed,ecg_derivated,ecg_squared,ecg_integ)


    Peaks_ecg_x,Peaks_ecg_y = findRPeaks(ecg_integ_end)

    Q_raw_x_end,Q_raw_y_end,R_raw_x_end,R_raw_y_end = FindQ_R_raw(ecg_integ_end,Peaks_ecg_x)

    #Применяем функцию для нахождения положений Q и R на отфильтрованном сигнале
    Q_x_end,R_x_end,S_x_end=AccurateQ_R_S(ecg_bandpassed_end,Q_raw_x_end,R_raw_x_end)


    #Находим амплитуды пиков для построения точек
    # Q_y_end,R_y_end,S_y_end=findAmplitudeQRS(Q_x_end,R_x_end,S_x_end,ecg_bandpassed_end)


    #Формируем данные по необходимому участку для отправки

    outputECG = ecg_bandpassed_end[startPoint:endPoint]
    outputQ_x = []
    outputR_x = []
    outputS_x = []
    for i in eachindex(Q_x_end)
        if Q_x_end[i]>startPoint && Q_x_end[i]<endPoint
            push!(outputQ_x,Q_x_end[i])
        end
        if R_x_end[i]>startPoint && R_x_end[i]<endPoint
            push!(outputR_x,R_x_end[i])
        end
        if S_x_end[i]>startPoint && S_x_end[i]<endPoint
            push!(outputS_x,S_x_end[i])
        end
    end

    outputQ_y, outputR_y, outputS_y = findAmplitudeQRS(outputQ_x, outputR_x, outputS_x, ecg_bandpassed_end)
    


    return outputECG, outputQ_x, outputR_x, outputS_x, outputQ_y, outputR_y, outputS_y

end

# TODO !!!!! НЕ ЗАБЫВАТЬ ЗАСИНХРОНИТЬ ВСЕ СТГНАЛЫ ПО ВРЕМЕНИ ОБРАТНО





# outputECG, outputQ_x, outputR_x, outputS_x, outputQ_y, outputR_y, outputS_y = getECGData("Мельникова_Елизавета_Дмитриевна_21-04-22_11-43-20_", 1, 2000)


# mRR,SDRR,MSD,rMSSD,pNN50=variabilityR_R(R_x_end)
# variabilityReachTime(R_x_end,Mins_x_updt)


################################################ Конец основной части, далее построение графиков и пр.

    #=Попытки построить АЧХ
    Hf = freqz(DF,0:1:50,fs)
    Hh = freqz(HH,0:1:50,fs)
    DF= digitalfilter(ftype, Butterworth(2))
    HH=digitalfilter(ftype2, Butterworth(2))
    M=20*log10.(abs.(Hf))
    L=20*log10.(abs.(Hh))

    numerator_coefs = coefb(HH)
    denominator_coefs = coefa(HH) =#

    #отбрасываем неинформативный период 2с


#  plotly()
# plot([ecg1_end[1:10000],ecg_lowpassed_end[1:10000],ecg_bandpassed_end[1:10000],ecg_derivated_end[1:10000],ecg_squared_end[1:10000],ecg_integ_end[1:10000]],layout=(6,1),legend=false)
# plot([ecg1_end[1:10000],ecg_bandpassed_end[1:10000]],layout=(2,1),legend=false)
# plot(data,[y1[1:15001]],title="Зарегистрированный сигнал ФПГ",xlabel="Время, с",layout=(1,1),legend=false)
#  plot([outputECG],layout=(1,1),legend=false)
#  scatter!(outputR_x, outputR_y)
#  scatter!(outputQ_x, outputQ_y)
#  scatter!(outputS_x, outputS_y)

# E=(length(ecg_bandpassed_end)-1)/fs
# data=collect(0:1/fs:15)

# Q_x=Q_x_end./fs
# R_x=R_x_end./fs
# S_x=S_x_end./fs

# plot([ecg_bandpassed_end],layout=(1,1),legend=false)
# scatter!(Q_x[1:6],Q_y_end[1:6])
# scatter!(R_x_end,R_y_end)
# scatter!(S_x[1:6],S_y_end[1:6])
#=
Peaks_ecg_y_filt=fill(0.0,length(Peaks_ecg_x))
index=0
for i in Peaks_ecg_x
Peaks_ecg_y_filt[index+1]=ecg_bandpassed_end[i]
index+=1
end

variabilityReachTime(R_x_end,Mins_x_updt)
@info S_x_Test
=#
# TODO Обратить внимание!!!
# R_x_Test=R_x_end.+2000 # Прибавляем координаты т.к убирали неинформативный участок и отнимаем задержку фильтра
# R_x_Test=R_x_Test.-14
# Q_x_Test=Q_x_end.+2000 # Прибавляем координаты т.к убирали неинформативный участок и отнимаем задержку фильтра
# Q_x_Test=Q_x_Test.-14
# S_x_Test=S_x_end.+2000 # Прибавляем координаты т.к убирали неинформативный участок и отнимаем задержку фильтра
# S_x_Test=S_x_Test.-14

# ab=DataFrame(ECG=ecg_bandpassed_end[1:10000])
# CSV.write("ECG_signal3_short",ab)

