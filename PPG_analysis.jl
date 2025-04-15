
include("utility_functions/functions.jl")
include("functions_for_main_processing/PPG_processing/exported_PPG_processing functions.jl")


# TODO сделать функцию formHttpResponse, которая будет возвращать JSON

function formHttpResponsePPG(req::HTTP.Request)

    data = JSON.parse(String(req.body))
    
    outputPPG, outputPPG_X, outputPPGPeaks_x, outputPPGPeaks_y, outputPPGMins_x, outputPPGMins_y =   getPPGData(data["path"], data["startPoint"], data["endPoint"])

    dataToSend = Dict("outputPPG" => outputPPG, "outputPPG_X" => outputPPG_X, "outputPPGPeaks_x" => outputPPGPeaks_x, "outputPPGPeaks_y" => outputPPGPeaks_y, "outputPPGMins_x" => outputPPGMins_x, "outputPPGMins_y" => outputPPGMins_y)

    json_data = JSON.json(dataToSend)
    return json_data
end


# TODO  Убрать неиспользованные переменные
function getPPGData(path::String, startPoint::Int64, endPoint::Int64)
    
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
    ap = named_channels.FPrsNorm1 ./ 1000 # давление
    ecg = named_channels.LR ./ 1000
    fs=1000

    #Реализуем фильтрацию ФПГ
    originSignal= ir[3:end]

    originalLowpassed = lowpassPPG(originSignal)

    ppgFullFiltered = highpassPPG(originalLowpassed)

    ppgFormatted = formatPPGSignal(ppgFullFiltered)

    ssfSignal = convertToSSF(ppgFormatted,64)

    Peaks_x, Peaks_y, Mins_x, Mins_y = detectPpgPeaks(ssfSignal,ppgFormatted)


    #Формируем данные по необходимому участку для отправки

    outputPPG = ppgFormatted[startPoint:endPoint]
    outputPPGPeaks_x = []
    outputPPGPeaks_y = []
    outputPPGMins_x = []
    outputPPGMins_y = []
     for i in eachindex(Peaks_x)
        if Peaks_x[i]>startPoint && Peaks_x[i]<endPoint
            push!(outputPPGPeaks_x,Peaks_x[i])
            push!(outputPPGPeaks_y,Peaks_y[i])
        end
        if Mins_x[i]>startPoint && Mins_x[i]<endPoint
            push!(outputPPGMins_x,Mins_x[i])
            push!(outputPPGMins_y,Mins_y[i])
        end
    end

    Xcoordinate = range(startPoint,length=length(outputPPG),step=1)

    return outputPPG, Xcoordinate, outputPPGPeaks_x, outputPPGPeaks_y, outputPPGMins_x, outputPPGMins_y

end













# P_x=Peaks_x./fs            # Переходим от отсчетов к секундам, строим графики
# Min_x=Mins_x./fs
# E=(length(y_end_modif)-1)/fs
# data=collect(0.3:1/fs:5.3)
# plotly()
# plot(ppgFormatted,title="Сигнал с выделением макс. и мин. ПВ",xlabel="Время, с",layout=(1,1),legend=false)

# scatter!(Peaks_x,Peaks_y)
# scatter!(Mins_x,Mins_y)
# # Построим ритмограммы вариабельности Peak-Peak и времени распространения волны
# variabilityReachTime(R_x_end,Mins_x_updt)
# variabilityPeaks(Peaks_x_updt)

# notchesX, notchesY = detectApDecroticNotch(originSignal,ppgFormatted,Peaks_x,Mins_x)






#############################################################################################

#Помещаем результаты в CSV
# Peaks_x_Test=Peaks_x_updt.+2000
# Peaks_x_Test=Peaks_x_Test.-15
# Mins_x_Test=Mins_x_updt.+2000
# Mins_x_Test=Mins_x_Test.-15

# ad=DataFrame(Peaks_coordinates= Peaks_x_Test,Mins_coordinates=Mins_x_Test)
# CSV.write("Мельникова_11-43_PPG-Test",ad)

# filesdata = JSON.parsefile(raw"D:\Juliawork\Мельникова_Елизавета_Дмитриевна_21-04-22_13-02-11_.json")

# t_R = filesdata["HR"]["dataHR"]["t"]

# t_SAD = filesdata["BP"]["dataSAD"]["t"]
# R_Test=R_x_Test[1:end]
# Test_PPG=Peaks_x_Test[213:end]
# Test_ap=ap_Peaks_x_updt_end[155:end]



# t_start = filesdata["bpMeasureSettings"]["adjTime"]

# Sense_PPV(R_Test,t_R,0.007)
# Sense_PPV(Test_PPG,t_SAD,0.007)
# Sense_PPV(Test_ap,t_SAD,0.007)

# Sense_PPV(R_Test,t_R,0.010)
# Sense_PPV(Test_PPG,t_SAD,0.010)
# Sense_PPV(Test_ap,t_SAD,0.010)

# Sense_PPV(R_Test,t_R,0.015)
# Sense_PPV(Test_PPG,t_SAD,0.015)
# Sense_PPV(Test_ap,t_SAD,0.015)

# Peaks_x_Test
# t_SAD
# @info t_SAD
# @info Test_PPG