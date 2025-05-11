include("utility_functions/functions.jl")
include("functions_for_main_processing/ECG_processing/exported_ECG_processing functions.jl")
include("functions_for_main_processing/PPG_processing/exported_PPG_processing functions.jl")


function formHttpResponseRrIntervals(req::HTTP.Request)

    data = JSON.parse(String(req.body))

    RrIntervals, RrIntervalsX,mRR, SDRR, MSD, rMSSD, pNN50  = getRrIntervals(data["path"])

    dataToSend = Dict("RrIntervals" => RrIntervals, "RrIntervalsX" => RrIntervalsX,"mRR"=>mRR, "SDRR"=>SDRR, "MSD"=>MSD, "rMSSD"=>rMSSD, "pNN50"=>pNN50)
    json_data = JSON.json(dataToSend)
    return json_data
end

function getRrIntervals(path::String)

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
   _,R_x_end,_=AccurateQ_R_S(ecg_bandpassed_end,Q_raw_x_end,R_raw_x_end)

  R_R = variabilityPeaks(R_x_end)
  Xcoordinate = range(1,length=length(R_R),step=1)
  
  mRR,SDRR,MSD,rMSSD,pNN50 = variabilityR_R(R_R)

 return R_R, Xcoordinate,mRR,SDRR,MSD,rMSSD,pNN50
end


function formHttpResponsePulseWaveReachTime(req::HTTP.Request)

    data = JSON.parse(String(req.body))

    PulseWaveReachTime, Intervals_X = getPulseWaveReachTime(data["path"])

    dataToSend = Dict("PulseWaveReachTime" => PulseWaveReachTime, "Intervals_X" => Intervals_X)
    json_data = JSON.json(dataToSend)
    return json_data
end

function getPulseWaveReachTime(path::String)

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
   _,R_x_end,_=AccurateQ_R_S(ecg_bandpassed_end,Q_raw_x_end,R_raw_x_end)

    PulseWaveReachTime, Intervals_X = variabilityReachTime(R_x_end,Mins_x)

    # Xcoordinate = range(1,length=length(PulseWaveReachTime),step=1)

    return PulseWaveReachTime, Intervals_X
end



function formHttpResponseHeartVolume(req::HTTP.Request)

    data = JSON.parse(String(req.body))

    HeartVolume, HeartVolumePeaksCoordinates = getHeartVolume(data["path"])

    dataToSend = Dict("HeartVolume" => HeartVolume, "HeartVolumePeaksCoordinates" => HeartVolumePeaksCoordinates)
    json_data = JSON.json(dataToSend)
    return json_data
    
end




function getHeartVolume(path)
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


    ap_Peaks_x_updt,ap_Mins_x_updt = findPeaks(ssfSignal,ap_bandpassed)

    notchesXCoordinates= detectApDecroticNotch(ap0, ap_bandpassed, ap_Peaks_x_updt, ap_Mins_x_updt)

    ap_Peaks_x_end = ap_Peaks_x_updt .+ 12449 # ВАЖНО! Переводим координаты пиков в координаты исходного сигнала
    # notchesXCoordinates = notchesXCoordinates .+ 12449 Делается внутри функции detectApDecroticNotch!!!!

    # Форматируем AP сигнал 

    ap0 = ap0[2000:end]
    ap_Peaks_x_end = ap_Peaks_x_end .- 2000
    notchesXCoordinates = notchesXCoordinates .- 2000

    apPeaks_x, apPeaks_y, apMins_x, apMins_y = apPeaksCorrection(ap_Peaks_x_end, ap0)

    # TODO Передавать значения: 907.92 площадь сечения аорты для девушки (34 мм диаметр), 12 см ширина манжеты (подобрано)
    strokeVolumes, HeartVolumePeaksCoordinates = calculateStrokeVolume(907.92, 120.0, ap0, notchesXCoordinates, apPeaks_x, apMins_x)

    return strokeVolumes, HeartVolumePeaksCoordinates

end