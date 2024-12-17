


include("functions/functions.jl")
filepath = raw"signals\Мельникова_Елизавета_Дмитриевна2_21-04-22_13-02-11_.hdr"
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

function lowpassPPG(originSignal::Vector{Float64})
    #Первый фильтр через уравнение (lowpass)
    range = collect(31:length(originSignal))
    originalLowpassed=fill(0.0,length(originSignal))
    for i in range
    originalLowpassed[i]= 2*originalLowpassed[i-1]-originalLowpassed[i-2]+originSignal[i]-2*originSignal[i-15]+originSignal[i-30]
    end
return originalLowpassed end


originalLowpassed = lowpassPPG(originSignal)

function highpassPPG(originalLowpassed::Vector{Float64})
    #Второй фильтр через уравнение (highpass)
    range2=collect(775:length(originalLowpassed))
    ppgFullFiltered=fill(0.0,length(originalLowpassed))
    for i in range2
    ppgFullFiltered[i]= ppgFullFiltered[i-1]-(1/774)*originalLowpassed[i]+originalLowpassed[i-387]-originalLowpassed[i-388]+(1/774)*originalLowpassed[i-774]
    end
return ppgFullFiltered end

ppgFullFiltered = highpassPPG(originalLowpassed)

function formatPPGSignal(ppgFullFiltered::Vector{Float64})
    ppgFullFiltered=ppgFullFiltered[388:end] #Для выравнивания с сигналом ЭКГ потом также откинуть первые 2000 значений
    ppgFullFiltered=ppgFullFiltered[2000:end]
    ppgFullFiltered=ppgFullFiltered.+500
    
end
ppgFormatted = formatPPGSignal(ppgFullFiltered)

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

ssfSignal = convertToSSF(ppgFormatted,64)

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
 
Peaks_x, Peaks_y, Mins_x, Mins_y = detectPpgPeaks(ssfSignal,ppgFormatted)


P_x=Peaks_x./fs            # Переходим от отсчетов к секундам, строим графики
Min_x=Mins_x./fs
# E=(length(y_end_modif)-1)/fs
data=collect(0.3:1/fs:5.3)
plotly()
plot(originSignal,title="Сигнал с выделением макс. и мин. ПВ",xlabel="Время, с",layout=(1,1),legend=false)

scatter!(Peaks_x,Peaks_y)
scatter!(Mins_x,Mins_y)
# # Построим ритмограммы вариабельности Peak-Peak и времени распространения волны
# variabilityReachTime(R_x_end,Mins_x_updt)
# variabilityPeaks(Peaks_x_updt)

notchesX, notchesY = detectApDecroticNotch(originSignal,ppgFormatted,Peaks_x,Mins_x)






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