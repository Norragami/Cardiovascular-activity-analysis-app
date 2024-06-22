using Plots
using DSP
using Statistics
using JSON
include("D:\\Juliawork\\readbin.jl")
include("D:\\Juliawork\\slide_mean_filter.jl")
include("D:\\Juliawork\\variabilities.jl")
include("D:\\Juliawork\\Sensitivity_and_PPV.jl")
filepath = raw"D:\Juliawork\Мельникова_Елизавета_Дмитриевна2_21-04-22_13-02-11_.hdr"
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

y1= ir[3:end]
#Первый фильтр через уравнение (lowpass)
range = collect(31:length(y1))
y=fill(0.0,length(y1))
for i in range
y[i]= 2*y[i-1]-y[i-2]+y1[i]-2*y1[i-15]+y1[i-30]
end

#Второй фильтр через уравнение (highpass)
range2=collect(775:length(y))
y_end=fill(0.0,length(y))
for i in range2
y_end[i]= y_end[i-1]-(1/774)*y[i]+y[i-387]-y[i-388]+(1/774)*y[i-774]
end

y_end=y_end[388:end] #Для выравнивания с сигналом ЭКГ потом также откинуть первые 2000 значений
y_end=y_end[2000:end]
y_end=y_end.+500
# Преобразуем сигнал в форму SSF с окном 64 отсчета
Awind = 64
SSF= fill(0.0, length(y_end))
S=0.0
k=0
range3 = collect(1:length(y_end))
for i in range3
    #=for j in collect(i-Awind:i)
        if y_end[j]>0
            S+=y_end[j]
        end
    end
    SSF[i]= S
    S=0.0=#

    if k== Awind
        SSF[i]= S
        k=0
        S=0.0
    end
        
    if y_end[i] > 0
        S+=y_end[i]
    end
    k+=1
end


#Находим расположение пиков и записываем их в вектор peaks_detected
SIG = SSF[1:end]
Learning=SSF[1:3000] #Обучающий участок для определения первоначального порога
SIG_LEV=0.7*maximum(Learning)
NOISE_LEV=0.5*mean(SSF)
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
            #=insert!(peaks_detected,ind+1,[i])
            SIG_LEV = 0.125*SIG[i]+0.875*SIG_LEV 
            k=0
            ind+=1 =#
        end 
    
    end 
    if (SIG[i] < THR_SIG) && (SIG[i]!=0)
        NOISE_LEV= 0.125*SIG[i]+0.575*NOISE_LEV
    end

    k+=1
   THR_SIG = 0.325*NOISE_LEV + 0.25*(SIG_LEV - NOISE_LEV)  
end

Peaks_x= collect(Iterators.flatten(peaks_detected)) # Переходим от "вектора векторов" к простому вектору значений


#f=reduce(vcat, peaks_detected)


Peaks_x= collect(Iterators.flatten(peaks_detected)) # Переходим от "вектора векторов" к простому вектору значений

#Уточняем положения найденных ранее пиков
    ind6=0
    Peaks_x_updt=Vector{Int64}[]
    for i in Peaks_x
        if i-200<0
            M=maximum(y_end[1:i+200])
            rangePeak = collect(1:i+200)
        elseif i+200 > length(y_end)
            M=maximum(y_end[i-200:length(y_end)])
            rangePeak = collect(i-200:length(y_end))
        else
            M=maximum(y_end[i-200:i+200])
            rangePeak = collect(i-200:i+200)
        end
        for j in rangePeak
            if y_end[j]==M
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
        Peaks_y_updt[index1+1]=y_end[i]
         index1+=1
    end
# Peaks_x_updt и Peaks_y_updt координаты и значения пиков соответственно




ind7=0
Mins_x_updt=Vector{Int64}[]
    for i in Peaks_x_updt
        
        if i+700 > length(y_end)
            M=minimum(y_end[i:length(y_end)])
            rangeMin = collect(i:length(y_end))
        else
            M=minimum(y_end[i:i+700])
            rangeMin = collect(i:i+700)
        end
        for j in rangeMin
            if y_end[j]==M
                if j+200 > length(y_end)
                    M1=minimum(y_end[j:length(y_end)])
                    rangeMin1 = collect(j:length(y_end))
                elseif j-200 < 0 
                    M1=minimum(y_end[1:j])
                    rangeMin1 = collect(1:j)
                else
                    M1=minimum(y_end[j-200:j+200])
                    rangeMin1 = collect(j-200:j+200)
                end
                for n in rangeMin1
                    if y_end[n] == M1
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
        Mins_y_updt[index2+1]=y_end[i]
        index2+=1
    end


# Mins_x_updt и Mins_y_updt координаты и значения минимумов соответственно


#=
P_x=Peaks_x_updt./fs            # Переходим от отсчетов к секундам, строим графики
Min_x=Mins_x_updt./fs
E=(length(y_end_modif)-1)/fs
data=collect(0.3:1/fs:5.3)
plotly()
plot(data, y_end[300:5300],title="Сигнал с выделением макс. и мин. ПВ",xlabel="Время, с",layout=(1,1),legend=false)

scatter!(P_x[1:6],Peaks_y_updt[1:6])
scatter!(Min_x[1:6],Mins_y_updt[1:6])
# Построим ритмограммы вариабельности Peak-Peak и времени распространения волны
variabilityReachTime(R_x_end,Mins_x_updt)
variabilityPeaks(Peaks_x_updt)
=#









#Помещаем результаты в CSV
Peaks_x_Test=Peaks_x_updt.+2000
Peaks_x_Test=Peaks_x_Test.-15
Mins_x_Test=Mins_x_updt.+2000
Mins_x_Test=Mins_x_Test.-15

ad=DataFrame(Peaks_coordinates= Peaks_x_Test,Mins_coordinates=Mins_x_Test)
CSV.write("Мельникова_11-43_PPG-Test",ad)

filesdata = JSON.parsefile(raw"D:\Juliawork\Мельникова_Елизавета_Дмитриевна_21-04-22_13-02-11_.json")

t_R = filesdata["HR"]["dataHR"]["t"]

t_SAD = filesdata["BP"]["dataSAD"]["t"]
R_Test=R_x_Test[1:end]
Test_PPG=Peaks_x_Test[213:end]
Test_ap=ap_Peaks_x_updt_end[155:end]



t_start = filesdata["bpMeasureSettings"]["adjTime"]

Sense_PPV(R_Test,t_R,0.007)
Sense_PPV(Test_PPG,t_SAD,0.007)
Sense_PPV(Test_ap,t_SAD,0.007)

Sense_PPV(R_Test,t_R,0.010)
Sense_PPV(Test_PPG,t_SAD,0.010)
Sense_PPV(Test_ap,t_SAD,0.010)

Sense_PPV(R_Test,t_R,0.015)
Sense_PPV(Test_PPG,t_SAD,0.015)
Sense_PPV(Test_ap,t_SAD,0.015)

Peaks_x_Test
t_SAD
@info t_SAD
@info Test_PPG