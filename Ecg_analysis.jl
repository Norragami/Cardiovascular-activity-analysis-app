using Plots
using DSP
using Statistics
using JSON
using FFTW
using DataFrames
using CSV
include("D:\\Juliawork\\diplom\\readbin.jl")
include("D:\\Juliawork\\diplom\\slide_mean_filter.jl")
include("D:\\Juliawork\\diplom\\Variabilities.jl")
include("D:\\Juliawork\\diplom\\Sensitivity_and_PPV.jl")
filepath = raw"D:\Juliawork\diplom\signals\Евстафьева_Анастасия_Васильевна_28-04-22_10-44-28_.hdr"
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
ecg_test=ecg[3:1000]
ecg1=ecg[3:end]
#=
#Интегрирование с помощью скользящего окна
N=150 #Ширина интегрирующего окна 150 
S=0.0 # Вспомогательная переменная для подсчета значения в скобках на каждом шаге

range4=collect(N+1:length(ecg_squared))
range4_1= collect(1:N)
ecg_integrated=fill(0.0,length(ecg_squared))
for i in range4
    for j in range4_1
        S+=ecg_squared[i-(N-j)]
    end
ecg_integrated[i]= (1/N)*S
S=0.0
end

=#




#Задаем фильтры для ЭКГ
fc=25
fc2=1
ftype= Lowpass(fc/fs)
ftype2=Highpass(fc2/fs)
df = digitalfilter(ftype, Butterworth(2)) |> DF2TFilter
dh=digitalfilter(ftype2, Butterworth(2)) |> DF2TFilter

ecg_lowpassed=filt(df,ecg1)
ecg_bandpassed=filt(dh,ecg_lowpassed)


#=Вспомогательный сигнал с подавленной Т-волной
rangeT1= collect(13:length(ecg1))
ecg_vspom1=fill(0.0,length(ecg1))
for i in rangeT1
    ecg_vspom1[i]= 2*ecg_vspom1[i-1]-ecg_vspom1[i-2]+ecg1[i]-2*ecg1[i-6]+ecg1[i-12]
end
ecg_vspom1=ecg_vspom1[7:end]

rangeT2 = collect(33:length(ecg_vspom1))
ecg_vspom2=fill(0.0,length(ecg_vspom1))
for i in rangeT2
    ecg_vspom2[i]= 32*ecg_vspom1[i-16]-(ecg_vspom2[i-1]+ecg_vspom1[i]-ecg_vspom1[i-32])
end
ecg_vspom2 = ecg_vspom2[17:end]



b = [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, -32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
a = [1, -1]
ecg_vspom3=filt(PolynomialRatio(b,a),ecg_vspom1)
=#
#Производная
function Derivate(signal::Vector{Float64})

    range3=collect(3:length(signal)-2)
    ecg_derivated=fill(0.0,length(signal))
    for i in range3
    ecg_derivated[i]= (1/8)*(-signal[i-2]-2*signal[i-1]+2*signal[i+1]+signal[i+2])
    end
return ecg_derivated end

ecg_derivated=Derivate(ecg_bandpassed)

#Возведение в квадрат
ecg_squared = ecg_derivated.^2
#Скользящее среднее
ecg_integ=Slide_Mean(ecg_squared,0.150)

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
ecg1_end = ecg1[2000:end]
ecg_lowpassed_end = ecg_lowpassed[2000:end]
ecg_bandpassed_end = ecg_bandpassed[2000:end]
ecg_derivated_end = ecg_derivated[2000:end]
ecg_squared_end = ecg_squared[2000:end]
ecg_integ_end = ecg_integ[1998:end] 


#Работа с пороговыми значениями на интегрированном сигнале

range_THR=collect(1:length(ecg_integ_end))
#Обучающий интервал 2с
SPKI = 0.25* maximum(ecg_integ_end[1:2000]) 
NPKI = 0.5*mean(ecg_integ_end)
THRESHOLDI1 = NPKI + 0.25*(SPKI - NPKI)
k=200 #Зона нечувствительности
ind=0 #Счетчик индексов
skp=0 #Переменная для пропуска, если вдруг поймает Т-волну
Pred_peak=0
peaks_ecg_detected=Vector{Int64}[]
for i in range_THR
    if k>=200
        if ecg_integ_end[i] >= THRESHOLDI1
            if (k < 500) && (peaks_ecg_detected != []) #Проверка, не поймало ли Т-волну
                if i+10< length(ecg_integ_end)
                    tgi = abs((ecg_integ_end[i+10]-ecg_integ_end[i])/10)
                    tg_pred = abs((ecg_integ_end[Pred_peak+10]-ecg_integ_end[Pred_peak])/10)
                    if tgi < 0.5*tg_pred
                        skp=1
                    end
                else
                    skp=1
                end
            end

            if skp == 0
                insert!(peaks_ecg_detected,ind+1,[i])
                SPKI = 0.125*ecg_integ_end[i]+0.875*SPKI 
                k=0
                ind+=1
                Pred_peak=i
            end  
        end
         
    end  
    
    if ecg_integ_end[i] < THRESHOLDI1
        NPKI = 0.125*ecg_integ_end[i] + 0.875*NPKI
    end
    skp=0   
    k+=1   
    THRESHOLDIl = NPKI + 0.25*(SPKI - NPKI)
end

#Переходим от вектора векторов к вектору значений
Peaks_ecg_x= collect(Iterators.flatten(peaks_ecg_detected))
Peaks_ecg_y=fill(0.0,length(peaks_ecg_detected))
index=0
for i in Peaks_ecg_x
    Peaks_ecg_y[index+1]=ecg_integ_end[i]
    index+=1
end

#Функция для нахождения грубых положений Q и R
function FindQ_R_raw(ecg_integ_end::Vector{T} where T,Peaks_ecg_x::Vector{T} where T)
    ind1=0
    Q_raw=Vector{Int64}[]
    for i in Peaks_ecg_x
    if i+100 > length(ecg_integ_end)
        V=minimum(ecg_integ_end[i-200:length(ecg_integ_end)])
        range_Q = collect(i-200:length(ecg_integ_end))
    
    elseif i-200 < 0
        V=minimum(ecg_integ_end[1:i+100])
        range_Q = collect(1:i+100)
    else
        V=minimum(ecg_integ_end[i-200:i+100])
        range_Q = collect(i-200:i+100)
    end
   
    for j in range_Q  
         if ecg_integ_end[j] >= 0.3*ecg_integ_end[i]
            #ecg_integ_end[j]>(25*V)
            insert!(Q_raw,ind1+1,[j])
            ind1+=1
            break
         end
      
        end
    
    end
    Q_raw_x_end= collect(Iterators.flatten(Q_raw))
    Q_raw_y_end=fill(0.0,length(Q_raw_x_end))
    index=0
    for i in Q_raw_x_end
    Q_raw_y_end[index+1]=ecg_integ_end[i]
    index+=1
    end

    R_raw=Vector{Int64}[]
    ind2=0
    for i in Peaks_ecg_x
    if i == 1
        rangeR=collect(i+1:i+100)
    elseif i+100 > length(ecg_integ_end)
    rangeR=collect(i:length(ecg_integ_end)-1)
    else 
        rangeR=(i:i+100)
    end
    for j in rangeR
         if (round(ecg_integ_end[j],sigdigits=2)==round(ecg_integ_end[j-1],sigdigits=2))&&(round(ecg_integ_end[j],sigdigits=2)==round(ecg_integ_end[j+1],sigdigits=2))
                insert!(R_raw,ind2+1,[j])
                ind2+=1
                break
         end
        end
    end
    
    
    R_raw_x_end= collect(Iterators.flatten(R_raw))
    R_raw_y_end=fill(0.0,length(R_raw_x_end))
    index=0
    for i in R_raw_x_end
    R_raw_y_end[index+1]=ecg_integ_end[i]
    index+=1
    end

    return Q_raw_x_end,Q_raw_y_end,R_raw_x_end,R_raw_y_end
end

Q_raw_x_end,Q_raw_y_end,R_raw_x_end,R_raw_y_end = FindQ_R_raw(ecg_integ_end,Peaks_ecg_x)



#Функция для нахождения положений Q, R и S на отфильтрованном сигнале
function AccurateQ_R_S(ecg_bandpassed_end::Vector{Float64},Q_raw_x_end::Vector{Int64},R_raw_x_end::Vector{Int64}) 
    R_x_end=Vector{Int64}[]
    Q_x_end=Vector{Int64}[]
    S_x_end=Vector{Int64}[]
    ind3=0
    for i in Q_raw_x_end
        if i-100<0
            M=minimum(ecg_bandpassed_end[1:i+30])
            rangeQ_end = collect(1:i+30)
        elseif i+30 > length(ecg_bandpassed_end)
            M=minimum(ecg_bandpassed_end[i-100:length(ecg_bandpassed_end)])
            rangeQ_end = collect(i-100:length(ecg_bandpassed_end))
        else
            M=minimum(ecg_bandpassed_end[i-100:i+30])
            rangeQ_end = collect(i-100:i+30)
        end
        for j in rangeQ_end
            if ecg_bandpassed_end[j]==M
                insert!(Q_x_end,ind3+1,[j])
                ind3+=1
                break
            end
        end
    end

    ind4=0
    for i in R_raw_x_end
        if i-50 < 0
            L=maximum(ecg_bandpassed_end[1:i+50])
            rangeR_end = collect(1:i+50)
        
        elseif i+50 > length(ecg_bandpassed_end)
            L=maximum(ecg_bandpassed_end[i-50:length(ecg_bandpassed_end)])
            rangeR_end = collect(i-50:length(ecg_bandpassed_end))
        else
            L=maximum(ecg_bandpassed_end[i-50:i+50])
            rangeR_end = collect(i-50:i+50)
        end
        for j in rangeR_end
            if ecg_bandpassed_end[j]==L
                insert!(R_x_end,ind4+1,[j])
                ind4+=1
                break
            end
        end
    end
    R_x_end1= collect(Iterators.flatten(R_x_end))
    ind5=0
    for i in R_x_end1
        if i+100 > length(ecg_bandpassed_end)
            H=minimum(ecg_bandpassed_end[i:length(ecg_bandpassed_end)])
            rangeS_end = collect(i:length(ecg_bandpassed_end))
        else
            H=minimum(ecg_bandpassed_end[i:i+100])
            rangeS_end = collect(i:i+100)
        end
        for j in rangeS_end
            if ecg_bandpassed_end[j]==H
                insert!(S_x_end,ind5+1,[j])
                ind5+=1
                break
            end
        end
    end
    Q_x_end1= collect(Iterators.flatten(Q_x_end))
    R_x_end1= collect(Iterators.flatten(R_x_end))
    S_x_end1=collect(Iterators.flatten(S_x_end))
    return Q_x_end1,R_x_end1,S_x_end1
end

#Применяем функцию для нахождения положений Q и R на отфильтрованном сигнале
Q_x_end,R_x_end,S_x_end=AccurateQ_R_S(ecg_bandpassed_end,Q_raw_x_end,R_raw_x_end)

#Находим значения сигнала, для построения скаттерограммы
Q_y_end=fill(0.0,length(Q_x_end))
index=0
for i in Q_x_end
Q_y_end[index+1]=ecg_bandpassed_end[i]
index+=1
end

R_y_end=fill(0.0,length(R_x_end))
index=0
for i in R_x_end
R_y_end[index+1]=ecg_bandpassed_end[i]
index+=1
end

S_y_end=fill(0.0,length(S_x_end))
index=0
for i in S_x_end
S_y_end[index+1]=ecg_bandpassed_end[i]
index+=1
end

mRR,SDRR,MSD,rMSSD,pNN50=variabilityR_R(R_x_end)
variabilityReachTime(R_x_end,Mins_x_updt)


################################################ Конец основной части, далее построение графиков и пр.


plotly()
plot([ecg1_end[1:10000],ecg_lowpassed_end[1:10000],ecg_bandpassed_end[1:10000],ecg_derivated_end[1:10000],ecg_squared_end[1:10000],ecg_integ_end[1:10000]],layout=(6,1),legend=false)
plot([ecg1_end[1:10000],ecg_bandpassed_end[1:10000]],layout=(2,1),legend=false)
plot(data,[y1[1:15001]],title="Зарегистрированный сигнал ФПГ",xlabel="Время, с",layout=(1,1),legend=false)
plot([ecg_bandpassed_end[1:10000]],layout=(1,1),legend=false)
scatter!(R_x_end[1:12],R_y_end[1:12])
scatter!(Q_x_end[1:12],Q_y_end[1:12])
scatter!(S_x_end[1:12],S_y_end[1:12])

E=(length(ecg_bandpassed_end)-1)/fs
data=collect(0:1/fs:15)

Q_x=Q_x_end./fs
R_x=R_x_end./fs
S_x=S_x_end./fs

plot([ecg_bandpassed_end],layout=(1,1),legend=false)
scatter!(Q_x[1:6],Q_y_end[1:6])
scatter!(R_x_end,R_y_end)
scatter!(S_x[1:6],S_y_end[1:6])
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
R_x_Test=R_x_end.+2000 # Прибавляем координаты т.к убирали неинформативный участок и отнимаем задержку фильтра
R_x_Test=R_x_Test.-14
Q_x_Test=Q_x_end.+2000 # Прибавляем координаты т.к убирали неинформативный участок и отнимаем задержку фильтра
Q_x_Test=Q_x_Test.-14
S_x_Test=S_x_end.+2000 # Прибавляем координаты т.к убирали неинформативный участок и отнимаем задержку фильтра
S_x_Test=S_x_Test.-14

ab=DataFrame(ECG=ecg_bandpassed_end[1:10000])
CSV.write("ECG_signal3_short",ab)

