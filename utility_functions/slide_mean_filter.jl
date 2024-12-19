# Фильтр скользящего среднего

mutable struct SlideMeanFilter{T}
    buf::Vector{T}     # Кольцевой буфер
    k::Int      # Состояние фильтра
    need_restart::Bool # Маркер инициализации фильтра

    # конструктор объектов типа SlideMeanFilter - требуется
    # задать только размер окна
    function SlideMeanFilter{T}(window::Int) where T
        new(fill(T(0), window-1), 1, true) # буфер заполинтся нулями
    end
end

# obj - аргумент типа SlideMeanFilter
# x - один отсчет(точка) сигнала

function exe(obj::SlideMeanFilter{T}, x::Number) where T
    buf, k = obj.buf, obj.k
    # exe(obj::SlideMeanFliter{T}, x::T) - функция вызова фильтра (execute)
    if obj.need_restart # инициализация на первой точке
        # заполняем буфер первой точкой
        fill!(buf, x)
        # отмечаем, что инициализация уже не требуется
        obj.need_restart = false
    end

    sum_x = x
    # сумма всех элементов в буфере + 1
    # новая точка

    for xi in buf
        sum_x += xi
    end

    # расчет среднего

    window = length(buf) + 1
    y = sum_x/window

    # записываем в буфер новую точку

    buf[k] = x
    k += 1
    # проверка, не кончился ли буфер
    if k > length(buf)
    # возвращаемся в его начало
        k = 1
    end
    # фиксируем состояние фильтра
    obj.k = k

    return y
end

#Функция для конкретного фильтра скольpяшего среднего с окном 0.115 с
function Slide_Mean(ecg::Vector{Float64}, window::Float64, fs::Int64)
    ECG = ecg   # Убираем 2 первых отсчета, ломающих график  (Тогда на вход можем подать только начальный сигнал)
    out = fill(0.0, size(ECG)) # Заполняем массив нулями, размера ECG
    wind= window*fs # Окно в отсчетах
    range6 = collect(1:length(ECG))
    slide_flt = SlideMeanFilter{Float64}(Int64(wind)) # Создаем и применяем фильтр с заданной длиной окна
       for i in range6
           x = ECG[i]
           y = exe(slide_flt, x)
           out[i] = y
       end
       return out   # На выходе получаем сглаженный сигнал в переменной out
   end

   