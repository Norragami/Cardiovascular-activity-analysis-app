
using Dates
#=filepath = raw"D:\Juliawork\Мельникова_Елизавета_Дмитриевна_21-04-22_11-43-20_.hdr"
num_ch, fs, ibeg, iend, timestart, names, lsbs, units, type = readhdr(filepath)
named_channels, fs, timestart, units = readbin(filepath,ibeg:iend)

#Для получения канала просто запрашиваем его через точку:
named_channels.LR
#все каналы
keys(named_channels)
ir = named_channels.Ir ./ 1000 # ИК сигнал
red = named_channels.Red ./ 1000# красный сигнал
ap = named_channels.FPrsNorm1 ./ 1000 # давление
ecg = named_channels.LR ./ 1000 
=#




# const string2datatype = Dict{String, DataType}(
#     "int8"    => Int8,
#     "uint8"   => UInt8,
#     "int16"   => Int16,
#     "uint16"  => UInt16,
#     "int32"   => Int32,
#     "uint32"  => UInt32,
#     "int64"   => Int64,
#     "uint64"  => UInt64,
#     "float"   => Float32,
#     "float32" => Float32,
#     "double"  => Float64,
#     "float64" => Float64
# )

"""
������ hdr-����� ���������
"""
function readhdr(filepath::AbstractString)

    lines = readlines(filepath) #, enc"windows-1251") # read and decode from windows-1251 to UTF-8 string
    lines = rstrip.(lines)

    delim = (' ', '\t')
    ln = split(lines[1], delim)
    num_ch, fs, lsb = parse(Int, ln[1]), parse(Float64, ln[2]), parse(Float64, ln[3])
    type = Int32
    if (length(ln) > 3) # optional field
        type = Int32
        #type = string2datatype[ln[4]]
    end

    ln = split(lines[2], delim)
    ibeg, iend = parse(Int, ln[1]), parse(Int, ln[2])
    timestart = parse(DateTime, ln[3])

    names = String.(split(lines[3], delim))
    strs = split(lines[4], delim)
    lsbs = parse.(Float64, strs[strs.!=""])
    units = String.(split(lines[5], delim))

    if num_ch != length(names) # ����, ���� � ������ ������� �������� ���-�� �������
        num_ch = length(names)
    end

    length(names) == length(lsbs) == length(units) || error("������ ���������� �����")

    return num_ch, fs, ibeg, iend, timestart, names, lsbs, units, type
end

"""
������ bin-����� � ��������, ����� ������ ������ hdr-����
"""
function readbin(filepath::AbstractString, range_::Union{Nothing, UnitRange{Int}} = nothing)
    # ������ �� ������
    fpath, ext = splitext(filepath)
    hdrpath = fpath * ".hdr"
    binpath = fpath * ".bin"

    num_ch, fs, _, _, timestart, names, lsbs, units, type = readhdr(hdrpath)

    offset = (range_ !== nothing) ? range_.start - 1 : 0
    
    elsize = num_ch * sizeof(type)
    byteoffset = offset * elsize # 0-based
    maxlen = (filesize(binpath) - byteoffset) / elsize #0-based
    len = (range_ !== nothing) ? min(maxlen, length(range_)) : maxlen

    if len <= 0
        data = Matrix{type}(undef, num_ch, 0)
    else
        data = Matrix{type}(undef, num_ch, Int(len))
        open(binpath, "r") do io
            #seek(io, byteoffset)
            read!(io, data)
        end
    end

    channels = [data[ch, :] .* lsbs[ch] for ch in 1:num_ch] |> Tuple # matrix -> vector of channel vectors
    sym_names = Symbol.(names) |> Tuple # column names: String -> Symbol 
    
    named_channels = NamedTuple{sym_names}(channels)
    return named_channels, fs, timestart, units
end

