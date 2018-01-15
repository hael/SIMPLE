#!/usr/bin/env julia
using dotsimple_filehandling

# set constants
const max_n_segments  = 20
const n_vars_head_seg = 5
const n_bytes_header  = max_n_segments * n_vars_head_seg * 8

function xplus2( x )
    x + 2
end

# input
fname = ARGS[1]
println("Name of input file: $fname")



show(xplus2(200))

# read raw header bytes

fhandle = open(fname, "r")

byte_array   = Array{UInt8}(n_bytes_header)
n_bytes_read = readbytes!(fhandle, byte_array, n_bytes_header)
println("# bytes header: $n_bytes_header")
println("# bytes read: $n_bytes_read")
# n_bytes_read != n_bytes_header, throw exception
# println("# bytes in header: $n_bytes_header")
# println("# bytes read: $n_bytes_read")
# re-shape the header & calculate # data bytes
tmp          = reinterpret(Int64, byte_array)
cnt          = 0
n_bytes_data = 0
for i in 1:5:max_n_segments * n_vars_head_seg
    cnt += 1
    header[cnt] = file_header_segment(tmp[i],tmp[i+1],tmp[i+2],tmp[i+3],tmp[i+4])
    n_bytes_data += header[cnt].n_bytes_per_record *  header[cnt].n_records
end


# show(header)

# read raw bytes of data section
byte_array = Array{UInt8}(n_bytes_data)
n_bytes_read = readbytes!(fhandle, byte_array, n_bytes_data)
close(fhandle) # done with I/O stream
println("# bytes of data: $n_bytes_data")
println("# bytes read: $n_bytes_read")

show(transcode(String, byte_array))

# transcode into one string per line
for i in 1:max_n_segments
    if header[i].n_records > 0
        first_byte = header[i].first_data_byte - n_bytes_header
        for j in 1:header[i].n_records
            last_byte = first_byte + header[i].n_bytes_per_record - 1

            # show(reinterpret(Char, byte_array[first_byte:last_byte]))

            # show(byte_array[first_byte:last_byte])

            # show(transcode(String, byte_array[first_byte:last_byte])
            # println("first_byte: $first_byte, last_byte: $last_byte")
            first_byte = last_byte + 1
        end


        # last_byte  = first_byte + header[i].n_bytes_per_record *  header[i].n_records - 1
        # println("first_byte: $first_byte, last_byte: $last_byte")


        # show(transcode(String, byte_array[first_byte : last_byte])
    end
end

module dotsimple_filehandling

export file_header_segment, header

# define a struct for the header
struct file_header_segment
    from               :: Int64
    to                 :: Int64
    n_bytes_per_record :: Int64
    n_records          :: Int64
    first_data_byte    :: Int64
end
header = Array{file_header_segment}(max_n_segments)

end # module dotsimple_filehandling
