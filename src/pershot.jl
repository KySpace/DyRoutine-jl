function subtract_corner_mean(arr::AbstractMatrix, corner_height::Integer, corner_width::Integer)
    corner_height > 0 || throw(ArgumentError("corner_height must be positive."))
    corner_width > 0 || throw(ArgumentError("corner_width must be positive."))

    height, width = size(arr)
    2 * corner_height <= height || throw(ArgumentError("corner_height is too large for array height $height."))
    2 * corner_width <= width || throw(ArgumentError("corner_width is too large for array width $width."))

    tl = @view arr[1:corner_height, 1:corner_width]
    tr = @view arr[1:corner_height, width-corner_width+1:width]
    bl = @view arr[height-corner_height+1:height, 1:corner_width]
    br = @view arr[height-corner_height+1:height, width-corner_width+1:width]

    corner_mean = (sum(tl) + sum(tr) + sum(bl) + sum(br)) / (4 * corner_height * corner_width)
    return arr .- corner_mean
end

function calc_number(arr::AbstractMatrix)
    return sum(arr)
end
