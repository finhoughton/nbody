using PartialFunctions

# these functions are written generically to save and read a vector of any type

const delimiter::String = ";" # can't use comma because vectors use commas

split_delim(l) = split(l, delimiter)

function save!(file::IOStream, items::Vector{T}, data::Vector{Int64})::Nothing where {T}
    push!(data, length(items))
    # append the number of items to the data

    join(file, data, delimiter)
    # write the data to the file, joined by the delimiter

    write(file, "\n")
    # write a newline between data and the items

    fields::Tuple{Vararg{Symbol}} = fieldnames(Particle)
    for item ∈ items
        join(file, map(getfield $ item, fields), delimiter)
        # write each attribute of the item to the file, joined by delimiters
        # have to use map instead of broadcast
        # because `length` may not be defined for `item`

        write(file, "\n")
        # write a newline between each item
    end
    nothing
end

function read!(stream::IOStream, T::DataType)::Tuple{Tuple{Vararg{Int64}}, Vector}
    (data..., num_items) = stream |> readline |> split_delim |> (x -> parse.(Int64, x)) |> Tuple
    # reading the data at the top of the file, the last number is always the number of items.
    # steps to read the numbers:
    # - read the first line of the file, where the numbers are, to a string
    # - split the string on delimiters
    # - parse each substring  to an Int64
    # - convert to a tuple
    # - unpack those into data and the number of items

    items::Vector{T} = Vector{T}(undef, num_items)
    # create an empty vector `num_items` long

    for (idx, line) ∈ enumerate(eachline(stream))
        item_data::Vector = [(data |> Meta.parse |> eval |> type) for (type, data) ∈ (line |> split_delim |> zip $ T.types)]
        # `line |> split_delim |> zip $ T.types` splits the line on delimiters
        # then matches up the attributes of `T` to the attributes that were written to the file
        # example output: `[(Int64, "22"), (String, "John"), (Float64, "89.5")]`

        # `Meta.parse` then parses the string to an expression,
        # which is then evaluated by `eval`
        # and then conveted to the type from T.types

        items[idx] = T(item_data...)
        # item_data then conveted to the type passed in by the caller and put in the vector.
    end
    return (data, items)
end

function get_untitled_filename_number()::Int
    files::Vector{String} = filter(n -> endswith(n, ".txt") && startswith(n, "untitled"), readdir("data"))
    # get a list of files in the data/ directory, and filter by txt files starting with "untitled"
    if isempty(files)
        return 1
    end
    filenums::Vector{Int} = map((parse $ Int) ∘ String ∘ collect ∘ (Iterators.takewhile $ isdigit) ∘ (Iterators.dropwhile $ !isdigit), files)
    # extract the numbers from the untitled files
    return max(filenums...) + 1
end
