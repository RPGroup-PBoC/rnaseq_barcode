using FASTX, PyCall, LightXML, CodecZlib, BioSequences, DataFrames, CSV, Distributed
@everywhere begin 
    using FASTX, PyCall, LightXML, CodecZlib, BioSequences, DataFrames, CSV

    using PyCall
    git = pyimport("git")
    glob = pyimport("glob")

    function get_seq(record)
        return identifier(record), sequence(record)
    end
    
    # Define sequence next to GFP barcode
    ref_seq = LongDNASeq("TTATTTGTACAGTT")

    function find_gfp_bc(seq; ref_seq=ref_seq, ind_filter=nothing)

        ind = findfirst(ref_seq, seq)
        if isnothing(ind)
            gfp_idx = missing
            gfp_bc = missing
        else
            gfp_idx = ind[1]-5
            if ~isnothing(ind_filter)
                if gfp_idx != ind_filter
                    gfp_idx = missing
                    gfp_bc = missing
                else
                    gfp_bc = seq[ind[1]-5:ind[1]-1]
                end
            else
                gfp_bc = seq[ind[1]-5:ind[1]-1]
            end     
        end

        return gfp_idx, gfp_bc
    end

    # Date for sequencing run of the library to map
    DATE = 20210908
    # Description to be attached to folder names
    DESCRIPTION = "_lacI_titration/"

    # Find home directory for repo
    repo = git.Repo("./", search_parent_directories=true)
    homedir = repo.working_dir
    processed_dir = homedir * "/data/processed_sequencing/$(DATE)$(DESCRIPTION)"
    datadir = homedir * "/data/raw_sequencing/$(DATE)$(DESCRIPTION)"

    outputdir = homedir * "/data/barcodes/$(DATE)$(DESCRIPTION)"
    if ~isdir(outputdir)
        mkdir(outputdir)
    end

    # Find sequencing files
    files = glob.glob(processed_dir * "*fastq.gz")

end
Threads.@threads for file in files
    println("Reading file $file...")
    r = FASTQ.Reader(GzipDecompressorStream(open(file)))
    seq_list = map(x -> get_seq(x), collect(r))#[1:n_samples])
    
    println("Extracting barcodes...")
    # Initialize dataframe to save sequences
    Names = [:id, :seq]
    df_seq = DataFrame(map(idx -> getindex.(seq_list, idx), eachindex(first(seq_list))), Names)

    # Add index and sequence length to dataframe
    df_seq[!, "seq_len"] = length.(seq_list)

    # Extract barcode
    df_seq[!, "op_bc"] = map(x->x[1:20], df_seq[!, "seq"])

    println("Extracting GFP barcodes...")
    df_seq[!, "gfp_bc"] = map(x->find_gfp_bc(x, ind_filter=27)[2], df_seq[!, "seq"])
    df_seq = dropmissing(df_seq, :gfp_bc)


    df_gfp = CSV.read("./gfp_barcode.csv", DataFrame)

    # Define if barcode is present
    bc_bool = [String(x) in df_gfp[!, "gfp_barcode_revcomp"] for x in df_seq[!, "gfp_bc"]]
    df_seq_filt = df_seq[bc_bool, :]
    
    println("Writing file...")
    output_file_name = split(split(file, "/")[end], ".")[1]
    CSV.write(outputdir * output_file_name * ".csv", df_seq_filt)
end
