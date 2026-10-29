# 最后一次修订时间: 2025-12-10 19:12

using Pkg

# 定义需要的包列表
required_packages = [
    "CSV", "DataFrames", "Plots", 
    "StatsBase", "Distributions", "DelimitedFiles", 
    "MultipleTesting" 
    ]

# 检查并加载每个包
for pkg in required_packages
    try
        eval(Meta.parse("using $pkg"))
        println("$pkg 已加载")
    catch
        println("$pkg 未安装，正在安装...")
        Pkg.add(pkg)
        eval(Meta.parse("using $pkg"))
        println("$pkg 安装并加载成功")
    end
end


"""
    process_genomic_data(temp_path::String, file_names::Vector{String})::DataFrame

Encapsulate the genomic data processing workflow into a single function.

This function reads multiple genomic data files, filters chromosomes, and combines 
the data into a single DataFrame with positional information and coverage data.

# Arguments
- `temp_path::String`: The path to the directory containing the data files
- `file_names::Vector{String}`: A vector of file names to process

# Returns
- `DataFrame`: A combined DataFrame containing:
    - Positional information (Chromosome, BinLowEdge, BinUpEdge)
    - Coverage data from multiple samples

# Examples
```julia
temp_path = "C:/Users/yghe_/Documents/Work/Julia/DuanAnQi/293细胞每1KbRD值/代码重构/"
file_names = ["E9-E10_his_rd_p_1000_allchr_new.csv", ...]
combineData = process_genomic_data(temp_path, file_names)
```

# Exceptions
- Throws ArgumentError if temp_path is empty
- Throws ArgumentError if file_names is empty
"""
function process_genomic_data(temp_path::String, file_names::Vector{String})::DataFrame
    # 输入验证
    if isempty(temp_path)
        throw(ArgumentError("temp_path must be a non-empty string"))
    end
    
    if isempty(file_names)
        throw(ArgumentError("file_names must be a non-empty vector"))
    end

    # 读取所有数据文件
    dfs = Vector{DataFrame}(undef, length(file_names))
    for i in 1:length(file_names)
        full_path = joinpath(temp_path, file_names[i])
        dfs[i] = CSV.read(full_path, DataFrame)
    end

    # 定义染色体字符串用于过滤
    chr_strings = [
        "chr1", "chr2", "chr3", "chr4", "chr5", 
        "chr6", "chr7", "chr8", "chr9", "chr10", 
        "chr11", "chr12", "chr13", "chr14", "chr15", 
        "chr16", "chr17", "chr18", "chr19", "chr20", 
        "chr21", "chr22" 
    ]

    # 为具有确定染色体的行创建标志
    n_rows = size(dfs[1], 1)
    chr_flag = zeros(n_rows, 1)

    for i in 1:n_rows
        chr_flag[i] = (sum(dfs[1].Chromosome[i] .== chr_strings) > 0)
    end

    # 查找位置完全确定的索引
    defined_indices = findall(chr_flag[:] .== 1)

    # 创建包含位置信息的组合DataFrame
    combineData = DataFrame(
        Chromosome = dfs[1][defined_indices, :Chromosome],
        BinLowEdge = dfs[1][defined_indices, :BinLowEdge],
        BinUpEdge = dfs[1][defined_indices, :BinUpEdge],
    )

    # 添加其他文件的数据
    for i in 1:length(dfs)
        # 从文件名生成列名
        column_prefix = split(file_names[i], "_")[1]
        column_name = replace(column_prefix, "-" => "")
        column_symbol = Symbol("df" * column_name)
        
        combineData[!, column_symbol] = dfs[i][defined_indices, :Content]
    end

    return combineData
end

"""
    moving_average_histogram(depth::Vector, Xbar::Vector, windowSize::Int)

计算给定向量的移动平均值并绘制直方图，在指定位置添加垂直线标记

# 参数说明:
- `depth`: 输入数据向量
- `Xbar`: 需要在直方图上标记的位置向量
- `windowSize`: 移动平均窗口大小

# 返回结果:
- `p`: 绘制的图表对象
- `moving_avg`: 计算得到的移动平均值向量

# 异常处理:
- 当窗口大小超过向量长度或小于1时抛出错误
"""
function moving_average_histogram(depth::Vector, Xbar::Vector, windowSize::Int)
    # 参数验证
    if windowSize > length(depth)
        throw(ArgumentError("窗口大小不能超过向量长度"))
    end
    
    if windowSize < 1
        throw(ArgumentError("窗口大小必须为正整数"))
    end
    
    # 计算移动平均
    n = length(depth)
    moving_avg = Vector{Float64}(undef, n - windowSize + 1)
    
    for i in 1:(n - windowSize + 1)
        window = depth[i:(i + windowSize - 1)]
        moving_avg[i] = mean(window)
    end
    
    # 创建直方图
    p = histogram(moving_avg, 
                  alpha=0.7,
                  fillcolor=:lightblue)
    
    # 获取当前绘图区域的y轴范围
    ylims_current = ylims()
    max_height = ylims_current[2]
    
    # 在Xbar位置画红色竖线，长度为当前窗口高度
    for x_val in Xbar
        # 绘制从y轴底部到顶部的竖线
        plot!([x_val, x_val], [0, max_height], 
              color=:red, 
              linewidth=2,
              legend = false)
    end
    
    # 显示图形
    display(p)
    
    # 返回结果
    #return p, moving_avg
end

"""
    paramFineTuning(param, reference_esti, target_esti, refComponent_idx, targetCompnent_idx)

手工微调参数，修正样本带来的参数估计不准确

调整均值以及标准差，不改变组分权重
调整方式：在参考值基础上进行平移和缩放

# 参数说明:
- `param`: 分布参数，来自函数fitting_gmm()
- `reference_esti`: 整数，指出标准数据在param中的位置
- `target_esti`: 整数，指出待修正数据在param中的位置
- `refComponent_idx`: 向量，参考的组分位置，不得少于2个
- `targetCompnent_idx`: 向量，待修正组分的位置

# 返回结果:
- `param`: 微调后的参数

# 异常处理:
- 如果输入参数类型不正确，抛出ArgumentError
- 如果refComponent_idx少于2个元素，抛出ArgumentError
"""
function paramFineTuning(param, reference_esti, target_esti, refComponent_idx, targetCompnent_idx)
    try
        # 输入验证
        if !isa(param, Vector)
            throw(ArgumentError("param必须是Vector类型"))
        end
        
        if !isa(reference_esti, Integer)
            throw(ArgumentError("reference_esti必须是整数类型"))
        end
        
        if !isa(target_esti, Integer)
            throw(ArgumentError("target_esti必须是整数类型"))
        end
        
        if !isa(refComponent_idx, Vector) || length(refComponent_idx) < 2
            throw(ArgumentError("refComponent_idx必须是至少包含2个元素的向量"))
        end
        
                
        # 线性最小二乘法，计算缩放参数和平移参数
        # 使用反斜杠运算符, 计算y = a * x + b, a为缩放参数，b为平移参数
        x = param[reference_esti].means[refComponent_idx]
        y = param[target_esti].means[refComponent_idx]
        
        num_means = length(refComponent_idx)
        X = hcat(x, ones(num_means))
        coefficients = X \ y  # a = coefficients[1]，b = coefficients[2]
        
        # 计算means和stds新值,取代旧值
        param[target_esti].means[targetCompnent_idx] = 
            coefficients[1] * param[reference_esti].means[targetCompnent_idx] .+ coefficients[2]
            
        param[target_esti].stds[targetCompnent_idx] = 
            coefficients[1] * param[reference_esti].stds[targetCompnent_idx] .+ coefficients[2]
        
        return param
        
    catch e
        @error "在paramFineTuning函数中发生错误: " * string(e)
        rethrow(e)
    end
end

"""
    plotCNVregionMulti(metaDataDF, CNVscaffold, pairedCNV, sample_pair; layout = (1,1), positive_cnv_content = 0.8, chr = "all")

绘制多个样本对的CNV区域图

# 参数说明:
- `metaDataDF`: 经过过滤且标准化后的数据框架，来自函数filtering_normalizing_seqDepth()
- `CNVscaffold`: CNV片段框架矩阵
- `pairedCNV`: 配对样本的CNV鉴别情况，3D数组，来自函数copyNumber_comparison()
- `sample_pair`: 一个具有两个数值元素的向量，表示样本对，如[1,2]表示第一个样本对和第二个样本对
- `layout`: 布局参数，默认值为(1, 1)
- `plotsize`: 图形大小，默认值为(1000, 500)
- `positive_cnv_content`: CNV区域内被判定为CNV的基因组片段覆盖率最低值，默认值为0.8
- `chr`: 要绘制的染色体，默认为"all"表示所有染色体

# 返回结果:
- 数据框架，所有染色体的全部可绘制CNV区域
- 图形对象，可直接显示或进一步操作

# 依赖库:
- Plots: 绘图库
"""
function plotCNVregionMulti_RD(
    metaDataDF,
    CNVscaffold,
    pairedCNV,
    sample_pair;
    layout = (1, 1),
    xlims_pre = [0 0],
    ylims_pre = [0 0],
    plotsize = (1000, 500),
    positive_cnv_content = 0.8,
    chr = "all"
)
    try
        # 输入验证
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF必须是DataFrame类型"))
        end
     
        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold必须是n行3列的矩阵"))
        end
     
        if !isa(pairedCNV, Array) || length(size(pairedCNV)) != 3
            throw(ArgumentError("pairedCNV必须是3D的数组"))
        end
        if !isa(sample_pair, Matrix) || size(sample_pair, 2) != 2
            throw(ArgumentError("sample_pair必须是n行2列的矩阵"))
        end
     
        num_pairs = size(sample_pair, 1)
        # 创建绘图
        plt = plot(layout = layout, size = plotsize)
        xlims_current = zeros(num_pairs, 2)
        ylims_current = zeros(num_pairs, 2)
        plotcnv_regionDF = DataFrame()
        for pair_i = 1:num_pairs
            # 数据准备
            cnv_regionsDF_samplepair = cnv_region(
                metaDataDF,
                CNVscaffold,
                pairedCNV[:, sample_pair[pair_i, 1], sample_pair[pair_i, 2]],
                positive_cnv_content
            )
         
            metaDataDF_samplepair = select(
                metaDataDF,
                1, 2, 3,
                3 + sample_pair[pair_i, 1],
                3 + sample_pair[pair_i, 2]
            )
         
            depths_metaData = Matrix(metaDataDF_samplepair[:, 4:end])
            # 计算CNV框架片段深度均值
            num_fragments = size(CNVscaffold, 1)
            depths_mean = zeros(num_fragments, 2)
            fragments_chr = String[]
            fragments_chr_startposi = Int64[]
            fragments_chr_endposi = Int64[]
         
            for i = 1:num_fragments
                push!(fragments_chr, String(metaDataDF_samplepair[CNVscaffold[i, 1], 1]))
                push!(fragments_chr_startposi, Int64(metaDataDF_samplepair[CNVscaffold[i, 1], 2]))
                push!(fragments_chr_endposi, Int64(metaDataDF_samplepair[CNVscaffold[i, 2], 3]))
             
                # 添加边界检查
                start_idx = max(1, CNVscaffold[i, 1])
                end_idx = min(size(depths_metaData, 1), CNVscaffold[i, 2])
                if start_idx <= end_idx
                    depths_mean[i, 1] = mean(depths_metaData[start_idx:end_idx, 1])
                    depths_mean[i, 2] = mean(depths_metaData[start_idx:end_idx, 2])
                else
                    depths_mean[i, 1] = 0
                    depths_mean[i, 2] = 0
                end
            end
          
            # 为每个CNV区域计算mean_depth_A, mean_depth_B, depth_diff
            if !isempty(cnv_regionsDF_samplepair)
                mean_depth_A = Float64[]
                mean_depth_B = Float64[]
                depth_diff = Float64[]
                for j in 1:nrow(cnv_regionsDF_samplepair)
                    start_bin = cnv_regionsDF_samplepair[j, :start_pos_metaData]
                    end_bin = cnv_regionsDF_samplepair[j, :end_pos_metaData]
                    if start_bin <= end_bin && start_bin >= 1 && end_bin <= size(depths_metaData, 1)
                        mean_A = mean(depths_metaData[start_bin:end_bin, 1])
                        mean_B = mean(depths_metaData[start_bin:end_bin, 2])
                        push!(mean_depth_A, mean_A)
                        push!(mean_depth_B, mean_B)
                        push!(depth_diff, mean_B - mean_A)
                    else
                        push!(mean_depth_A, 0.0)
                        push!(mean_depth_B, 0.0)
                        push!(depth_diff, 0.0)
                    end
                end
                cnv_regionsDF_samplepair[!, :mean_depth_A] = mean_depth_A
                cnv_regionsDF_samplepair[!, :mean_depth_B] = mean_depth_B
                cnv_regionsDF_samplepair[!, :depth_diff] = depth_diff
            end
         
            # 记录所有染色体可绘制区域
            sample1 = repeat([sample_pair[pair_i,1]], size(cnv_regionsDF_samplepair, 1))
            sample2 = repeat([sample_pair[pair_i,2]], size(cnv_regionsDF_samplepair, 1))
            temp_DF = DataFrame(pair_sampleA = sample1, pair_sampleB = sample2)
            temp_DF = hcat(temp_DF, cnv_regionsDF_samplepair)
            plotcnv_regionDF = vcat(plotcnv_regionDF, temp_DF)
         
            # 挑选绘制片段
            if chr == "all"
                # 去除不合要求片段
                selected_idx = findall(fragments_chr_endposi .- fragments_chr_startposi .> 0)
                # 绘制所有染色体
                temp_x = 1:length(selected_idx)
                temp_y1 = depths_mean[selected_idx, 1]
                plot!(
                    plt,
                    subplot = pair_i,
                    temp_x,
                    temp_y1,
                    linewidth = 0.5,
                    color = "red",
                    legend = false,
                    label = "Sample $(sample_pair[pair_i, 1])"
                )
             
                temp_y2 = depths_mean[selected_idx, 2]
                plot!(
                    plt,
                    subplot = pair_i,
                    temp_x,
                    temp_y2,
                    linewidth = 0.5,
                    color = "blue",
                    legend = false,
                    label = "Sample $(sample_pair[pair_i, 2])"
                )
             
                # 如果有CNV区域数据，绘制CNV区域
                if size(cnv_regionsDF_samplepair, 1) > 0
                    # 这里需要根据实际数据结构调整
                    # 示例：绘制黑色线条表示CNV区域
                    # plot!(plt, subplot=pair_i, temp_xCNV', temp_yCNV', linewidth = 1, color = "black", label = "CNV Regions")
                end
            else
                # 绘制指定染色体
                selected_idx = findall((fragments_chr .== chr) .&& (fragments_chr_endposi .- fragments_chr_startposi .> 0))
                if !isempty(selected_idx)
                    temp_x = hcat(fragments_chr_startposi[selected_idx], fragments_chr_endposi[selected_idx])
                    temp_y1 = hcat(depths_mean[selected_idx, 1], depths_mean[selected_idx, 1])
                    plot!(
                        plt,
                        subplot = pair_i,
                        temp_x',
                        temp_y1',
                        linewidth = 1,
                        color = "red",
                        legend = false,
                        title = "Sample $(sample_pair[pair_i, :])"
                    )
                 
                    temp_y2 = hcat(depths_mean[selected_idx, 2], depths_mean[selected_idx, 2])
                    plot!(
                        plt,
                        subplot = pair_i,
                        temp_x',
                        temp_y2',
                        linewidth = 1,
                        color = "blue",
                        legend = false,
                    )
                 
                    # 绘制CNV区域
                    if isa(cnv_regionsDF_samplepair, DataFrame) && size(cnv_regionsDF_samplepair, 1) > 0 && size(cnv_regionsDF_samplepair, 2) >= 3
                        chr_cnv_rows = cnv_regionsDF_samplepair[:, :chr] .== chr
                        if any(chr_cnv_rows)
                            selected_cnv_idx = findall(chr_cnv_rows)
                            # 过滤掉 |depth_diff| < 15 的行
                            valid_mask = abs.(cnv_regionsDF_samplepair[selected_cnv_idx, :depth_diff]) .>= 15
                            selected_cnv_idx = selected_cnv_idx[valid_mask]
                            if !isempty(selected_cnv_idx)
                                temp_xCNV = hcat(
                                    cnv_regionsDF_samplepair[selected_cnv_idx, :start_pos],
                                    cnv_regionsDF_samplepair[selected_cnv_idx, :end_pos]
                                )
                                temp_yCNV = ones(length(selected_cnv_idx), 2) * 1 # 固定Y值用于可视化
                                plot!(
                                    plt,
                                    subplot = pair_i,
                                    temp_xCNV',
                                    temp_yCNV',
                                    linewidth = 2,
                                    color = "black",
                                    legend = false
                                )
                            end
                        end
                    end
                else
                    @warn "未找到染色体 $chr 的数据"
                end
            end
        end
     
        #调整显示范围
        for pair_i = 1:num_pairs #获取当前参数
            xlims_current[pair_i, :] = [xlims(plt[pair_i])[1] xlims(plt[pair_i])[2]]
            ylims_current[pair_i, :] = [ylims(plt[pair_i])[1] ylims(plt[pair_i])[2]]
        end
        min_xlims = minimum(xlims_current[:, 1]) #取极值
        min_ylims = minimum(ylims_current[:, 1])
        max_ylims = maximum(ylims_current[:, 2])
        max_xlims = maximum(xlims_current[:, 2])
     
        if xlims_pre == [0 0]
            for pair_i = 1:num_pairs #调整画面
                xlims!(plt[pair_i], min_xlims, max_xlims)
            end
        else
            for pair_i = 1:num_pairs #调整画面
                xlims!(plt[pair_i], xlims_pre[1], xlims_pre[2])
            end
        end
        if ylims_pre == [0 0]
            for pair_i = 1:num_pairs #调整画面
                ylims!(plt[pair_i], min_ylims, max_ylims)
            end
        else
            for pair_i = 1:num_pairs #调整画面
                ylims!(plt[pair_i], ylims_pre[1], ylims_pre[2])
            end
        end
         
        # 显示图形
        display(plt)
       
        # 创建过滤后的 subplotcnv_regionDF
        if hasproperty(plotcnv_regionDF, :depth_diff) && !isempty(plotcnv_regionDF)
            valid_filter = abs.(coalesce.(plotcnv_regionDF.depth_diff, 0.0)) .>= 20
            subplotcnv_regionDF = plotcnv_regionDF[valid_filter, :]
        else
            subplotcnv_regionDF = DataFrame()
        end
       
        return plotcnv_regionDF, subplotcnv_regionDF, plt
     
    catch e
        @error "在plotCNVregionMulti函数中发生错误: " * string(e)
        rethrow(e)
    end
end








"""
    plotCNVregion(metaDataDF, CNVscaffold, pairedCNV, sample_pair; positive_cnv_content = 0.8, chr = "all")

绘制CNV区域图

# 参数说明:
- `metaDataDF`: 经过过滤且标准化后的数据框架,来自函数filtering_normalizing_seqDepth()
- `CNVscarfold`: CNV片段框架矩阵
- `pairedCNV`: 配对样本的CNV鉴别情况,3D数组,来自函数copyNumber_comparison()
- `sample_pair`: 一个具有两个数值元素的向量，表示样本对，如[1,2]表示第一个样本对和第二个样本对
- `positive_cnv_content`: CNV区域内被判定为CNV的基因组片段覆盖率最低值，默认值为0.8
- `chr`: 要绘制的染色体，默认为"all"表示所有染色体

# 返回结果:
- 无返回值，直接显示图形

# 依赖库:
- Plots: 绘图库
"""
function plotCNVregion(metaDataDF, CNVscarfold, pairedCNV, sample_pair; positive_cnv_content = 0.8, chr = "all")
    
    try
        # 输入验证
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF必须是DataFrame类型"))
        end
        
        if !isa(CNVscarfold, Matrix) || size(CNVscarfold, 2) != 3
            throw(ArgumentError("CNVscarfold必须是n行3列的矩阵"))
        end               
        
        if !isa(pairedCNV, Array) || length(size(pairedCNV)) != 3
            throw(ArgumentError("pairedCNV必须是3D的数组"))
        end

        if !isa(sample_pair, Vector) || length(sample_pair) != 2
            throw(ArgumentError("sample_pair必须是长度为2的向量"))
        end  
        

        # 数据准备        
        cnv_regionsDF_samplepair = cnv_region(metaDataDF, CNVscarfold, pairedCNV[:, sample_pair[1], sample_pair[2]], positive_cnv_content)
        metaDataDF_samplepair = select(metaDataDF, 1, 2, 3, 3 + sample_pair[1], 3 + sample_pair[2])        
        depths_metaData = Matrix(metaDataDF_samplepair[:, 4:end])

        # 计算CNV框架片段深度均值
        num_fragments = size(CNVscarfold, 1)
        depths_mean = zeros(num_fragments, 2)
        fragments_chr = String[]
        fragments_chr_startposi = Int64[]
        fragments_chr_endposi = Int64[]
        
        for i = 1:num_fragments
            push!(fragments_chr, String(metaDataDF_samplepair[CNVscarfold[i, 1], 1]))
            push!(fragments_chr_startposi, Int64(metaDataDF_samplepair[CNVscarfold[i, 1], 2]))
            push!(fragments_chr_endposi, Int64(metaDataDF_samplepair[CNVscarfold[i, 2], 3]))
            # 添加边界检查
            start_idx = max(1, CNVscarfold[i, 1])
            end_idx = min(size(depths_metaData, 1), CNVscarfold[i, 2])
            if start_idx <= end_idx
                depths_mean[i, 1] = mean(depths_metaData[start_idx:end_idx, 1])
                depths_mean[i, 2] = mean(depths_metaData[start_idx:end_idx, 2])
            else
                depths_mean[i, 1] = 0
                depths_mean[i, 2] = 0
            end
        end

        # 创建绘图
        plt = plot()
        
        # 挑选绘制片段
        if chr == "all"
            # 去除不合要求片段
            selected_idx = findall(fragments_chr_endposi .- fragments_chr_startposi .> 0)

            # 绘制所有染色体
            temp_x = 1:length(selected_idx)
            temp_y1 = depths_mean[selected_idx, 1]
            plot!(plt, temp_x, temp_y1, linewidth = 0.5, color = "red", legend = false)
            temp_y2 = depths_mean[selected_idx, 2]
            plot!(plt, temp_x, temp_y2, linewidth = 0.5, color = "blue", legend = false)
            
            # 如果有CNV区域数据，绘制CNV区域
            if size(cnv_regionsDF_samplepair, 1) > 0
                # 这里需要根据实际数据结构调整
                # 示例：绘制红色线条表示CNV区域
                # plot!(plt, temp_xCNV', temp_yCNV', linewidth = 1, color = "red", label = "CNV Regions")
            end
        else
            # 绘制指定染色体
            selected_idx = findall((fragments_chr .== chr) .&& (fragments_chr_endposi .- fragments_chr_startposi .> 0))
            if !isempty(selected_idx)
                temp_x = hcat(fragments_chr_startposi[selected_idx], fragments_chr_endposi[selected_idx])
                temp_y1 = hcat(depths_mean[selected_idx, 1], depths_mean[selected_idx, 1])
                plot!(plt, temp_x', temp_y1', linewidth = 1, color = "red", legend = false)
                temp_y2 = hcat(depths_mean[selected_idx, 2], depths_mean[selected_idx, 2])
                plot!(plt, temp_x', temp_y2', linewidth = 1, color = "blue", legend = false)
                
                # 绘制CNV区域
                chr_cnv_rows = cnv_regionsDF_samplepair[:, 1] .== chr
                if any(chr_cnv_rows) && size(cnv_regionsDF_samplepair, 2) >= 3
                    selected_cnv_idx = findall(chr_cnv_rows)
                    temp_xCNV = hcat(
                        cnv_regionsDF_samplepair[selected_cnv_idx, 2],
                        cnv_regionsDF_samplepair[selected_cnv_idx, 3]
                    )
                    temp_yCNV = ones(length(selected_cnv_idx), 2) * 1  # 固定Y值用于可视化
                    plot!(plt, temp_xCNV', temp_yCNV', linewidth = 2, color = "black", legend = false )
                end
            else
                @warn "未找到染色体 $chr 的数据"
            end
        end
        
        # 显示图形
        display(plt)
        return plt
        
    catch e
        @error "在plotCNVregion函数中发生错误: " * string(e)
        rethrow(e)
    end
end


"""
    cnv_region(metaDataDF, CNVscarffold, fragment_status, threshold_proportion)

根据片段状态合成CNV区域

# 参数说明:
- `metaDataDF`: 经过过滤且标准化后的数据框架，来自函数filtering_normalizing_seqDepth()
- `CNVscarffold`: 片段信息矩阵，来自函数fragment_status_matrix()
- `fragment_status`: 片段CNV状态，元素对应CNVscarffold的行，可以来自copyNumber_comparison
- `threshold_proportion`: 片段CNV的含量阈值，例如 0.8

# 返回结果:
- `coordinate_cnv_regions`: CNV区域坐标信息

# 依赖函数:
- synth_cnv: CNV合成函数

# 异常处理:
- 如果输入参数类型不正确，抛出ArgumentError
- 如果参数维度不匹配，抛出DimensionMismatch
"""
function cnv_region(metaDataDF, CNVscarffold, fragment_status, threshold_proportion)
    try
        # 输入验证
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF必须是DataFrame类型"))
        end
        
        if !isa(CNVscarffold, Matrix) || size(CNVscarffold, 2) != 3
            throw(ArgumentError("CNVscarffold必须是n行3列的矩阵"))
        end
        
        if !isa(fragment_status, Vector)
            throw(ArgumentError("fragment_status必须是Vector类型"))
        end
        
        if !isa(threshold_proportion, Number) || threshold_proportion < 0 || threshold_proportion > 1
            throw(ArgumentError("threshold_proportion必须是[0,1]范围内的数值"))
        end
        
        # 检查DataFrame是否有足够的列数（至少3列：chr, binLowEdge, binUpEdge）
        if size(metaDataDF, 2) < 3
            throw(DimensionMismatch("metaDataDF至少需要3列数据"))
        end
        
        # 检查维度兼容性
        if length(fragment_status) != size(CNVscarffold, 1)
            throw(DimensionMismatch("fragment_status长度必须与CNVscarffold行数匹配"))
        end
        
        # 构建片段坐标矩阵
        # 第1列为染色体编号，第2列为起始位置，第3列为终止位置，第4列为片段长度
        fragment_coords = hcat(
            metaDataDF[CNVscarffold[:, 1], 1],  # 染色体编号
            CNVscarffold[:, 1],                 # 起始位置
            CNVscarffold[:, 2],                 # 终止位置
            CNVscarffold[:, 2] - CNVscarffold[:, 1] .+ 1  # 片段长度
        )
        
        # 使用synth_cnv函数合成CNV区域
        temp_regions = synth_cnv(fragment_coords, fragment_status, threshold_proportion)
        
        # 如果synth_cnv返回错误码，则直接返回
        if temp_regions == -1
            @error "CNV合成过程中发生错误"
            return []
        end
        

        # 提取CNV区域坐标信息
        chr = String[]
        start_pos = Int64[]
        end_pos = Int64[]
        num_regions = size(temp_regions, 1)
        start_pos_metaData = Int64[]
        end_pos_metaData = Int64[]         
        if !isempty(temp_regions)
            for i in 1:num_regions
                temp_chr = String(metaDataDF[CNVscarffold[temp_regions[i, 1], 1], 1])
                temp_start_pos = Int64(metaDataDF[CNVscarffold[temp_regions[i, 1], 1], 2])
                temp_end_pos = Int64(metaDataDF[CNVscarffold[temp_regions[i, 2], 2], 3])
                if temp_end_pos >= temp_start_pos # 确保索引不越界                  
                    push!(chr, temp_chr)
                    push!(start_pos, temp_start_pos)
                    push!(end_pos, temp_end_pos)
                    push!(start_pos_metaData, CNVscarffold[temp_regions[i, 1], 1])
                    push!(end_pos_metaData, CNVscarffold[temp_regions[i, 2], 2])
                end
            end
        end

        coordinate_cnv_regions = DataFrame()
        coordinate_cnv_regions[!,:chr] = chr
        coordinate_cnv_regions[!,:start_pos] = start_pos
        coordinate_cnv_regions[!,:end_pos] = end_pos
        coordinate_cnv_regions[!,:length] = end_pos .- start_pos
        coordinate_cnv_regions[!,:start_pos_metaData] = start_pos_metaData
        coordinate_cnv_regions[!,:end_pos_metaData] = end_pos_metaData

        return coordinate_cnv_regions
        
    catch e
        @error "在cnv_region函数中发生错误: " * string(e)
        rethrow(e)
    end
end


"""
    copyNumber_comparison(metaDataDF, CNVscaffold, copynumber_info, gmmParam, windowSize; pvalue_theSame = 0.05, fdr_threshold = 0.2)

鉴别样本间拷贝数差异

# 参数说明:
- `metaDataDF`: 经过过滤且标准化后的数据框架，来自函数filtering_normalizing_seqDepth()
- `CNVscaffold`: CNV考察区间框架矩阵，由genome_segment函数生成
- `copynumber_info`: 拷贝数信息矩阵，由copynumber_determination函数生成
- `gmmParam`: GMM参数数组，由fitting_gmm函数生成
- `windowSize`: 估计GMM参数时片段窗口大小
- `pvalue_theSame`: 拷贝数一致性检验的p值，默认值0.05
- `fdr_threshold`: FDR阈值，宽松假阳性，降低假阴性，默认值0.2

# 返回结果:
- `cnv_state_info`: 三维表，表示样本间CNV状态信息
    - 第1维: CNV片段
    - 第2维: 样本j
    - 第3维: 样本k

# 依赖函数:
- prob_copy_differenceM: 计算拷贝数差异概率的自定义函数
- upgrade_state_vector: 更新状态向量的自定义函数

# 异常处理:
- 如果输入参数类型不正确，抛出ArgumentError
- 如果参数维度不匹配，抛出DimensionMismatch
"""
function copyNumber_comparison(metaDataDF, CNVscaffold, copynumber_info, gmmParam, windowSize; pvalue_theSame = 0.05, fdr_threshold = 0.1)
    try
        # 输入验证
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF必须是DataFrame类型"))
        end
        
        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold必须是n行3列的矩阵"))
        end
        
        if !isa(copynumber_info, Array) || ndims(copynumber_info) != 3
            throw(ArgumentError("copynumber_info必须是三维数组"))
        end
        
        if !isa(gmmParam, Vector)
            throw(ArgumentError("gmmParam必须是Vector类型"))
        end
        
        if !isa(pvalue_theSame, Number) || pvalue_theSame <= 0 || pvalue_theSame >= 1
            throw(ArgumentError("pvalue_theSame必须是(0,1)范围内的数值"))
        end
        
        if !isa(fdr_threshold, Number) || fdr_threshold <= 0 || fdr_threshold >= 1
            throw(ArgumentError("fdr_threshold必须是(0,1)范围内的数值"))
        end
        
        # 检查DataFrame是否有足够的列数（至少4列：chr, binLowEdge, binUpEdge, depth数据）
        if size(metaDataDF, 2) < 4
            throw(DimensionMismatch("metaDataDF至少需要4列数据"))
        end
        
        # 检查维度兼容性
        num_CNVfragment = size(CNVscaffold, 1)
        num_sample = size(copynumber_info, 1)
        
        if size(copynumber_info, 2) != num_CNVfragment
            throw(DimensionMismatch("copynumber_info的第二维必须与CNVscaffold的行数匹配"))
        end
        
        if length(gmmParam) != num_sample
            throw(DimensionMismatch("gmmParam长度必须与样本数量匹配"))
        end
        
        # ****************鉴别样本间拷贝数差异*****************
        subData = Matrix(metaDataDF[:, 4:end])
        num_fragment = size(subData, 1)
        num_samples = size(subData, 2)
        
        # 检查CNVscaffold的索引是否在有效范围内
        if !isempty(CNVscaffold)
            max_index = maximum(CNVscaffold[:, [1, 2]])
            if max_index > num_fragment
                throw(BoundsError("CNVscaffold中的索引超出了subData的行数范围"))
            end
            
            min_index = minimum(CNVscaffold[:, [1, 2]])
            if min_index < 1
                throw(BoundsError("CNVscaffold中的索引不能小于1"))
            end
        end
    
        
        cnv_state_info = zeros(num_CNVfragment, num_sample, num_sample)  # 三维表
        cnv_prob_info = zeros(num_CNVfragment, num_sample, num_sample)   # 三维表
        
        
        for j in 1:num_sample
            for k in 1:num_sample
                for i in 1:num_CNVfragment                    
                    # 两个片段正确鉴别且拷贝数不同
                    if (copynumber_info[j, i, 1] > 0) & (copynumber_info[k, i, 1] > 0) & 
                       (copynumber_info[j, i, 1] != copynumber_info[k, i, 1])
                        try
                            temp_star = CNVscaffold[i, 1]
                            temp_end = CNVscaffold[i, 2]
                            # 确保索引在有效范围内
                            frag_start = max(1, temp_star)
                            frag_end = min(num_fragment, temp_end)

                            # 计算拷贝数差异概率
                            cnv_prob_info[i, j, k] = prob_copy_differenceM(
                                [subData[frag_start:frag_end, j] subData[frag_start:frag_end, k]], 
                                [copynumber_info[j, i, 1] copynumber_info[k, i, 1]], 
                                gmmParam[[j k]],
                                windowSize
                            )                            
                            # 拷贝数一致概率低则标记为存在差异
                            if (cnv_prob_info[i, j, k] < pvalue_theSame) & (copynumber_info[j, i, 1] != copynumber_info[k, i, 1])
                                cnv_state_info[i, j, k] = 1
                            end 
                        catch e
                            @warn "计算片段$i 样本 $j 与 $k 的拷贝数差异时出错: " * string(e)
                        end
                    else
                         cnv_prob_info[i, j, k] = 1
                    end                    
                end
                #println("样本 $j 与 $k 的拷贝数差异比较完成: ", sum(cnv_state_info[:, j, k]))
                #println("样本 $j 与 $k 的差异概率 = 0 数量: ", sum(cnv_prob_info[:, j, k] .== 0))
            end
        end
        
        # *************扩展拷贝数差异区域，宽松假阳性，降低假阴性*****************
        # 以下处理会造成不必要的困扰，暂时放弃，待完善
        
        for j in 1:num_sample
            for k in 1:num_sample
                try
                    #println("样本 $j 与 $k 的差异: ", sum((cnv_prob_info[:, j, k] .<= fdr_threshold) .& (cnv_state_info[:, j, k]  .== 0)))
                    #println("样本 $j 与 $k 的差异: ", sum((cnv_prob_info[:, j, k] .<= fdr_threshold))) 
                    #println("样本 $j 与 $k 的差异: ", sum(cnv_state_info[:, j, k]  .== 0)) 
                    #println("Size : ", size(cnv_state_info[:, j, k]))                     
                    while true
                        temp_state_info = upgrade_state_vector(cnv_prob_info[:, j, k], cnv_state_info[:, j, k], fdr_threshold)
                        #println("样本 $j 与 $k 的拷贝数差异: ", sum(temp_state_info .!= cnv_state_info[:, j, k]))
                        #println("样本 $j 与 $k 的从1到0的差异: ", sum((temp_state_info .== 1) .& (cnv_state_info[:, j, k] .== 0)))
                        #println("样本 $j 与 $k 的概率小于阈值的数量: ", sum(cnv_prob_info[:, j, k] .< fdr_threshold))
                        if sum(temp_state_info .!= cnv_state_info[:, j, k]) > 0                        
                            cnv_state_info[:, j, k] = temp_state_info
                        else
                            break
                        end
                    end
                catch e
                    @warn "更新样本$j 与$k 的状态向量时出错: " * string(e)
                end
            end
        end
                
        return cnv_state_info #, cnv_prob_info
        
    catch e
        @error "在copyNumber_comparison函数中发生错误: " * string(e)
        rethrow(e)
    end
end



function copyNumber_comparison_filter(metaDataDF, CNVscaffold, copynumber_info, gmmParam, windowSize; pvalue_theSame = 0.05, fdr_threshold = 0.1 , highRDthreshold = 100)
    try
        # 输入验证
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF必须是DataFrame类型"))
        end
        
        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold必须是n行3列的矩阵"))
        end
        
        if !isa(copynumber_info, Array) || ndims(copynumber_info) != 3
            throw(ArgumentError("copynumber_info必须是三维数组"))
        end
        
        if !isa(gmmParam, Vector)
            throw(ArgumentError("gmmParam必须是Vector类型"))
        end
        
        if !isa(pvalue_theSame, Number) || pvalue_theSame <= 0 || pvalue_theSame >= 1
            throw(ArgumentError("pvalue_theSame必须是(0,1)范围内的数值"))
        end
        
        if !isa(fdr_threshold, Number) || fdr_threshold <= 0 || fdr_threshold >= 1
            throw(ArgumentError("fdr_threshold必须是(0,1)范围内的数值"))
        end
        
        # 检查DataFrame是否有足够的列数（至少4列：chr, binLowEdge, binUpEdge, depth数据）
        if size(metaDataDF, 2) < 4
            throw(DimensionMismatch("metaDataDF至少需要4列数据"))
        end
        
        # 检查维度兼容性
        num_CNVfragment = size(CNVscaffold, 1)
        num_sample = size(copynumber_info, 1)
        
        if size(copynumber_info, 2) != num_CNVfragment
            throw(DimensionMismatch("copynumber_info的第二维必须与CNVscaffold的行数匹配"))
        end
        
        if length(gmmParam) != num_sample
            throw(DimensionMismatch("gmmParam长度必须与样本数量匹配"))
        end
        
        # ****************鉴别样本间拷贝数差异*****************
        subData = Matrix(metaDataDF[:, 4:end])
        num_fragment = size(subData, 1)
        num_samples = size(subData, 2)
        
        # 检查CNVscaffold的索引是否在有效范围内
        if !isempty(CNVscaffold)
            max_index = maximum(CNVscaffold[:, [1, 2]])
            if max_index > num_fragment
                throw(BoundsError("CNVscaffold中的索引超出了subData的行数范围"))
            end
            
            min_index = minimum(CNVscaffold[:, [1, 2]])
            if min_index < 1
                throw(BoundsError("CNVscaffold中的索引不能小于1"))
            end
        end
    
        
        cnv_state_info = zeros(num_CNVfragment, num_sample, num_sample)  # 三维表
        cnv_prob_info = zeros(num_CNVfragment, num_sample, num_sample)   # 三维表
        
        
        for j in 1:num_sample
            for k in 1:num_sample
                for i in 1:num_CNVfragment                    
                    # 两个片段正确鉴别且拷贝数同（仅保留概率判定，不直接用deltaRDmean阈值）
                    if (copynumber_info[j, i, 1] > 0) & (copynumber_info[k, i, 1] > 0) & 
                       (copynumber_info[j, i, 1] == copynumber_info[k, i, 1])
                        # 本分支不再直接设置cnv_state_info=1（移到后续扩展）
                        continue
                    elseif (copynumber_info[j, i, 1] > 0) & (copynumber_info[k, i, 1] > 0) & 
                       (copynumber_info[j, i, 1] != copynumber_info[k, i, 1])
                        try
                            temp_star = CNVscaffold[i, 1]
                            temp_end = CNVscaffold[i, 2]
                            # 确保索引在有效范围内
                            frag_start = max(1, temp_star)
                            frag_end = min(num_fragment, temp_end)

                            # 计算拷贝数差异概率（核心判定条件）
                            cnv_prob_info[i, j, k] = prob_copy_differenceM(
                                [subData[frag_start:frag_end, j] subData[frag_start:frag_end, k]], 
                                [copynumber_info[j, i, 1] copynumber_info[k, i, 1]], 
                                gmmParam[[j k]],
                                windowSize
                            )
                                                       
                            if (cnv_prob_info[i, j, k] < pvalue_theSame)
                                cnv_state_info[i, j, k] = 1
                            end 
                        catch e
                            @warn "计算片段$i 样本 $j 与 $k 的拷贝数差异时出错: " * string(e)
                        end
                    elseif (copynumber_info[j, i, 1] < 0) | (copynumber_info[k, i, 1] < 0)
                        # 本分支不再直接设置cnv_state_info=1（移到后续扩展）
                        continue
                    else
                         cnv_prob_info[i, j, k] = 1
                    end                    
                end
            end
        end
        
        # *************扩展拷贝数差异区域（宽松假阳性，降低假阴性）*****************
        # 先基于核心概率判定得到初始cnv_state_info
        # 再进行扩展：对于原始为0的片段，如果deltaRDmean > highRDthreshold 且相邻片段至少有一个为1，则扩展为1
        # 这里进行多次迭代（3次足够传播大多数连续区域），避免单次遍历不完全传播
        for iter in 1:3
            for j in 1:num_sample
                for k in 1:num_sample
                    for i in 1:num_CNVfragment
                        if cnv_state_info[i, j, k] == 0  # 只扩展原始为0的
                            temp_star = CNVscaffold[i, 1]
                            temp_end = CNVscaffold[i, 2]
                            frag_start = max(1, temp_star)
                            frag_end = min(num_fragment, temp_end)
                            
                            deltaRDmean = abs(mean(subData[frag_start:frag_end, j]) - mean(subData[frag_start:frag_end, k]))
                            
                            if deltaRDmean > highRDthreshold
                                has_neighbor_diff = false
                                if i > 1 && cnv_state_info[i-1, j, k] == 1
                                    has_neighbor_diff = true
                                end
                                if i < num_CNVfragment && cnv_state_info[i+1, j, k] == 1
                                    has_neighbor_diff = true
                                end
                                
                                if has_neighbor_diff
                                    cnv_state_info[i, j, k] = 1
                                end
                            end
                        end
                    end
                end
            end
        end
                
        return cnv_state_info #, cnv_prob_info
        
    catch e
        @error "在copyNumber_comparison函数中发生错误: " * string(e)
        rethrow(e)
    end
end


#先进行拷贝数差异鉴别，再运行以下函数，减少RD假阴性
#用法：pairedCNVT293FN = copyNumber_comparison_filter_fn(metaDataDF, T293CNVscaffold, copynumber_info, pairedCNVT293;  highRDthreshold =20);
#pairedCNVT293FN为初步鉴别配对样本拷贝数差异数据pairedCNVT293的补充版本
#增添的区域R:紧邻阳性区域，且RD差异大于阈值(不管拷贝数鉴定是否相同、不同、还是-1都包含在内)
#循环迭代几轮，直至无法延伸
function copyNumber_comparison_filter_fn(metaDataDF, CNVscaffold, copynumber_info, cnv_state_info;   highRDthreshold = 100)
    try
        # 输入验证
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF必须是DataFrame类型"))
        end
        
        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold必须是n行3列的矩阵"))
        end
        
        if !isa(copynumber_info, Array) || ndims(copynumber_info) != 3
            throw(ArgumentError("copynumber_info必须是三维数组"))
        end
       
        if !isa(cnv_state_info, Array) || ndims(copynumber_info) != 3
            throw(ArgumentError("copynumber_info必须是三维数组"))
        end
        
        
        # 检查DataFrame是否有足够的列数（至少4列：chr, binLowEdge, binUpEdge, depth数据）
        if size(metaDataDF, 2) < 4
            throw(DimensionMismatch("metaDataDF至少需要4列数据"))
        end
        
        # 检查维度兼容性
        num_CNVfragment = size(CNVscaffold, 1)
        num_sample = size(copynumber_info, 1)
        
        if size(copynumber_info, 2) != num_CNVfragment
            throw(DimensionMismatch("copynumber_info的第二维必须与CNVscaffold的行数匹配"))
        end
        
    
        
        # ****************鉴别样本间拷贝数差异*****************
        subData = Matrix(metaDataDF[:, 4:end])
        num_fragment = size(subData, 1)
        num_samples = size(subData, 2)
        
        # 检查CNVscaffold的索引是否在有效范围内
        if !isempty(CNVscaffold)
            max_index = maximum(CNVscaffold[:, [1, 2]])
            if max_index > num_fragment
                throw(BoundsError("CNVscaffold中的索引超出了subData的行数范围"))
            end
            
            min_index = minimum(CNVscaffold[:, [1, 2]])
            if min_index < 1
                throw(BoundsError("CNVscaffold中的索引不能小于1"))
            end
        end
    
        
        # 初始化：第一轮的输入为原始cnv_state_info
        cnv_state_info_prev = copy(cnv_state_info)
        # 标记是否继续循环（初始为true，至少执行一轮）
        has_change = true
        # 记录迭代轮数（可选，用于监控）
        iter_num = 0

        # 循环迭代：只要有更新就继续
        while has_change
            iter_num += 1
            @info "开始第 $iter_num 轮迭代"
            # 初始化本轮的新状态（基于上一轮的结果）
            cnv_state_info_new = copy(cnv_state_info_prev)
            # 标记本轮是否有更新
            current_change = false
            # 核心逻辑循环（与原代码一致，仅替换cnv_state_info为cnv_state_info_prev）
            for j in 1:num_samples
                for k in 1:num_samples 
                    for i in 2:(num_CNVfragment-1)
                        # 两个片段正确鉴别且拷贝数同
                        if (copynumber_info[j, i, 1] > 0) & (copynumber_info[k, i, 1] > 0)
                            try
                                temp_star = CNVscaffold[i, 1]
                                temp_end = CNVscaffold[i, 2]
                                # 确保索引在有效范围内
                                frag_start = max(1, temp_star)
                                frag_end = min(num_fragment, temp_end)
            
                                deltaRDmean = abs(mean(subData[frag_start:frag_end, j]) - mean(subData[frag_start:frag_end, k]))
                                if (deltaRDmean .> highRDthreshold) & ((cnv_state_info_prev[i-1, j, k] == 1) | (cnv_state_info_prev[i+1, j, k] == 1))
                                # 检查是否与原值不同，若不同则标记更新
                                   if (cnv_state_info_new[i, j, k] != 1)
                                       cnv_state_info_new[i, j, k] = 1
                                       current_change = true
                                   end
                                end 
                            catch e
                                @warn "计算片段$i 样本 $j 与 $k 的拷贝数差异时出错: " * string(e)
                            end
                
                        elseif (copynumber_info[j, i, 1] < 0) | (copynumber_info[k, i, 1] < 0)
                            try
                                temp_star = CNVscaffold[i, 1]
                                temp_end = CNVscaffold[i, 2]
                                # 确保索引在有效范围内
                                frag_start = max(1, temp_star)
                                frag_end = min(num_fragment, temp_end)
                                deltaRDmean = abs(mean(subData[frag_start:frag_end, j]) - mean(subData[frag_start:frag_end, k]))
                                if (deltaRDmean .> highRDthreshold) & 
                                    ((cnv_state_info_prev[i-1, j, k] == 1) | (cnv_state_info_prev[i+1, j, k] == 1))
                                # 检查是否与原值不同，若不同则标记更新
                                    if (cnv_state_info_new[i, j, k] != 1)
                                        cnv_state_info_new[i, j, k] = 1
                                        current_change = true
                                    end
                                end
                            catch e
                                @warn "计算片段$i 样本 $j 与 $k 的拷贝数差异时出错: " * string(e)
                            end
                        end                    
                    end
        
                end
            end

    # 更新：将本轮结果作为下一轮的输入
            cnv_state_info_prev = copy(cnv_state_info_new)
    # 更新循环标记：本轮有更新则继续，否则终止
            has_change = current_change

            @info "第 $iter_num 轮迭代完成，是否有更新: $has_change"
        end

        # 最终结果存储在 cnv_state_info_prev 中（或 cnv_state_info_new）
        @info "迭代结束，总迭代轮数: $iter_num"
        final_cnv_state_info = cnv_state_info_prev
        
        # *************扩展拷贝数差异区域，宽松假阳性，降低假阴性*****************
        # 以下处理会造成不必要的困扰，暂时放弃，待完善
        
                
        return final_cnv_state_info #, cnv_prob_info
        
    catch e
        @error "在copyNumber_comparison函数中发生错误: " * string(e)
        rethrow(e)
    end
end


"""
    copynumber_determination(metaDataDF, CNVscaffold, gmmParam, flanking_length; thresholdBFactor = 2.0)

识别基因组片段的拷贝数

# 参数说明:
- `metaDataDF`: 经过过滤且标准化后的数据框架，来自函数filtering_normalizing_seqDepth()
- `CNVscaffold`: CNV考察区间框架矩阵，由genome_segment函数生成
- `gmmParam`: 模型参数，来自函数fitting_gmm()
- `flanking_length`: 片段的左右边界长度，与调用genome_segment()时保持一致   
- `thresholdBFactor`: 贝叶斯因子的阈值，默认典型值为2.0

# 返回结果:
- `copyNumberInfo`: 拷贝数信息，三维表，1D样品，2D片段，3D拷贝数的详细信息
  拷贝数详细信息：第一列为拷贝数，后续各列为拷贝数为k-1的概率

# 依赖函数:
- identify_copy_number: 识别拷贝数的自定义函数

# 异常处理:
- 如果输入参数类型不正确，抛出ArgumentError
- 如果参数维度不匹配，抛出DimensionMismatch
- 如果索引越界，抛出BoundsError
"""
function copynumber_determination(metaDataDF, CNVscaffold, gmmParam, flanking_length; thresholdBFactor = 2.0)
    try
        # 输入验证
        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF必须是DataFrame类型"))
        end
        
        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold必须是n行3列的矩阵"))
        end
        
        if !isa(gmmParam, Vector)
            throw(ArgumentError("gmmParam必须是Vector类型"))
        end
        
        if !isa(flanking_length, Integer) || flanking_length <= 0
            throw(ArgumentError("flanking_length必须是正整数"))
        end
        
        if !isa(thresholdBFactor, Number) || thresholdBFactor <= 0
            throw(ArgumentError("thresholdBFactor必须是正数"))
        end
        
        # 检查DataFrame是否有足够的列数（至少4列：chr, binLowEdge, binUpEdge, depth数据）
        if size(metaDataDF, 2) < 4
            throw(DimensionMismatch("metaDataDF至少需要4列数据"))
        end
        
        # ***************确定片段拷贝数*******************
        subData = Matrix(metaDataDF[:, 4:end])
        num_row = size(subData, 1)
        num_col = size(subData, 2)
        
        # 检查参数兼容性
        if flanking_length * 2 >= num_row
            throw(ArgumentError("flanking_length过大，与数据长度不兼容"))
        end
        
        if length(gmmParam) != num_col
            throw(DimensionMismatch("gmmParam长度必须与样本数量匹配"))
        end
        
        # 检查CNVscaffold的索引是否在有效范围内
        if !isempty(CNVscaffold)
            max_index = maximum(CNVscaffold[:, [1, 2]])
            if max_index > num_row
                throw(BoundsError("CNVscaffold中的索引超出了subData的行数范围"))
            end
            
            min_index = minimum(CNVscaffold[:, [1, 2]])
            if min_index < 1
                throw(BoundsError("CNVscaffold中的索引不能小于1"))
            end
        end
        
        AllMeanSub = zeros(size(subData))
        for i in (flanking_length + 1):(num_row - flanking_length)
            for j in 1:num_col
                start_idx = max(1, i - flanking_length)
                end_idx = min(num_row, i + flanking_length)
                AllMeanSub[i, j] = mean(subData[start_idx:end_idx, j])
            end
        end
        
        depthFrag = zeros(size(CNVscaffold, 1), num_col)
        for i in 1:size(CNVscaffold, 1)
            for j in 1:num_col
                temp_star = CNVscaffold[i, 1]
                temp_end = CNVscaffold[i, 2]
                # 确保索引在有效范围内
                frag_start = max(1, temp_star)
                frag_end = min(num_row, temp_end)
                if frag_start <= frag_end
                    depthFrag[i, j] = mean(AllMeanSub[frag_start:frag_end, j])
                else
                    depthFrag[i, j] = 0.0
                end
            end
        end
        
        # 三维表，1D样品，2D片段，3D拷贝数的详细信息
        copyNumberInfo = zeros(num_col, size(CNVscaffold, 1), length(gmmParam[1].means) + 1)
        for i in 1:num_col
            for j in 1:size(CNVscaffold, 1)
                try
                    # 使用自写函数识别拷贝数
                    copyNumberInfo[i, j, :] = identify_copy_number(depthFrag[j, i], gmmParam[i], thresholdBFactor)
                    # 解决0拷贝问题
                    if depthFrag[j, i] < (gmmParam[i].means[1] - 2 * gmmParam[i].stds[1])
                        copyNumberInfo[i, j, 1] = 0
                    end
                catch e
                    @warn "样本 $i, 片段 $j 的拷贝数识别失败: " * string(e)
                    copyNumberInfo[i, j, 1] = -1  # 标记为无法确定
                end
            end
        end
        
        return copyNumberInfo
        
    catch e
        @error "在copynumber_determination函数中发生错误: " * string(e)
        rethrow(e)
    end
end


"""
    genome_segment(inputDataDF, givenDistParam, breaks_threshold, flanking_length, minimum_length)

识别基因组片段断点并构建CNV考察区间

# 参数说明:
- `inputDataDF`: 经过过滤且标准化后的数据框架，来自函数filtering_normalizing_seqDepth()
- `givenDistParam`: 存储所有基因片段的分布参数，来自函数fitting_gmm()
- `breaks_threshold`: 断点阈值，值域[0, 1]，例如0.001
- `flanking_length`: 被考察片段左右侧翼边界长度，应与GMM拟合时片段长度一致，例如10
- `minimum_length`: 最小片段长度，例如10

# 依赖函数:
- rate_coincident: 计算数据与分布参数的一致性概率
- run_length_nonzero: 计算非零游程

# 返回值:
- `CNVscaffold`: CNV考察区间框架矩阵

# 异常处理:
- 如果输入参数类型不正确，抛出ArgumentError
- 如果参数值超出有效范围，抛出ArgumentError
"""
function genome_segment(inputDataDF, givenDistParam, breaks_threshold, flanking_length, minimum_length)
    try
        # 输入验证
        if !isa(inputDataDF, DataFrame)
            throw(ArgumentError("inputDataDF必须是DataFrame类型"))
        end
        
        if !isa(givenDistParam, Vector)
            throw(ArgumentError("givenDistParam必须是Vector类型"))
        end
        
        if !isa(breaks_threshold, Number) || breaks_threshold < 0 || breaks_threshold > 1
            throw(ArgumentError("breaks_threshold必须是[0, 1]范围内的数值"))
        end
        
        if !isa(flanking_length, Integer) || flanking_length <= 0
            throw(ArgumentError("flanking_length必须是正整数"))
        end
        
        if !isa(minimum_length, Integer) || minimum_length <= 0
            throw(ArgumentError("minimum_length必须是正整数"))
        end
        
        # 检查DataFrame是否有足够的列数（至少4列：chr, binLowEdge, binUpEdge, depth数据）
        if size(inputDataDF, 2) < 4
            throw(DimensionMismatch("inputDataDF至少需要4列数据"))
        end
        
        # ***************找出所有候选断点(CNV片段交界)*******************
        subData = Matrix(inputDataDF[:, 4:end])
        row_num = size(subData, 1)
        sample_num = size(subData, 2)
        
        # 检查参数兼容性
        if flanking_length * 2 >= row_num
            throw(ArgumentError("flanking_length过大，与数据长度不兼容"))
        end
        
        if length(givenDistParam) != sample_num
            throw(DimensionMismatch("givenDistParam长度必须与样本数量匹配"))
        end
        
        LeftMeanSub = zeros(size(subData))
        RightMeanSub = zeros(size(subData))
        for i in (flanking_length + 1):(row_num - flanking_length)
            for j in 1:sample_num
                LeftMeanSub[i, j] = mean(subData[(i - flanking_length):(i - 1), j])
                RightMeanSub[i, j] = mean(subData[(i + 1):(i + flanking_length), j])
            end
        end
        
        temp_breaks = zeros(size(subData))
        for i in (flanking_length + 1):(row_num - flanking_length)
            for j in 1:sample_num
                temp = [LeftMeanSub[i, j]; RightMeanSub[i, j]]
                # 通过自编函数返回的概率确定候选断点
                if (rate_coincident(temp, givenDistParam[j]) < breaks_threshold) || 
                    (inputDataDF[i + flanking_length, 3] - inputDataDF[i - flanking_length, 3] > 100*10^3) # 测序深度差异大或者物理间隔过长，设置为断点
                    temp_breaks[i, j] = 1
                end
            end
        end
        
        temp = sum(temp_breaks, dims = 2) # 每个片段在各个样本的累积断裂次数
        
        # **************考虑断点间游程********************
        temp_runs = run_length_nonzero(temp) # 合并所有样本的断点
        
        # 将相邻断裂块的区间作为基本考量单位，区间包括边界在内
        # 计算区间
        if size(temp_runs, 1) <= 1
            @warn "未检测到足够的断点游程"
            return zeros(0, 3)
        end
        
        temp_runs_region = zeros(size(temp_runs, 1) - 1, 3)
        for i in 1:(size(temp_runs, 1) - 1)
            temp_star = temp_runs[i, 3]
            temp_end = temp_runs[i + 1, 2]
            temp_len = temp_end - temp_star
            temp_runs_region[i, :] = [temp_star temp_end temp_len] 
        end
        
        # **将片段长度大于minimum_length的区间挑选出来，作为最终CNV考察区间****
        temp_indices = findall(temp_runs_region[:, 3] .>= minimum_length)
        
        if isempty(temp_indices)
            @warn "未找到满足最小长度要求的CNV区间"
            return zeros(0, 3)
        end
        
        # 滤除跨染色体的片段，形成基本考察框架
        tempCNVscaffold = Int64.(temp_runs_region[temp_indices, :]) 
        temp_indices = findall(inputDataDF[tempCNVscaffold[:, 1], 1] .== inputDataDF[tempCNVscaffold[:, 2], 1])
        CNVscaffold = tempCNVscaffold[temp_indices, :]
        
        return CNVscaffold
        
    catch e
        @error "在genome_segment函数中发生错误: " * string(e)
        rethrow(e)
    end
end


"""
    fitting_gmm(inputData, pre_mean, windowSize)

使用高斯混合模型(GMM)拟合数据

# 参数说明:
- `inputData`: M×N 矩阵，N 是样本数量
- `pre_mean`: 高斯混合模型中每个组件的预定义均值
- `windowSize`: 窗口大小

# 返回值:
- `EstimateParam`: 包含每个样本GMM参数估计结果的数组

# 异常处理:
- 如果输入参数类型不正确，抛出ArgumentError
- 如果inputData维度与windowSize不匹配，抛出DimensionMismatch
"""
function fitting_gmm(inputData, pre_mean, windowSize)
    try
        # 输入验证
        if !isa(inputData, Matrix)
            throw(ArgumentError("inputData必须是Matrix类型"))
        end
        
        if !isa(pre_mean, Vector)
            throw(ArgumentError("pre_mean必须是Vector类型"))
        end
        
        if !isa(windowSize, Integer) || windowSize <= 0
            throw(ArgumentError("windowSize必须是正整数"))
        end
        
        # 检查维度兼容性
        M, N = size(inputData)
        if windowSize > M
            throw(DimensionMismatch("windowSize不能大于inputData的行数"))
        end
        
        # 初始化窗口均值矩阵
        windowMean = zeros(M - windowSize + 1, N)
        
        # 计算滑动窗口均值
        for i in 1:(M - windowSize + 1)
            for j in 1:N
                # 确保索引不越界
                startIdx = max(1, i)
                endIdx = min(M, i + windowSize - 1)
                windowMean[i, j] = mean(inputData[startIdx:endIdx, j])
            end
        end
        
        # 获取组件数量
        num_components = length(pre_mean)
        
        # 如果组件数量不匹配预定义均值数量，给出警告
        if num_components <= 0
            @warn "预定义均值数量为零"
            return []
        end
        
        # 参数估计
        EstimateParam = []
        for i in 1:N
            try
                # 使用estimate_gmm函数进行参数估计
                gmm_result = estimate_gmm(windowMean[:, i], num_components; means = pre_mean)
                EstimateParam = [EstimateParam; gmm_result]
                print("处理样本 $i/$N\n")
            catch e
                @warn "样本 $i 的参数估计失败: " * string(e)
            end
        end
        
        return EstimateParam
        
    catch e
        @error "在fitting_gmm函数中发生错误: " * string(e)
        rethrow(e)
    end
end


"""
    filtering_normalizing_seqDepth(dataframe_raw_sequencingdepth, depth_filter1, depth_filter2)

过滤和标准化测序深度数据

此函数的主要目的是过滤掉由于各种原因导致测序深度过低或者过高的片段，
这些片段的数据可靠性存疑，并对剩余数据进行标准化处理。

# 参数说明:
- `dataframe_raw_sequencingdepth`: DataFrame对象，包含测序深度信息
    - 第1列: 染色体编号(chr)
    - 第2列: 区间起始位置(binLowEdge)
    - 第3列: 区间终止位置(binUpEdge)
    - 第4列及之后: 每个样本在该区间的覆盖深度
- `depth_filter1`: 长度为2的向量，分别表示标注化前允许的最小和最大总深度阈值
    - 第1个元素: 最小深度阈值
    - 第2个元素: 最大深度阈值
- `depth_filter1`: 长度为2的向量，分别表示标注化后允许的最小和最大总深度阈值
    - 第1个元素: 最小深度阈值
    - 第2个元素: 最大深度阈值

# 返回值:
- `normalSubCombineData`: 标准化并过滤后的测序深度数据

# 异常处理:
- 如果输入参数类型不正确，抛出ArgumentError
- 如果depth_filter长度不为2，抛出DimensionMismatch
"""
function filtering_normalizing_seqDepth(dataframe_raw_sequencingdepth, depth_filter1, depth_filter2)
    try
        # 输入验证
        if !isa(dataframe_raw_sequencingdepth, DataFrame)
            throw(ArgumentError("dataframe_raw_sequencingdepth必须是DataFrame类型"))
        end
        
        if (!isa(depth_filter1, Vector) || length(depth_filter1) != 2) && (isa(depth_filter2, Vector) || length(depth_filter2) != 2)
            throw(DimensionMismatch("depth_filter必须是长度为2的向量"))
        end
        
        if depth_filter1[1] >= depth_filter1[2]
            throw(ArgumentError("depth_filter1中的最小深度必须小于最大深度"))
        end

        if depth_filter2[1] >= depth_filter2[2]
            throw(ArgumentError("depth_filter2中的最小深度必须小于最大深度"))
        end
        
        # 检查DataFrame是否有足够的列数
        if size(dataframe_raw_sequencingdepth, 2) < 4
            throw(DimensionMismatch("dataframe_raw_sequencingdepth至少需要4列数据"))
        end
        
        # 查看总覆盖深度
        tempDepth = sum(Matrix(dataframe_raw_sequencingdepth[:, 4:end]), dims = 2)
        
        # 仅提取测序深度在合适范围的区域
        # 提取总测序深度在指定范围内的片段
        tempidx = findall((tempDepth[:] .> depth_filter1[1]) .& (tempDepth[:] .< depth_filter1[2]) .== 1)
        
        # 如果没有满足条件的行，返回空DataFrame
        if isempty(tempidx)
            @warn "没有满足深度过滤条件的片段"
            return DataFrame()
        end
        
        subCombineData = dataframe_raw_sequencingdepth[tempidx, :]
        
        # 将平均测序深度标准化到100
        tempMean = mean(Matrix(subCombineData[:, 4:end]), dims = 1)
        normalSubCombineData = hcat(subCombineData[:, 1:3], Float64.(subCombineData[:, 4:end]))
        
        # 对每一行进行标准化处理
        for i in 1:size(normalSubCombineData, 1)
            # 防止除零错误
            if all(tempMean .> 0)
                normalSubCombineData[i, 4:end] = (Vector{Float64}(normalSubCombineData[i, 4:end])' ./ tempMean) .* 100
            else
                @warn "平均深度包含零值，跳过标准化"
                break
            end
        end

        #取出每行最小最大值满足要求的的子集并返回
        miniDepth = minimum(Matrix(normalSubCombineData[:,4:end]), dims = 2)
        maxiDepth = maximum(Matrix(normalSubCombineData[:,4:end]), dims = 2)
        tempidx = findall((miniDepth[:] .> depth_filter2[1]) .& (maxiDepth[:] .< depth_filter2[2]))
        
        return normalSubCombineData[tempidx,:]
        
    catch e
        @error "在filtering_normalizing_seqDepth函数中发生错误: " * string(e)
        rethrow(e)
    end
end
