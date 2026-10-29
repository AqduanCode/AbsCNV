# 最后一次修订时间: 2025-12-10 19:12

function synth_cnv(fragment_coords, fragment_status, threshold_proportion)
    # CNV合成函数，将CNV片段合成为CNV区域
    #参数:
    #- fragment_coords: n行4列的数据框架
    #    第1列: 染色体（字符）编号
    #    第2列: 当前scarffold片段起始位置
    #    第3列: 当前scarffold片段终止位置
    #    第4列: 当前scarffold片段长度
    #- fragment_status: 长度为n的数值{0,1}向量
    #    0表示当前scarffold片段为阴性
    #    1表示当前scarffold片段为阳性片段
    #    -1表示当前scarffold片段为状态不定

    #返回值:
    #- m行2列的矩阵
    #    第1列: CNV区域的起始片段位置在输入参数矩阵中的行号
    #    第2列: CNV区域的终止片段位置在输入参数矩阵中的行号
    #- 如果发生错误，返回-1


    # 步骤1: 检查参数1是否已按照染色体编号和起始位置排序
    # 检查染色体编号是否有序
    chromosomes = fragment_coords[:, 1]
    start_positions = fragment_coords[:, 2]
    
    # 检查是否按染色体和起始位置排序
    for i in 2:size(fragment_coords, 1)
        # 如果染色体编号变小，则未排序
        if chromosomes[i] < chromosomes[i-1]
            println("错误: CNV片段坐标矩阵未按染色体编号和起始位置排序")
            return -1
        # 如果在同一染色体上但起始位置变小，则未排序
        elseif chromosomes[i] == chromosomes[i-1] && start_positions[i] < start_positions[i-1]
            println("错误: CNV片段坐标矩阵未按染色体编号和起始位置排序")
            return -1
        end
    end
    
    # 步骤2: 从参数1（scarffold片段位置信息）中挑选出参数2（拷贝差异状态）值为1的元素对应的行号，组成CNV表
    # 获取阳性片段(值为1)的行号
    positive_indices = findall(fragment_status .== 1)
    
    # 如果没有阳性片段，直接返回空矩阵
    if isempty(positive_indices)
        println("错误: 无阳性片段")
        return -1
    end

    # 步骤3: 根据拟合并后的阳性片段含量进行CNV区域合并    
    # CNV表初始化，第1列和第2列都是阳性片段的行号
    #= 用一个k*2的矩阵模拟一个具有删除功能的链表，记录片段合并。
    第一列为片段起始位置，第二列为片段终止位置。
    片段i与i+1合并时（进行删除操作），将第i行第2列改为-1，i+1行第一列改为-1。
    对于整个矩阵，第一列第j个非-1元素记录的起始位置，始终与第二列第j个非-1元素形成逻辑上的(start_posi, end_posi)配对，记录一个片段。
    输出最终结果时进行重整，将各列所有非-1元素写入对应列，形成m*2输出矩阵。
    =#
    cnv_table = hcat(positive_indices, positive_indices)
    
    k = size(cnv_table, 1)
    if k <= 1
        # 只有一个CNV片段，无需合并
        return cnv_table
    end    

    # 评估拟合并后阳性片段含量，指导进行合并
    # 每轮评估中，满足合并要求的，全部进行合并
    # 进行多轮评估，直到无法再进行合并（没有满足要求的片段）

    merged = true
    while merged # 反复多轮评估，直到无法再进行合并
        merged = false
        #启动一通评估
        for i in 2:size(cnv_table, 1) 
            # 在原始CNV表中找到这两行代表的片段索引
            last_idx_fragment_coords = cnv_table[i-1, 1]
            current_idx_fragment_coords = cnv_table[i, 2]            

            # 计算CNV阳性片段含量
            temp_lengths = fragment_coords[last_idx_fragment_coords:current_idx_fragment_coords, 4]
            possible_proportion = sum(temp_lengths .* (fragment_status[last_idx_fragment_coords:current_idx_fragment_coords] .== 1)) / 
                                  sum(temp_lengths)           

            if (possible_proportion >= threshold_proportion) & 
                (fragment_coords[last_idx_fragment_coords, 1] == fragment_coords[current_idx_fragment_coords, 1])
                # 满足合并要求，进行合并
                cnv_table[i, 1] = cnv_table[i - 1, 1] #合并
                cnv_table[i - 1, 1] = -1 # 删除
                cnv_table[i - 1, 2] = -1 # 删除
                merged = true #如进行了合并，则需要进行下一通评估
            end            
        end
        #每通评估与合并完成后，升级 cnv_table
        temp_cnvtable = hcat(cnv_table[findall(cnv_table[:, 1] .> 0), 1], cnv_table[findall(cnv_table[:, 2] .> 0), 2])
        cnv_table = temp_cnvtable
    end     
    
    return cnv_table
end


function upgrade_state_vector(prob_vector, 
                             state_vector, 
                             fdr_threshold)
# 参数
#- `prob_vector::Vector{Float64}`: 概率向量，每个元素为[0,1]范围内的概率值
#- `state_vector::Vector{Int}`: 状态向量，取值为0（阴性）或1（阳性）
#- `fdr_threshold::Float64`: FDR控制阈值，如0.05
# 返回
#- `Vector{Int}`: 升级后的状态向量
# 功能流程
#1. 识别与阳性状态相邻的阴性位置索引
#2. 提取对应位置的概率值
#3. 使用Benjamini-Hochberg方法进行FDR校正
#4. 根据显著性结果更新状态向量
    
    # 参数验证
    @assert length(prob_vector) == length(state_vector) #"概率向量和状态向量长度必须一致"
    @assert all(0 .<= prob_vector .<= 1) #"概率值必须在[0,1]范围内"
    @assert all(state_vector .== 0 .|| state_vector .== 1) #"状态向量值必须为0或1"
    @assert 0 <= fdr_threshold <= 1 #"FDR阈值必须在[0,1]范围内"

    
    
    n = length(state_vector)
    adjacent_zero_indices = Int[]
    
    # 步骤1: 查找与阳性状态相邻的阴性位置
    for i in 1:n
        if state_vector[i] == 1
            # 检查左侧邻居
            if i > 1 && state_vector[i-1] == 0
                push!(adjacent_zero_indices, i-1)
            end
            # 检查右侧邻居
            if i < n && state_vector[i+1] == 0
                push!(adjacent_zero_indices, i+1)
            end
        end
    end

    

    # 去除重复索引
    adjacent_zero_indices = unique(adjacent_zero_indices)
    #println("相邻阴性位置索引长度: ", length(adjacent_zero_indices))
   
    
    # 步骤2: 提取相邻阴性位置的概率值
    if isempty(adjacent_zero_indices)
        return state_vector  # 无相邻阴性状态时直接返回
    end
    
    adjacent_probs = prob_vector[adjacent_zero_indices]    

    # 步骤3: FDR校正和显著性判断
    adjusted_pvalues = adjust(adjacent_probs, BenjaminiHochberg())
    
    significant_indices = findall(adjusted_pvalues .< fdr_threshold)
    #println("相邻显著阴性位置索引长度: ", length(significant_indices)) 
    #println("相邻显著阴性位置概率和: ",  -log10.(adjacent_probs'))
    #println("差异: ", sum((state_vector[adjacent_zero_indices]  .== 0)))

    if isempty(significant_indices)
        return state_vector  # 无满足阈值要求的相邻阴性状态时直接返回
    end
    
    # 更新状态向量
    for i in 1:length(significant_indices)
        original_idx = adjacent_zero_indices[significant_indices[i]]
        state_vector[original_idx] = 1
    end
    
    # 步骤4: 返回升级后的状态向量
    return state_vector
end


#***********************提升拷贝数变化的检验****************************
function prob_copy_differenceM(frag_depths, num_copies, givenParameters, windowSize)
    # frag_depth 是一个n行2列的矩阵，每列分别表示样本1和样本2该区域各短片段的测序深度
    # num_copies 是样本1和样本2该片段的各自拷贝数
    # givenParameters 是估计得到的参数，是样本1和样本2的estimate_gmm()的返回值
    # windowSize 是估计参数givenParameters时用的片段长度
    # 返回，片段在两个样本中不存在拷贝数变化的概率

    #依赖库 Distributions.jl

    #组分数
    len_frag = size(frag_depths)[1]
    num_copies = Int64.(num_copies)
    #pesudo_len = len_frag/windowSize
    pesudo_len = sqrt(len_frag/windowSize)
    #pesudo_len = 1

    meanwithin =  mean(frag_depths, dims = 1) # 1*2矩阵，列是不同样本
    meantotal = [givenParameters[1].means[num_copies]; givenParameters[2].means[num_copies]] # 2*2矩阵，行是不同片段，列是不同样本
    varwithin =  [var(frag_depths[:,1]) var(frag_depths[:,2])] ./ windowSize # 1*2矩阵，列是不同样本
    vartotal = [givenParameters[1].stds[num_copies] .^ 2; givenParameters[2].stds[num_copies] .^ 2] # 2*2矩阵，行是不同片段，列是不同样本

    # varbetween_s1 = max.(10^-16, reshape(vartotal[1,:],1,:) .- varwithin) #组间方差，必然大于0
    # varbetween_s2 = max.(10^-16, reshape(vartotal[2,:],1,:) .- varwithin)

    varbetween_s1 = max.(10^-16, reshape(vartotal[1,:],1,:) .- varwithin) #1*2矩阵，不同组同属片段1的组间方差，必然大于0
    varbetween_s2 = max.(10^-16, reshape(vartotal[2,:],1,:) .- varwithin) #1*2矩阵，不同组同属片段2的组间方差，必然大于0   

    chi_dist_s1 = (meanwithin .- meantotal[1, :]) .^ 2 ./ (varwithin ./ pesudo_len + varbetween_s1) #1*2矩阵，不同组样本同属片段1的卡方值
    chi_dist_s2 = (meanwithin .- meantotal[2, :]) .^ 2 ./ (varwithin ./ pesudo_len + varbetween_s2) #1*2矩阵，不同组样本同属片段2的卡方值


    chi_dist = Chisq(2)       # 定义自由度为2的中心卡方分布
    p_value_g1 = 1 - cdf(chi_dist, sum(chi_dist_s1)) 
    p_value_g2 = 1 - cdf(chi_dist, sum(chi_dist_s2)) 

    output = (p_value_g1 * givenParameters[1].weights[num_copies[1]] + p_value_g2 * givenParameters[2].weights[num_copies[2]]) /
             (givenParameters[1].weights[num_copies[1]] + givenParameters[2].weights[num_copies[2]])

    return output
end

#*****************以下函数用于生成马尔可夫链******************************
function simulate_markov_chain(
    M::Matrix{Float64},  # 转移概率矩阵(n×n)
    C::Int,              # 时间步长
    initial_state::Int   # 初始状态(1到n)
)::Vector{Int}
    #此函数根据转移概率矩阵M，生成长度为C的向量
    #输出向量chain的元素为所在时刻的状态
    #输出状态位于1至n之间，用于CNV模拟时需要考虑是否矫正为 0至n-1
    #依赖包“StatsBase”
    
    # 参数校验
    @assert size(M,1) == size(M,2) "转移矩阵必须是方阵"
    @assert all(0 .<= M .<= 1) "概率值必须在[0,1]范围内"
    @assert all(sum(M, dims=2) .≈ 1) "每行概率和必须为1"
    @assert 1 <= initial_state <= size(M,1) "初始状态越界"
    
    states = Vector{Int}(undef, C)
    states[1] = initial_state 
    
    # 状态转移模拟
    for t in 2:C
        states[t] = sample(1:size(M,2), Weights(M[states[t-1], :]))
    end

    return states
end

#*****************以下函数用于生成干扰因素如GC含量，MAPPING错误等******************************
function generate_multinomial_samples(
    p::Vector{Float64},  # 概率向量(各事件概率)
    n_trials::Int,       # 试验次数
    sample_size::Int     # 生成样本数
)::Vector{Vector{Int}}
    #依赖包“Distributions”

    # 参数校验
    @assert all(0 .<= p .<= 1) "概率值必须在[0,1]范围内"
    @assert sum(p) ≈ 1 "概率向量和必须为1"
    @assert n_trials > 0 "试验次数必须为正整数"
    
    # 创建多项分布
    d = Multinomial(n_trials, p)
    
    # 生成样本(每个样本是长度为length(p)的向量)
    return [rand(d) for _ in 1:sample_size]
end

#****************************************************************************************
function estimate_gmm(data, n_components; means=0, max_iter=400, tol=1e-6)
    # 估计混合正态分布的参数
    # 此版本消除了数值计算精度不足带来的'除以0'错误

    # 初始化参数
    n = length(data)
    weights = fill(1/n_components, n_components)
    if length(means) != n_components
          means = quantile(data, range(0, 1, length=n_components+2)[2:end-1])
    end
    #means = quantile(data, range(0, 1, length=n_components+2)[2:end-1])
    stds = fill(std(data)/n_components, n_components)
    
    # EM算法
    for iter in 1:max_iter        
        # E步：计算后验概率
        responsibilities = zeros(n, n_components)
        for k in 1:n_components            
            responsibilities[:,k] = weights[k] * pdf.(Normal(means[k], stds[k]), data)
            temp_posi = findall(responsibilities[:,k] .< 10^-300) 
            responsibilities[temp_posi,k] .= 10^-300  #限制最小值，避免出现'除以0'错误
        end
        responsibilities ./= sum(responsibilities, dims=2) 
        # M步：更新参数
        new_means = zeros(n_components)
        new_stds = zeros(n_components)
        for k in 1:n_components
            resp_k = responsibilities[:,k]            
            new_means[k] = sum(resp_k .* data) / sum(resp_k)
            new_stds[k] = sqrt(sum(resp_k .* (data .- new_means[k]).^2) / sum(resp_k))
            weights[k] = mean(resp_k)
        end
        
        # 检查收敛
        if maximum(abs.(new_means - means)) < tol && 
           maximum(abs.(new_stds - stds)) < tol
            break
        end        

        if iter == max_iter
            max_dev_mean = maximum(abs.(new_means - means))
            max_dev_std = maximum(abs.(new_stds - stds))
            println("Maximum iteration reached at iter = $iter")
			println("$max_dev_mean")
			println("$max_dev_std")
			
        end
        
        means, stds = new_means, new_stds
    end
    
    return (means=means, stds=stds, weights=weights)
end

#************************计算左右两翼具有相同拷贝数的概率*************************************
function rate_coincident(obs_data, givenParameters)
    #计算两个观察值各自拷贝数的概率矩阵
    #obs_data是长度为2的一个向量，记录左右两侧的测序深度观察值
    # givenParameters 是一个元组，来自函数function estimate_gmm
   #依赖库 Distributions.jl
    
    #组分数
   #num_components = length(givenParameters.means)
   normal_obs_data = (obs_data' .-  givenParameters.means) ./ givenParameters.stds

   d = Normal(0, 1)       # 定义标准正态分布（μ=0, σ=1）
   #cdf_value = cdf(d, 1.96) 

   pdf_value = pdf(d, normal_obs_data) 

   pdr_value = (pdf_value .* givenParameters.weights) ./(sum(pdf_value .* givenParameters.weights, dims = 1))

   coincident_rate = sum(pdr_value[:, 1] .* pdr_value[:, 2])

   return coincident_rate
end

#***********************统计向量中非0值的游程***************************************
function run_length_nonzero(v)
    # 初始化结果矩阵
    result = Matrix{Int}(undef, 0, 3)
    n = length(v)
    n == 0 && return result
    
    # 初始化游程变量
    current_val = v[1]
    start_idx = 1
    
    for i in 2:n
        if v[i] != current_val
            # 当值变化时，记录前一个游程（如果是非零）
            if current_val != 0
                result = vcat(result, [current_val start_idx i-1])
            end
            current_val = v[i]
            start_idx = i
        end
    end
    
    # 处理最后一个游程
    if current_val != 0
        result = vcat(result, [current_val start_idx n])
    end
    
    return result
end

#***********************确定聚合片段的拷贝数***************************************
function identify_copy_number(frag_depth, givenParameters, threshold)
    # frag_depth 是一个向量，每片段的拷贝数
    # givenParameters 是估计得到的参数，是estimate_gmm()的返回值
    # threshold 为阈值，
    # 返回一个矩阵, 第一列为拷贝数（-1表示无法确定），随后列为取拷贝数1,2,3，...的概率

    #依赖库 Distributions.jl

    #组分数
    num_components = length(givenParameters.means)
    normal_obs_data = (frag_depth .-  givenParameters.means) ./ givenParameters.stds

    temp = zeros(num_components)

    d = Normal(0, 1)       # 定义标准正态分布（μ=0, σ=1）
    pdf_value = pdf(d, normal_obs_data) 
    pdr_value = (pdf_value .* givenParameters.weights) ./(sum(pdf_value .* givenParameters.weights, dims = 1))
    
    tempBayesFactors = pdr_value ./ (1 .- pdr_value)
    (max_val, max_idx) = findmax(tempBayesFactors)
    if max_val > threshold
        copyNum = max_idx
    else
        copyNum = -1
    end
    
    return [copyNum pdr_value']
end

#***********************拷贝数变化的检验****************************
function prob_copy_difference(frag_depth, num_copies, givenParameters)
    # frag_depth 是一个长度为2的向量，样本1和样本2该片段的测序深度
    # num_copies 是样本1和样本2该片段的各自拷贝数
    # givenParameters 是估计得到的参数，是样本1和样本2的estimate_gmm()的返回值
    # 返回，片段在两个样本中不存在拷贝数变化的概率

    #依赖库 Distributions.jl

    #组分数
    num_copies = Int64.(num_copies)
    normal_frag_depths_g1 = (frag_depth .-  givenParameters[1].means[num_copies[1]]) ./ givenParameters[1].stds[num_copies[1]]
    normal_frag_depths_g2 = (frag_depth .-  givenParameters[2].means[num_copies[2]]) ./ givenParameters[2].stds[num_copies[2]]

    chi_dist = Chisq(2)       # 定义自由度为2的中心卡方分布
    p_value_g1 = 1 - cdf(chi_dist, sum(normal_frag_depths_g1.^2)) 
    p_value_g2 = 1 - cdf(chi_dist, sum(normal_frag_depths_g2.^2)) 

    output = (p_value_g1 * givenParameters[1].weights[num_copies[1]] + p_value_g2 * givenParameters[2].weights[num_copies[2]]) /
             (givenParameters[1].weights[num_copies[1]] + givenParameters[2].weights[num_copies[2]])

    return output
end

