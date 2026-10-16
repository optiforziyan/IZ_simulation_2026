"""Independent raster checks and result-grounded manuscript materials.
Run with Python 3 plus numpy, pandas and Pillow. No optimization is rerun here.
"""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
D, R, M = ROOT/'data', ROOT/'results', ROOT/'manuscript'
M.mkdir(exist_ok=True)
def raster(p): return np.asarray(Image.open(p), dtype=float).ravel()
S = np.column_stack([raster(D/'SDM'/f'Sp{j}_presence.tif') for j in range(1,9)])
P = pd.read_csv(D/'population_presence_matrix.csv').to_numpy()
F = pd.read_csv(D/'functional_presence_matrix.csv').to_numpy()
B = pd.read_csv(D/'branch_presence_matrix.csv').to_numpy()
bl = pd.read_csv(D/'branch_species_incidence.csv')['length'].to_numpy()
pm = pd.read_csv(D/'population_parameters.csv')
sp = pd.read_csv(D/'species_parameters.csv')
local = np.column_stack([raster(D/f'{k}.tif') for k in ['TD','PD','FD','GD']])
lock = raster(D/'Locked_in.tif')
assert np.array_equal(S.sum(1),local[:,0])
assert np.allclose(B@bl,local[:,1])
assert np.allclose(F.sum(1),local[:,2])
assert np.allclose(P@pm.pi.to_numpy(),local[:,3])
assert np.allclose(P.sum(0),pm.range_cells)
for j in range(8):
    assert np.array_equal(P[:,pm.species.eq(f'Sp{j+1}')].sum(1),S[:,j])

def components(v):
    a=v.reshape(24,24).astype(bool); seen=set(); n=0
    for i,j in zip(*np.where(a)):
        if (i,j) in seen: continue
        n+=1; stack=[(i,j)]; seen.add((i,j))
        while stack:
            y,x=stack.pop()
            for z,w in [(y-1,x),(y+1,x),(y,x-1),(y,x+1)]:
                if 0<=z<24 and 0<=w<24 and a[z,w] and (z,w) not in seen:
                    seen.add((z,w)); stack.append((z,w))
    return n
clusters=[components(S[:,j]) for j in range(8)]
assert clusters == [1]*8
table=pd.read_csv(R/'main_comparison.csv').set_index('scenario')
Z=(local[:,1:]-local[:,1:].min(0))/(local[:,1:].max(0)-local[:,1:].min(0))
for name,row in table.iterrows():
    x=raster(R/f'{name}.tif'); assert set(x).issubset({0,1})
    assert np.all(x[lock==1]==1)
    assert x.sum()==row.area
    assert np.allclose(x@local[:,1:],row[['local_PD_sum','local_FD_sum','local_GD_sum']].to_numpy(float))
    assert np.isclose(x@Z.mean(1),row.Q_sum)
    assert ((S*x[:,None]).sum(0)>=sp.range_cells.to_numpy()*sp.target_fraction.to_numpy()-1e-7).sum()==row.species_targets_met
    assert ((P*x[:,None]).sum(0)>=.2*P.sum(0)-1e-7).sum()==row.population_targets_met
    assert bl[(B*x[:,None]).sum(0)>0].sum()==row.network_PD
    assert ((F*x[:,None]).sum(0)>0).sum()==row.network_FD
pop_rows=[]
for name in table.index:
    x=raster(R/f'{name}.tif')
    for k,pop in pm.iterrows():
        held=float(x@P[:,k]); required=float(.2*P[:,k].sum())
        pop_rows.append({'scenario':name,'population':pop.population,'species':pop.species,
                         'range_cells':pop.range_cells,'target_cells':int(np.ceil(required-1e-9)),
                         'held_cells':int(held),'target_met':held>=required-1e-7})
pd.DataFrame(pop_rows).to_csv(R/'population_representation.csv',index=False)
agx=raster(R/'Aggregated_features_20pct.tif')
agg=pd.DataFrame({'feature':['TD','PD','FD','GD_proxy'],'landscape_total':local.sum(0),
                  'target_amount':.2*local.sum(0),'held_amount':agx@local})
agg['target_met']=agg.held_amount>=agg.target_amount-1e-6
assert agg.target_met.all()
agg.to_csv(R/'aggregate_feature_targets.csv',index=False)
equivalence=[]
for rho in [.25,1,4]:
    tag=f'{rho:g}'
    a=raster(R/f'Effective_cost_rho{tag}.tif'); b=raster(R/f'Soft_penalty_rho{tag}.tif')
    costs=1+rho*(1-Z.mean(1))
    equivalence.append({'rho':rho,'effective_objective':a@costs,'penalty_objective':b.sum()+rho*(b@(1-Z.mean(1))),
                        'objective_difference':a@costs-b@costs,'different_cells':int(np.sum(a!=b))})
eq=pd.DataFrame(equivalence); assert (eq.objective_difference.abs()<1e-5).all()
eq.to_csv(R/'cost_penalty_equivalence.csv',index=False)
solver=pd.read_csv(R/'solver_audit.csv')
assert solver.status.eq('Optimal').all() and solver.gap.eq(0).all()
rep=pd.read_csv(R/'replicate_results.csv')
piv=rep.pivot(index='seed',columns='scenario',values='Q_sum')
gain0=100*(piv.h0/piv.b-1); gain1=100*(piv.h1/piv.h0-1)
rep_summary=rep.groupby('scenario')[['area','species_targets_met','species_present','population_targets_met','Q_sum']].agg(['min','mean','max'])
rep_summary.to_csv(R/'replicate_summary.csv')
pd.DataFrame({'seed':piv.index,'Q_gain_same_area_percent':gain0,'Q_gain_extra_area_percent':gain1}).to_csv(R/'paired_gains.csv',index=False)
pd.DataFrame(local,columns=['TD','PD','FD','GD']).corr().to_csv(R/'local_index_correlations.csv')
audit={'raster_solutions_verified':len(table),'optimizations_verified':len(solver),'solver_status':'all Optimal','maximum_mip_gap':float(solver.gap.max()),'main_species_connected_components':clusters,'seed_sequence':[20260925,20260934],'all_checks_passed':True}
(R/'validation.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')

def fmt_table(df):
    headers=list(df.columns)
    rows=['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']
    for _,r in df.iterrows(): rows.append('| '+' | '.join(str(v) for v in r)+' |')
    return '\n'.join(rows)

base=table.loc['Species_min_area']; best=table.loc['Species_Q_max_at_Cstar']; h10=table.loc['Hier_equal_eps0.1']
gains=100*(best[['local_PD_sum','local_FD_sum','local_GD_sum']]/base[['local_PD_sum','local_FD_sum','local_GD_sum']]-1)
keys=['Aggregated_features_20pct','Aggregated_score_fixed_area','Species_min_area','Species_Q_max_at_Cstar','Hier_equal_eps0.1','Species_and_populations_Q']
names=['Aggregated features with 20% targets','Aggregated score at equal area','Species minimum area','Species plus Q at minimum area','Species plus Q at 10% tolerance','Species and population targets plus Q']
display=table.loc[keys,['area','species_targets_met','population_targets_met','Q_sum']].copy()
display.insert(0,'Scenario',names)
display.columns=['Scenario','Area cells','Species targets met / 8','Population targets met / 19','Local score Q']
display['Local score Q']=display['Local score Q'].map(lambda x:f'{x:.3f}')
for col in display.columns[1:4]: display[col]=display[col].astype(int)
main_table=fmt_table(display)

readme=f'''# 模拟重建与审稿回应

本套材料将“实体代表性”和“局地多维空间偏好”分开检验。所有结果来自运行后的文件，未筛选有利随机种子。原始文件保留在原目录，新材料全部位于 simulation_v2。

## 结论与推荐

推荐把真实成本保持为面积，将多维指数作为明确的次级空间偏好；主结果采用两阶段优化，而不是人为改变预算口径。先在物种硬约束下求最低面积 C*，再在 C(x) <= (1+epsilon)C* 下最大化局地加权得分 Q。epsilon=0 是严格不增加面积的比较；0.05、0.10、0.20 是展示代价的假设情景，不是推荐给实际保护项目的阈值。

成本修正与正向软惩罚在本例中数学等价：

`min sum[c_i + rho * sum_d w_d(1-z_id)] x_i`

等价于

`min sum c_i*x_i + rho * sum_i sum_d w_d(1-z_id)*x_i`。

二者使用相同物种硬约束。rho=0.25、1、4 的成对求解均核对了最优目标值一致。它们是同一模型的不同实现，不能写成一方表现优于另一方。有效成本是偏好修正的目标系数，不能当作实际货币开支，也不能用其总量重新定义另一套“30%面积预算”。

另一个重要区别：当选取面积可变时，`sum(1-Q_i)*x_i` 与 `-sum(Q_i*x_i)` 并不等价，前者还惩罚选取单元数。只有固定面积时，两者仅差常数。

## 已完成的模拟

- 24 × 24 共576个等面积单元，8个物种，3个广布种、3个中间范围物种、2个狭域特有种。主数据所有物种均为单一连通分布斑块。
- 广布种分布占42.7%—52.1%；两个特有种占2.95%和2.43%。同时输出连续适宜度和二值范围，未把合成适宜度称为实测拟合SDM。
- 物种约束全部为硬约束；现有保护区在同一景观的全部情景中一致。所有模拟成本都是每个单元1个面积单位。
- 从同一套物种、树、功能实体和19个种群分区生成局地TD、PD、FD和GD代理指标，并另行评价不重复的网络分支、功能实体及种群目标。
- 5组权重 × 4个epsilon、3组rho成对比较、3组物种目标政策，以及10个预先确定的连续随机种子。
- 检验最低面积最优解的Q值上下界，避免只选一个较差的基准解来夸大增益。
- {len(solver)}个保存的优化问题均返回Optimal、MIP gap=0；独立Python复核了{len(table)}个主结果栅格及其评价指标。

## 主要结果

{main_table}

在相同110个单元下，相对于一个最低面积基准解，最优局地偏好解的PD总和提高{gains.iloc[0]:.2f}%，FD总和变化{gains.iloc[1]:.2f}%，GD代理总和提高{gains.iloc[2]:.2f}%；Q提高{100*(best.Q_sum/base.Q_sum-1):.2f}%。这不是相对于所有最优基准解的固定收益：110单元可行解的Q范围是65.062—86.209，偏好优化选择其上界。

但种群目标达标数从16/19下降到15/19，种群出现数从19下降到18。增加10%面积后，种群目标仍为15/19。只有显式加入种群硬目标，才保证19/19达标，需要118个单元，比物种最低面积多7.27%。因此不能把GD代理指数总和的增长解释为种群或遗传谱系互补性的改善。

10个景观中，物种硬约束情景全部8/8达标；同面积聚合得分情景只有4—5个物种达标，并仅出现6—7个物种。同面积局地偏好相对于一个求解器选取的最低面积基准的Q增益平均{gain0.mean():.2f}%，范围{gain0.min():.2f}%—{gain0.max():.2f}%。这些是同一合成景观族内的描述性结果，不是经验普适效应或统计显著性检验。

## 文件入口

- `code/run_simulation.R`：完整生成、优化、评价、敏感性、重复模拟和绘图脚本。替代旧planA/planB的分析逻辑。
- `code/audit_results.R`：核对全部缓存求解的约束、目标值、预算和gap。
- `code/validate_and_document.py`：独立读取TIF验证主结果，生成本文档与英文材料。
- `data/SDM/Sp1_SDM.tif` 至 `Sp8_SDM.tif`：可直接用于规划的二值范围；同目录有presence别名和continuous suitability图层。
- `results/main_comparison.csv`：主要情景汇总；`species_representation.csv`：每个物种的目标与实际保留量。
- `results/weight_epsilon_sensitivity.csv`、`target_sensitivity.csv`、`replicate_results.csv`：全部敏感性和重复结果。
- `results/checkpoints/`：每次求解的输入、解、状态和目标值，支持中断恢复。
- `figures/00_conceptual_framework.svg`：含FD和假设数值说明的新概念图，可编辑矢量版；同名PNG为3600×2160像素。
- `figures/01_species_ranges.png` 至 `04_tradeoffs.png`：物种分布、局地指数、规划网络和代价图。
- `manuscript/Revised_Main_Text.md`：统一数学模型、模拟结论与适用边界的英文重写稿。
- `manuscript/Supplementary_Methods_and_Results.md`：可用于补充材料的完整方法与结果。
- `manuscript/Response_to_Reviewer.md`：8条审稿意见的英文逐条回复，按本套重写稿定位，尚未转换为原Word文件的修订痕迹或投稿页码。

## 运行

安装R以及terra、prioritizr、highs。当前实际使用的版本见results/sessionInfo.txt。在项目根目录运行：

```powershell
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' simulation_v2/code/run_simulation.R all
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' simulation_v2/code/audit_results.R
& 'C:/Program Files/R/R-4.6.1/bin/Rscript.exe' simulation_v2/code/make_concept_figure.R
```

最后用安装有numpy、pandas、Pillow的Python执行`simulation_v2/code/validate_and_document.py`生成独立核验和文稿。脚本无D盘硬编码，目录相对于脚本自动定位；主种子20260925，重复种子20260925—20260934。相同输入优先读取已验证的缓存。若需完全重新求解，可把checkpoints目录改名后重跑。

## 如何回答审稿人

1. 明确认可指数作为features的合法性，并用真正的聚合指数目标情景演示“满足指数目标不保证物种目标”。
2. 区分局地加总与网络不重复实体；报告种群下降这一反例，不宣称全面改善互补性。
3. 将贡献定位为概念与建模角色区分，承认最小集、线性惩罚和epsilon方法为既有优化思路。
4. 分开解释相对权重w、总体偏好强度rho、真实成本容忍度epsilon，并报告敏感性；所有数值均为透明的示例设定。
5. 提供可运行模拟、输入图层、完整结果、重复景观及局限。
6. 新概念图加入FD与假设值注释；正文公式统一编号；首次出现SDM给出全称。

## 不应作出的结论

本例没有实测遗传数据、真实种群过程、气候、连通性、物种存活概率或现实经济成本；FD是离散功能实体数，GD是模拟pi的局地加和代理。所有具有物种硬目标的网络都出现8个物种，因此网络分支并集PD固定为16、功能实体并集固定为4。不能把局地PD/FD加总变化写成网络独特PD/FD增加，也不能把这套模拟写成已经证明某方案在自然界普遍优越。
'''
(ROOT/'README_中文.md').write_text(readme,encoding='utf-8')

supplement=r'''# Supplementary methods and results

## S1 Purpose and scope

We use a fully synthetic, reproducible example to distinguish explicit representation requirements from local biodiversity preferences. The experiment illustrates consequences of feature definition; it is not an empirical validation of a new optimizer or a universally superior conservation strategy. No landscapes or runs were screened to retain favourable outcomes.

## S2 Synthetic species distributions

The landscape contains 24 by 24 square planning units, each assigned cost one. Coordinates span 0 to 24,000 in each direction and are schematic; the GeoTIFF projection supplies a common Cartesian grid, not a real study location. Eight species comprise three broad-ranging species, three intermediate-scale species and two range-restricted endemics. For species j, a rotated elliptical Gaussian suitability surface is multiplied by a low-amplitude, spatially smooth sinusoidal perturbation and rescaled to a maximum of one. Centres are shifted by uniform draws within ±0.025 landscape units of predefined centres; rotation angles are sampled between -0.5 and 0.5 radians. The Gaussian widths and baseline centres are in data/species_parameters.csv and the complete draw sequence is in the script. Binary occurrence is defined by normalized suitability >= 0.30. These surfaces are synthetic habitat suitability, not fitted species distribution models (SDMs). Binary occurrences enter optimization; continuous surfaces are also distributed.

The primary landscape has species range sizes of 296, 300, 246, 74, 93, 53, 17 and 14 cells. Each binary range is a single four-neighbour connected component. Main targets are 20% for species occupying at least 25% of the landscape, 30% for those occupying 10% to less than 25%, and 50% for those occupying less than 10%. Targets therefore depend on realized range size, not the descriptive broad/intermediate/endemic labels. Because decisions and presence are binary, effective required cell counts are the ceilings of target fraction times range size: 60, 60, 50, 23, 28, 27, 9 and 7. These are illustrative policy assumptions, not validated persistence thresholds. We also test uniform 30% targets and an increase of 10 percentage points in every range-based target.

Two compact circular reserves centred at (0.43, 0.48) and (0.69, 0.60), with radii 0.085 and 0.055, are locked in. The same reserves are used by all scenarios. Their areas count toward total cost and representation; reported area is total protected area, not additional area beyond existing reserves. No connectivity or boundary-length objective is imposed.

## S3 Entities and local biodiversity indices

TD is local species richness. A fixed rooted toy phylogeny contains eight terminal branches of length one, four sister-pair branches of length one and two four-species clade branches of length two. Local PD is the sum of lengths of branches with at least one descendant species present locally. Total possible distinct branch length is 16. Tree topology, lengths and branch incidence are distributed.

FD is the number of four discrete functional entities represented locally. The species-by-entity incidence matrix is distributed. This transparent richness measure is not functional dispersion, functional volume, or a general estimate of ecosystem function.

Each of the three broad-ranging species is partitioned into three spatial population strata; each other species is partitioned into two, giving 19 strata. Within each species, occupied cells are ordered along a randomly oriented spatial gradient and partitioned into nearly equal-sized groups. Every occupied cell belongs to exactly one population of that species. Each stratum receives a simulated pi value sampled uniformly from 0.002 to 0.020. The local GD proxy is the sum of simulated pi values for the populations present in a cell. It is an illustrative additive proxy; no sequences, alleles, demography or independently estimated genetic lineages are simulated. Different strata represent separate planning entities, not empirically established evolutionarily significant units.

PD, FD and GD are each min-max normalized over the complete landscape. Let z_id denote these normalized local values, w_d non-negative weights summing to one, and q_i = sum_d w_d z_id. The secondary network score Q is sum_i q_i x_i. Local indices can be correlated because they share species inputs; their correlations are provided in results/local_index_correlations.csv. Equal weighting gives equal coefficients after normalization, not equal independent information or ecological importance.

Separate network evaluations count each phylogenetic branch and functional entity once, and report both population presence and the number of populations with at least 20% of their occupied cells selected. Local summed GD is not interpreted as distinct genetic diversity retained. Likewise, sums of local PD and FD are not network union measures.

## S4 Optimization scenarios

Binary x_i specifies selection and r_ij is binary species occurrence. Species targets are hard constraints sum_i r_ij x_i >= alpha_j sum_i r_ij in every species-feature scenario. Solver feasibility is checked after every run.

1. Species minimum area: minimize C = sum_i x_i subject to species targets and lock-in. The resulting minimum is C*. There can be multiple optimum networks.
2. Minimum-area envelope: separately minimize and maximize Q with C = C*, retaining species targets and lock-in. This measures alternative-optimum variability instead of treating one baseline network as unique.
3. Aggregated features: minimize area with four features TD, PD, FD and GD, each required to achieve 20% of its landscape sum, plus lock-in. Species targets are intentionally absent. This is a valid index-target problem, with a different feasible set; its smaller area is not interpreted as a fair efficiency improvement over species-target planning.
4. Aggregated-score diagnostic: maximize Q with area exactly C* and lock-in, without species targets. The implementation supplies zero targets on the four index features. This diagnostic isolates entity omission at equal area and is not described as the same formulation as scenario 3.
5. Effective-cost formulation: minimize sum_i [1 + rho(1-q_i)] x_i subject to species targets and lock-in.
6. Soft-penalty formulation: minimize sum_i x_i + rho sum_i (1-q_i)x_i under identical constraints. Scenarios 5 and 6 are algebraically equivalent. We solve both at rho = 0.25, 1 and 4 and check equality of optimal objective values. Identical objective values do not generally require identical maps when alternative optima exist.
7. Hierarchical cost-tolerance formulation: first obtain C*, then maximize Q with species targets, lock-in and sum_i x_i <= floor((1+epsilon)C*). We use epsilon = 0, 0.05, 0.10 and 0.20. This is cost priority followed by a weighted local-index objective; it is not a claim to estimate the complete multidimensional Pareto frontier. At epsilon=0, any gain reflects selection among minimum-area networks. At positive epsilon, gains can also reflect additional area and are reported with that cost.
8. Explicit population targets: retain species constraints and add a hard 20% target for every population stratum, first minimizing area, then maximizing Q at that minimum area. This illustrates the extra requirements needed to guarantee the simulated within-species representation targets; it is not matched to the species-only feasible set.

We implement hard-target scenarios using prioritizr::add_min_set_objective or add_min_penalties_objective with explicit targets; no minimum-shortfall objective is used. Secondary maximization uses a negative linear penalty on q_i. A positive penalty on 1-q_i has an additional selection-count term, so it is not generally equivalent to maximizing q_i when total area varies. Real cost is always reported as selected area, regardless of objective coefficients.

## S5 Parameter choices and sensitivity

Equal weights (1/3, 1/3, 1/3) are the neutral illustrative reference. Additional sets are PD-focused (0.6, 0.2, 0.2), FD-focused (0.2, 0.6, 0.2), GD-focused (0.2, 0.2, 0.6), and (0.3, 0.3, 0.4), matching the earlier example. All five sets are crossed with four epsilon values. These coefficients are disclosed decision preferences, not fitted biological parameters. The overall penalty scale rho is kept separate from relative weights; in the original notation, lambda_d = rho w_d when conventional cost per unit is one. Rescaling physical costs would require reconsidering rho but would not change the meaning of a relative cost tolerance epsilon.

The three rho values explore a 16-fold range in strength; they are not chosen after examining gains. The epsilon grid explicitly shows no-cost, 5%, 10% and 20% additional-cost policies, with integer rounding. None is asserted to be universally appropriate. Target-policy sensitivity is reported independently.

The main seed is 20260925. Ten predeclared consecutive seeds, 20260925 through 20260934, repeat species minimum area, species-constrained Q maximization at epsilon=0 and 0.10, and equal-area aggregate-score maximization. The main landscape is included in these ten; it is not an additional independent replicate. Repeats vary species centres, smooth perturbations, rotations, population partition orientations and simulated pi values, while retaining the landscape family, tree, functional incidence, widths and reserve configuration. Results are descriptive robustness checks within this family, not independent empirical tests of generality.

## S6 Results

MAIN_TABLE

The species minimum area is 110/576 cells (19.10%). The Q envelope over minimum-area species-feasible networks is 65.062 to 86.209; the solver-returned baseline is 83.185. Thus the baseline is not a uniquely defined comparator. At the same 110 cells, the Q-optimal network increases local PD sum from 961 to 995, leaves local FD sum at 357, and increases the local GD proxy sum from 4.64064 to 5.02117. Its Q gain over the displayed baseline is 3.64%, but its population target count declines from 16 to 15, and one of the 19 populations is absent. The positive-cost-tolerance solutions retain all species targets but do not restore all population targets: at epsilon=0.10, 121 cells are selected and 15/19 population targets are achieved.

The aggregated-feature model achieves all four index targets using 71 cells, but only three species targets, omitting both endemic species entirely. The equal-area aggregate-score diagnostic achieves five species targets at 110 cells; Sp7 receives 8 of the required 9 cells and Sp8 is absent. This demonstrates a possible mismatch between index goals and entity goals, not an impossibility of using aggregated features.

Explicit species-plus-population requirements need at least 118 cells, 7.27% above species-only C*. At that area, Q optimization retains all eight species targets and all 19 population targets. The comparison with 121-cell species-only planning shows why changing the feature set can matter more for omitted entities than increasing the area allocated to a local proxy objective.

All species-target networks contain all eight species, and therefore share network union PD = 16 and union FD = 4. We do not claim their larger summed local indices increase unique branch or functional-entity coverage. The three pairs of effective-cost and soft-penalty problems have equal optimal objectives; in this main dataset all three strengths return the same reported evaluation metrics. The modest weight changes mostly have small effects here, and are not evidence that weights never matter in other landscapes.

REPLICATE_TEXT

## S7 Reproducibility and limitations

All saved optimization runs returned status Optimal and relative MIP gap zero. Single-thread HiGHS uses random_seed=1, a 120-second limit and zero requested gap; the script refuses to accept a positive reported gap. The completed run did not rely on a time-limited incumbent. Checkpoint files retain exact inputs, outputs, statuses and objectives. audit_results.R recomputes target satisfaction, lock-in, area limits and objective values. An independent Python audit reads the exported GeoTIFFs and feature matrices and verifies main result summaries and connected ranges. Package and runtime versions are in sessionInfo.txt.

The example omits real economic costs, uncertain SDM predictions, persistence processes, spatial connectivity requirements and environmental change. Endemic species are spatially separated in this landscape family by construction, so failures of aggregate-only selection are illustrative counterexamples, not estimates of how often failure occurs in nature. Target thresholds, population strata, tree lengths, functional entities and pi values are synthetic assumptions. Claims about robust management performance require broader landscape families and empirical evaluation. No optimization of an additive local score alone guarantees representation of unmodelled entities.

## Figure captions

Figure S1 (01_species_ranges.png). Synthetic binary distributions of eight species, with range sizes and hard representation targets. Background is unoccupied and teal is occupied. Each range is connected; the species overlap is intentional. Continuous synthetic suitability is supplied separately.

Figure S2 (02_local_indices.png). Local TD, toy-tree PD, functional-entity richness and the summed simulated-pi GD proxy. Panels use separate scales. These are local summaries, not unique network-level quantities.

Figure S3 (03_network_comparison.png). Six illustrative solutions. Light grey is unselected, teal newly selected and navy locked-in existing reserve. Titles report total selected cell count. Scenarios have different feature requirements as specified in S4; only the indicated fixed-area comparisons hold area constant.

Figure S4 (04_tradeoffs.png). Equal-weight secondary local score and population target counts across epsilon = 0, 0.05, 0.10 and 0.20. The first point is the best local score among minimum-area networks, not the arbitrary minimum-area baseline. Additional area raises the optimized score but does not guarantee previously omitted population targets.

## Software documentation

The implementation follows the documented hard-target and linear-penalty semantics of prioritizr:
- https://prioritizr.net/reference/add_min_set_objective.html
- https://prioritizr.net/reference/add_linear_penalties.html
- https://prioritizr.net/reference/add_min_penalties_objective.html
- https://prioritizr.net/reference/add_linear_constraints.html
'''
reptext=f'''Across ten synthetic landscapes, species minimum area ranges from 104 to 110 cells. All species-constrained solutions satisfy 8/8 species targets. Equal-area aggregate-score solutions meet only 4-5 targets and contain 6-7 species. The paired Q gain at unchanged area averages {gain0.mean():.2f}% (range {gain0.min():.2f}-{gain0.max():.2f}%) relative to each solver-returned minimum-area baseline. Allowing 10% extra cost increases Q by a further mean {gain1.mean():.2f}% (range {gain1.min():.2f}-{gain1.max():.2f}%) relative to the Q-optimal minimum-area network. Because integer budgets are floored, realized area increments can be slightly below 10%. Mean population-target counts are 15.1 for the arbitrary minimum-area baseline, 13.9 for Q optimization at minimum area, and 14.3 for Q optimization at 10% tolerance. Thus local-score improvements are not a reliable substitute for explicit population representation in these simulations.'''
supplement=supplement.replace('MAIN_TABLE',main_table).replace('REPLICATE_TEXT',reptext)
(M/'Supplementary_Methods_and_Results.md').write_text(supplement,encoding='utf-8')

main=r'''# From Multidimensional Biodiversity Patterns to Conservation Features in Systematic Conservation Planning

The Kunming-Montreal Global Biodiversity Framework has intensified interest in where conservation actions should be deployed and which biodiversity components they should represent. Systematic conservation planning (SCP) addresses these questions through explicit features, representation targets, costs and implementation constraints (Margules & Pressey, 2000; Kukkala & Moilanen, 2013). Increasing attention to phylogenetic (PD), functional (FD) and genetic (GD) dimensions complements taxonomic diversity (TD), but raises a modelling question: which aspects of biodiversity should define representation requirements, and which should influence spatial preferences?

Aggregated biodiversity indices can serve as conservation features when clearly defined, spatially explicit and assigned appropriate targets. Their inclusion is not intrinsically incompatible with SCP. However, achieving a target for an aggregated index does not necessarily ensure representation of the biological entities that generate it. Two sites can have the same local species richness, PD, FD or GD while differing in species, branches, functional entities or genetic lineages (Figure 1). Consequently, identical local values need not imply identical contributions to a network that already contains some of those entities. Evidence that local plant PD hotspots differ from areas selected to maximize cumulative PD illustrates this distinction (Tietje et al., 2023); it does not imply that every index-based optimization is invalid.

Complementarity is conditional on the entities, targets and other sites considered. Individual species distributions with explicit targets permit optimization to account for representation across the selected network. The same principle can be extended to phylogenetic branches, functional entities, populations or genetic units when suitable data exist (Carvalho et al., 2017; Paz-Vinas et al., 2018). Species-only targets do not automatically guarantee these additional dimensions: a network can meet every species target while omitting a population or genetic lineage within a species. The appropriate feature set therefore follows the conservation objective and available evidence.

Where only aggregated multidimensional surfaces are available, they can be retained as additional features or used to guide spatial preferences alongside explicit species targets. One established option is a linear penalty added to conventional planning-unit costs. Let x_i be binary selection, c_i the conventional cost, and z_id a normalized local value of dimension d in planning unit i. The latter is not the unit's marginal contribution to network-level diversity. Let w_d >= 0 sum to one, and rho >= 0 denote the overall preference strength. A minimum-set formulation is

\[
\min_x \left[\sum_i c_i x_i + \rho\sum_i\sum_d w_d(1-z_{id})x_i\right], \tag{1}
\]

subject to explicit representation requirements

\[
\sum_i r_{ij}x_i\geq\alpha_j\sum_i r_{ij},\quad\forall j,\qquad x_i\in\{0,1\}. \tag{2}
\]

In species distribution model (SDM)-based applications, r_ij can represent predicted suitable habitat amount for species j; alpha_j is the required proportion. Other implementation constraints, including existing reserves, can be retained. This formulation builds on established minimum-set and penalty-based planning (Ball et al., 2009; Watts et al., 2009), rather than introducing a new optimization method.

The same objective can be implemented using effective costs

\[
c_i^{\mathrm{eff}}=c_i+\rho\sum_d w_d(1-z_{id}). \tag{3}
\]

Thus effective-cost modification and a corresponding additive soft penalty are algebraically equivalent when applied to the same objective and constraints. Effective costs encode spatial preferences; they should not be interpreted as measured financial costs or used to redefine a supposedly unchanged physical-area budget. In the notation lambda_d, lambda_d = rho w_d. Separating rho from w_d distinguishes overall strength from relative dimensional priorities. Neither has a universal value: normalization, cost units and planning preferences determine their interpretation. Furthermore, the positive term (1-z_id)x_i penalizes selecting a unit as well as rewarding its index value relative to alternatives; it is not generally interchangeable with a negative benefit term -z_id x_i when selected area varies.

A transparent alternative gives conventional cost explicit priority. First find C*, the minimum conventional cost satisfying Eq. (2), then optimize a declared secondary objective under a cost tolerance:

\[
C^*=\min_{x\in\mathcal F}\sum_i c_i x_i,\qquad
\max_{x\in\mathcal F}\sum_i\sum_d w_dz_{id}x_i
\quad\mathrm{subject\ to}\quad\sum_i c_i x_i\leq(1+\varepsilon)C^*, \tag{4}
\]

where F contains the hard representation and implementation constraints. Epsilon expresses an acceptable proportional increase in conventional cost, not a biological constant. Epsilon=0 selects among minimum-cost alternatives; positive values explicitly permit a cost-benefit trade-off. This two-stage example still uses a weighted local-index objective. It makes the cost allowance explicit but does not eliminate value judgements or estimate the complete trade-off surface among biodiversity dimensions. We recommend reporting the chosen weights, cost tolerance and sensitivity, and testing actual entity representation separately.

To illustrate these distinctions, we constructed a reproducible 576-unit landscape with eight spatially clustered synthetic species ranges, including broad-ranging species and two restricted endemics (Supplementary Methods S1-S5; Figures S1-S3). Local indices were derived from a toy phylogeny, functional entities and 19 simulated population strata. Species targets were hard constraints. At the minimum species-feasible area of 110 units, optimizing the local index score increased summed local PD by 3.54%, left summed local FD unchanged, and increased a local GD proxy by 8.20% relative to one minimum-area solution. However, the number of population targets achieved declined from 16/19 to 15/19. These gains therefore cannot be interpreted as improved genetic complementarity. All species-target networks contained the same eight species and had identical unique network PD and functional-entity coverage.

An index-feature model that met 20% targets for each of four local indices satisfied only three of eight species targets and omitted both endemics; an equal-area aggregate-score comparison satisfied five species targets. Conversely, explicitly adding population targets achieved all 19 population targets at 118 units, 7.27% above the species-only minimum. Ten predeclared simulated landscapes, five weight sets, four cost tolerances, three penalty strengths and three species-target policies provide sensitivity checks (Supplementary Methods S5-S7). The examples show what can happen under specified assumptions; they do not establish that an index approach must fail, that weights never matter, or that local-index gains translate into unique network diversity.

The contribution of this framework is a conceptual separation of modelling roles. Features and targets specify the entities or quantities whose representation is required; costs and preferences influence where those requirements are achieved. Aggregated indices remain legitimate features for objectives defined in those terms. Where the intended objective is representation of distinct biological entities, those entities should be modelled explicitly whenever possible. Local multidimensional preferences can supplement this representation, but cannot guarantee it for entities absent from the feature set. This distinction supports transparent SCP while identifying the data and evaluation needed for empirical application.

![Conceptual framework](../figures/00_conceptual_framework.png)

Figure 1. Equal local multidimensional scores can conceal different biological identities. Units A and B share the same entities, whereas C contains different species and additional branches, functional entities and genetic lineages. Given A has already been selected, equal aggregated values alone cannot distinguish the additional entity representation of B and C. Aggregated indices can be valid conservation features; meeting their targets does not necessarily meet targets for the underlying entities. Hard constraints guarantee only the representation requirements explicitly included. All SR, PD, FD and GD values shown are hypothetical illustrations, not empirical estimates. SR denotes species richness.

## References

Ball, I. R., H. P. Possingham, and M. E. Watts. 2009. Marxan and Relatives: Software for Spatial Conservation Prioritization. In A. Moilanen, K. A. Wilson, and H. P. Possingham, eds., Spatial Conservation Prioritization, pp. 185-195. Oxford University Press.

Carvalho, S. B., G. Velo-Anton, P. Tarroso, et al. 2017. Spatial Conservation Prioritization of Biodiversity Spanning the Evolutionary Continuum. Nature Ecology & Evolution 1: 0151.

Kukkala, A. S., and A. Moilanen. 2013. Core Concepts of Spatial Prioritisation in Systematic Conservation Planning. Biological Reviews 88: 443-464.

Margules, C. R., and R. L. Pressey. 2000. Systematic Conservation Planning. Nature 405: 243-253.

Paz-Vinas, I., G. Loot, V. Hermoso, C. Veyssiere, N. Poulet, G. Grenouillet, and S. Blanchet. 2018. Systematic Conservation Planning for Intraspecific Genetic Diversity. Proceedings of the Royal Society B: Biological Sciences 285: 20172746.

Tietje, M., A. Antonelli, F. Forest, R. Govaerts, S. A. Smith, M. Sun, W. J. Baker, and W. L. Eiserhardt. 2023. Global Hotspots of Plant Phylogenetic Diversity. New Phytologist 240: 1636-1646.

Watts, M. E., I. R. Ball, R. S. Stewart, C. J. Klein, K. Wilson, C. Steinback, R. Lourival, L. Kircher, and H. P. Possingham. 2009. Marxan with Zones: Software for optimal conservation based land- and sea-use zoning. Environmental Modelling & Software 24: 1513-1521.
'''
(M/'Revised_Main_Text.md').write_text(main,encoding='utf-8')

response=r'''# Response to Reviewer 1

We thank the reviewer for identifying the need to distinguish local biodiversity summaries from representation of biological entities. We have clarified the scope and conceptual contribution, aligned the numerical implementation with hard representation constraints, and added a reproducible synthetic example, sensitivity analyses and explicit limitations. The example also shows that improving local multidimensional scores does not necessarily improve population representation; we report this result rather than interpret index gains as network-level complementarity.

## Major comment 1

**Comment summary:** Distinguish the limitations of aggregated indices for entity complementarity from the claim that they cannot be conservation features.

**Response:** We agree. The second paragraph of the revised main text explicitly states that clearly defined, spatially explicit indices with suitable targets can serve as conservation features. Our concern is the mismatch that can arise between a target for a scalar index and representation of the entities generating it, rather than the mathematical admissibility of the index. We removed the categorical hotspot-regression claim and qualified the concluding recommendation. We also implemented a valid aggregated-feature optimization with explicit 20% targets for TD, PD, FD and a GD proxy. It meets these index targets but, in the illustrative landscape, meets only three of eight species targets. We distinguish this different-target example from the equal-area aggregate-score diagnostic and do not treat its smaller area as superior efficiency. Figure 1 now states this distinction directly.

**Location:** Main text paragraphs beginning “Aggregated biodiversity indices” and “An index-feature model”; Figure 1; Supplement S4 and S6.

## Major comment 2

**Comment summary:** z_id is a local value, not marginal contribution to the selected network; explain how complementarity is addressed.

**Response:** We now explicitly define z_id as a normalized local index. Complementarity operates through representation of the entities specified in the feature matrix and their network-wide targets, while the index term supplies a spatial preference. Species targets alone do not guarantee within-species genetic or population representation. We therefore separately evaluate summed local indices, unique branch and functional-entity coverage, and population representation. At equal area in the main example, optimizing the local score raises the GD proxy sum but reduces population-target attainment from 16/19 to 15/19. All species-target solutions share the same network union PD and functional-entity coverage. Adding explicit population targets achieves 19/19 at a separately reported area cost. This identifies the scope of the representation guarantee and avoids attributing network-level complementarity to the local score.

**Location:** Main text third paragraph and paragraphs introducing and interpreting the example; Supplement S3, S4, S6 and S7; Tables in main_comparison.csv and species_representation.csv.

## Major comment 3

**Comment summary:** Clarify the methodological basis and contribution beyond existing optimization approaches.

**Response:** We cite established minimum-set and penalty-based planning formulations (Ball et al., 2009; Watts et al., 2009) and explicitly state that we do not introduce a new optimizer. The contribution is the conceptual separation of representation requirements from local spatial preferences. We also clarify that effective-cost modification and the corresponding additive penalty are algebraically identical under the same objective and constraints. Numerical paired checks at three strengths confirm equal optimal objective values. The revised Eq. (4) is a conventional two-stage cost-tolerance formulation, presented as a transparent decision structure rather than an algorithmic innovation. We no longer imply that this weighted secondary objective eliminates subjective choices or recovers an entire multidimensional Pareto frontier.

**Location:** Main text surrounding Eqs. (1), (3) and (4), and final paragraph; Supplement S4 and software documentation.

## Major comment 4

**Comment summary:** Justify the selection of lambda_d and epsilon.

**Response:** We separate lambda_d into overall strength rho and relative weights w_d, so lambda_d = rho w_d. This makes their roles and dependence on cost units explicit. Equal weights are used as an illustrative reference, with PD-, FD- and GD-focused weights and the previous 0.3/0.3/0.4 set as sensitivity cases. Rho values 0.25, 1 and 4 span a 16-fold range. Epsilon values 0, 0.05, 0.10 and 0.20 correspond to explicitly stated additional-area allowances above the minimum feasible cost, including a no-additional-cost reference. We also vary species target policies. These are transparent scenario assumptions, not biologically calibrated constants or universal recommendations. Their purpose is to show sensitivity and make decision preferences explicit. We acknowledge the relatively small weight effects in this particular landscape rather than infer general parameter insensitivity.

**Location:** Main text explanations of Eqs. (3)-(4); Supplement S5; weight_epsilon_sensitivity.csv and target_sensitivity.csv.

## Major comment 5

**Comment summary:** Add numerical examples or practical evaluation and clarify application scope.

**Response:** We added a reproducible synthetic landscape with 576 planning units, eight clustered species distributions including two restricted endemics, a toy phylogeny, four functional entities, and 19 population strata. All scenarios retain the same reserve lock-in and physical-area cost definition. Species-feature scenarios enforce all species targets as hard constraints. Comparisons include actual aggregated-feature targets, an equal-area aggregate-score diagnostic, species-only planning, equivalent effective-cost/soft-penalty formulations, a cost-tolerance analysis and explicit population targets. We evaluate alternative minimum-area solutions, perform parameter sensitivity, and repeat four core scenarios over ten predeclared seeds. All 79 saved optimization problems return Optimal with zero reported MIP gap; numerical and independent raster audits verify the outputs. We emphasize that these examples demonstrate mechanisms under controlled assumptions and do not validate a universal empirical advantage. The source scripts, rasters, complete results, software versions and verification records accompany the supplement.

**Location:** Main text numerical-example paragraphs; Supplement S1-S7; Figures S1-S4; code, data and results directories.

## Minor comment 1

**Comment summary:** Mark figure values as illustrative, consider including FD, and improve figure resolution.

**Response:** We replaced the conceptual diagram with a figure that explicitly includes SR, PD, FD and GD and labels all values as hypothetical rather than empirical. Both the figure and caption contain this qualification. The figure is supplied as editable SVG and a 3600 by 2160 pixel PNG. Its wording also reflects the corrected distinction between valid index features and entity-representation guarantees.

**Location:** Revised Figure 1 and its caption; figures/00_conceptual_framework.svg and .png.

## Minor comment 2

**Comment summary:** Number the equations.

**Response:** The main text now consistently numbers the penalized objective, hard representation constraints, effective-cost equivalence and two-stage cost-tolerance formulation as Eqs. (1)-(4), and refers to them in the accompanying explanation.

## Minor comment 3

**Comment summary:** Define SDM at its first occurrence.

**Response:** We now write “species distribution model (SDM)-based applications” at first occurrence in the main text. The supplement also distinguishes synthetic suitability surfaces from fitted SDMs, avoiding the implication that the synthetic distributions are empirical predictions.
'''
(M/'Response_to_Reviewer.md').write_text(response,encoding='utf-8')
files=[p for p in ROOT.rglob('*') if p.is_file() and p.suffix not in {'.log'} and p.name!='file_checksums.sha256']
(R/'file_checksums.sha256').write_text('\n'.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(ROOT)) for p in sorted(files)),encoding='utf-8')
print(json.dumps(audit,indent=2))
print('Wrote Chinese report and three English manuscript materials.')
