# 甘肃地质灾害与典型农村房屋损伤风险研究：代码与配套数据

本仓库对应最新版正文 1022 及其补充材料，包含原有研究程序、数值结果、
有限元输入文件和新增的复现核查脚本，不包含未发表的完整论文、私人照片、
账号凭据或 Abaqus 商业软件。

## 运行

在 Python 3.10 及以上环境中，从仓库根目录执行：

```bash
python -m pip install -e .
python scripts/audit_archive.py
python scripts/reproduce_annual_example.py --output results/annual
python scripts/reproduce_fragility_bounds.py --output results/fragility
python scripts/check_local_sensitivity.py
```

年度算例的解析 EAL 为 4974.742362 元/年，20 万年随机模拟均值为
4989.45 元/年。脚本保留原随机种子、参数范围、积分和损失状态。
三个减灾方案降低年度损失，但在当前成本与折现条件下 BCR 均小于 1。

## 与正文对应

D1：区域空间验证、地形因素贡献和 ALE；D2：源区—建筑条件到达；
D3：有限元响应及 117 个预设位置的质量状态；D4：容量区间和条件易损性；
D5：局部敏感性；D6：独立参数化年度风险算例；D7：三地完整市域的情景空间表达。
详细对应关系见 docs/MANUSCRIPT_MAPPING.md。

区域模型重训、2884 万条传播路径重算仍需原始 DEM、编目、网格数组和建筑轮廓。
这些原始第三方数据不在本次归档内，输入结构及来源说明见 docs/INPUTS.md。
完整有限元求解需自行获得 Abaqus 许可。归档的 50 份 INP 是合格案例输入，
一个输入可包含多个强度步骤，不能把 INP 数量等同于 117 个预设工况位置。

## 研究边界

- 天水区域易发性及传播为区域计算；三灾种有限元为九类原型的代表性荷载情景。
- 传播速度代理与九类房屋实际荷载的转换接口仍需工程标定。
- D_t、EDI 是局部材料或构件损伤指标，容量离散性不是现场实测标定结果。
- 地图中的建筑类型为情景赋值，未知真实属性保持未知；NC 单独保留。
- 年度 EAL 是独立算例，不赋给地图中的实际建筑。

代码采用 MIT 许可；自有非空间数值数据采用 CC BY 4.0；OSM/Overture 相关
数据库衍生内容保留 ODbL 及来源署名。详见 DATA_LICENSES.md。
DOI 只有在 Zenodo 正式发布并可访问后才写入引用，不以预留号码代替发布。
