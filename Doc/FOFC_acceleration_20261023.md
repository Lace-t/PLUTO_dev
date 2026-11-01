# FOFC 加速代码同步记录（2026-10-23）

本地源码 `/home/titan/Desktop/BHAC/PLUTO` 合入已在 GRE/r5_dc 验证的 P2 受影响行段回放，FOFC 默认最大重试次数设为 2。此前已同步的受影响 stencil 缓存、candidate 恢复、通量缓存回放时跳过 States 等改动继续保留。

## P2 行段回放

每次重试交换 mask 后，根据累计 stage_affected 计算各方向每行的最小和最大受影响索引。空行跳过 RHS 和守恒更新，非空行仅装载、计算包围段及所需相邻面；每个受影响单元仍按原方向顺序获得全部方向和局部源项贡献。

完整 stage U/V 恢复、缓存波速装载和全域 CFL 扫描继续执行。仅非 CT、DIV_CLEANING 路径启用；预先存在 FLAG_FLAT 或外层低阶回退时使用完整行回放。用户源项回调需为局部、无副作用函数；GRE/r5_dc 使用只依赖位置的引力。如果自定义回调依赖全行调用次数或副作用，应关闭该开关。

编译开关（均可在 definitions.h 或编译参数中覆盖）：

```c
#define PCP_STAGE_SPARSE_REPLAY YES /* NO 恢复完整行回放 */
#define PCP_STAGE_LOCAL_PASSES 2    /* 设为 4 可恢复此前最大重试次数 */
```

最大重试次数是数值算法参数。2 次重试更早进入 PressureFix，会改变其处理集合、dt 和后续轨迹；它通过了该算例的 200 步检查，不代表任意问题均与 4 次等价。

## 服务器验证证据

- P2：同一 allocation，5 节点/320 MPI，从 restart #2 推进 1687–1886 共 200 步，排除初始化和 restart。平均 5.654518 → 4.601422 秒/步，耗时降低 18.624%。逐步 time/dt、事件计数和最终 320 rank 的 U/V/flags 哈希一致，最终非法单元为零。这是一次固定顺序对照。
- 最大重试 4 → 2：同规模同起始 restart 的 200 步结果为 4.601422 → 3.861410 秒/步，耗时降低 16.082%。来自不同 allocation，数值轨迹不同；显式阶段 PressureFix 失败、外层重试和最终非法单元均为零。不能直接视为等物理时长加速。
- P1 归约合并、CheckStage 异步实验尚无整步加速证据，本次未合入；性能诊断及 rank 到达时间插桩也未合入。

原始报告位于 baseline/9_gre/fofc_p2_20261023/REPORT_zh.md 和 fofc_p2_pass2_20261023/REPORT_zh.md。本次备份、补丁、校验值和本地检查日志位于 baseline/9_gre/fofc_local_sync_20261023/；本地检查结果见其中 REPORT_zh.md。
