# ccc-skill — 细胞间通讯分析技能 🔗

[English](README_EN.md)

## 功能特性

- **工具选型决策树** —— 根据数据类型、分析目标、语言偏好自动推荐最佳 CCC 工具
- **5 大工具全覆盖** —— LIANA+、CellPhoneDB、CellChat、NicheNet、COMMOT 的完整工作流模板
- **Python + R 双语** —— AnnData/scanpy 和 Seurat 管线均有代码模板
- **空间转录组支持** —— Visium、MERFISH、Xenium 等空间 CCC 分析指导
- **多样本比较** —— tensor 分解、MOFA+、CellChat merged、MultiNicheNet 工作流
- **陷阱警告** —— 每个工具的常见错误和参数调优建议
- 当用户提出细胞间通讯、配体-受体分析、信号通路推断等需求时自动触发

## 多平台支持

| 平台 | 状态 | 说明 |
|------|------|------|
| **Claude Code** | ✅ 完全支持 | 原生 SKILL.md 加载 |
| **Codex** | ✅ 完全支持 | `agents/openai.yaml` 配置 |
| **OpenClaw** | ✅ 完全支持 | `metadata.openclaw` 命名空间 |

## 覆盖工具

| 工具 | 语言 | 核心场景 | 指南 |
|------|------|---------|------|
| **LIANA+** | Python | 多方法聚合、空间双变量、MISTy、tensor/MOFA+ | [guides/liana.md](guides/liana.md) |
| **CellPhoneDB** | Python | 高可信度人工策展 LR 数据库、置换检验 | [guides/cellphonedb.md](guides/cellphonedb.md) |
| **CellChat** | R | 信号通路层级、丰富可视化、条件比较 | [guides/cellchat.md](guides/cellchat.md) |
| **NicheNet** | R | 配体→TF→靶基因预测、MultiNicheNet 多样本 | [guides/nichenet.md](guides/nichenet.md) |
| **COMMOT** | Python | 最优传输空间 CCC、信号方向向量场 | [guides/commot.md](guides/commot.md) |

## 对比：有技能 vs 无技能

| 能力 | 原生 Agent | 本技能 |
|------|-----------|--------|
| 知道哪个工具适合哪种数据/目标 | 大致了解 | ✅ 精确决策树 |
| 正确的 API 调用和参数 | 经常出错 | ✅ 验证过的代码模板 |
| 知道常见陷阱（z-scale、counts_data 默认值等） | ❌ | ✅ 每个工具 5-10 条警告 |
| 多工具组合工作流 | ❌ | ✅ 完整管线模板 |
| 空间 CCC 参数调优（dis_thr、bandwidth 等） | ❌ | ✅ 按技术平台推荐 |

## 前置条件

Python 工具：
```bash
pip install liana cellphonedb commot scanpy
```

R 工具：
```r
devtools::install_github("jinworks/CellChat")
devtools::install_github("saeyslab/nichenetr")
# 多样本:
devtools::install_github("saeyslab/multinichenetr")
```

## 技能安装

### Claude Code

```bash
# 全局安装
git clone https://github.com/Agents365-ai/ccc-skill.git ~/.claude/skills/ccc-skill

# 项目级安装
git clone https://github.com/Agents365-ai/ccc-skill.git .claude/skills/ccc-skill
```

### Codex

```bash
git clone https://github.com/Agents365-ai/ccc-skill.git ~/.codex/skills/ccc-skill
```

### OpenClaw

```bash
git clone https://github.com/Agents365-ai/ccc-skill.git ~/.openclaw/skills/ccc-skill

# 项目级
git clone https://github.com/Agents365-ai/ccc-skill.git skills/ccc-skill
```

### 安装路径汇总

| 平台 | 全局路径 | 项目路径 |
|------|---------|---------|
| Claude Code | `~/.claude/skills/ccc-skill/` | `.claude/skills/ccc-skill/` |
| Codex | `~/.codex/skills/ccc-skill/` | N/A |
| OpenClaw | `~/.openclaw/skills/ccc-skill/` | `skills/ccc-skill/` |

## 使用方式

直接用自然语言描述你的 CCC 分析需求即可：

```
> 我有 Visium 空间转录组数据，想分析细胞间通讯，推荐用什么工具？

> 帮我写 LIANA+ rank_aggregate 的代码，数据是 scRNA-seq h5ad 格式

> CellChat 怎么比较两个条件之间的通讯差异？

> 我想预测哪些配体导致了受体细胞的转录变化，该用什么方法？

> COMMOT 的 dis_thr 参数怎么设置？我的数据是 MERFISH
```

技能会根据决策树选择合适的工具，参考对应的 guide 生成正确代码。

## 文件结构

```
ccc-skill/
├── SKILL.md              # 技能主文件 — 决策树 + 工具对比 + 通用原则
├── guides/
│   ├── liana.md          # LIANA+ 完整工作流指南
│   ├── cellphonedb.md    # CellPhoneDB 完整工作流指南
│   ├── cellchat.md       # CellChat 完整工作流指南
│   ├── nichenet.md       # NicheNet 完整工作流指南
│   └── commot.md         # COMMOT 完整工作流指南
├── agents/
│   └── openai.yaml       # Codex/OpenAI 平台配置
├── README.md             # 中文文档（本文件）
└── README_EN.md          # English documentation
```

- `SKILL.md` —— **核心文件**，所有平台以此为技能指令
- `guides/*.md` —— 各工具的详细代码模板和参数说明，SKILL.md 按需引用

## 常见问题

### LLM 已经知道这些工具了，为什么还需要 skill？

LLM 对这些工具有大致了解，但在以下方面经常出错：

1. **参数默认值** —— 如 CellPhoneDB 的 `counts_data` 默认是 `'ensembl'` 而非 symbol，LIANA+ 的 `use_raw=True` 会读 `adata.raw`
2. **工具选型** —— 不同数据类型/分析目标该用哪个工具，LLM 的推荐常不精确
3. **空间参数** —— `dis_thr`、`bandwidth`、`contact.range` 等需要根据技术平台调整
4. **API 变更** —— CellChat v2 vs v1、Seurat v5 vs v4 的 API 差异
5. **隐性陷阱** —— z-scale 破坏零值、数值型细胞名、missing `subsetData()` 等

Skill 将验证过的代码模板和陷阱警告前置注入到 agent 的上下文中，避免试错循环。

## 已知限制

- 代码模板基于各工具 2025-2026 年版本，API 可能随版本更新变化
- R 工具的代码模板使用 Seurat 管线；SCE 用户需自行适配
- 不包含所有参数的完整 API 文档，聚焦于常用工作流和关键参数

## License

MIT

## 支持

如果这个技能对你有帮助，欢迎打赏支持作者：

<table>
  <tr>
    <td align="center">
      <img src="https://raw.githubusercontent.com/Agents365-ai/images_payment/main/qrcode/wechat-pay.png" width="180" alt="微信支付">
      <br>
      <b>微信支付</b>
    </td>
    <td align="center">
      <img src="https://raw.githubusercontent.com/Agents365-ai/images_payment/main/qrcode/alipay.png" width="180" alt="支付宝">
      <br>
      <b>支付宝</b>
    </td>
    <td align="center">
      <img src="https://raw.githubusercontent.com/Agents365-ai/images_payment/main/qrcode/buymeacoffee.png" width="180" alt="Buy Me a Coffee">
      <br>
      <b>Buy Me a Coffee</b>
    </td>
  </tr>
</table>

## 作者

**Agents365-ai**

- Bilibili: https://space.bilibili.com/441831884
- GitHub: https://github.com/Agents365-ai
