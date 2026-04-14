# ccc-skill — Cell-Cell Communication Analysis Skill 🔗

[中文文档](README_CN.md)

## What it does

- **Tool selection decision tree** — automatically recommends the best CCC tool based on data type, analysis goal, and language preference
- **8 tools with full guides** — LIANA+, CellPhoneDB, CellChat, NicheNet, COMMOT, Squidpy, MEBOCOST, DIALOGUE
- **12 additional tools referenced** — CellOracle, SCENIC+, stLearn, MultiNicheNet, Scriabin, FlowSig, cell2cell, SpaTalk, NICHES, CrossTalkeR, NeuronChat, scDiffCom
- **Python + R** — code templates for both AnnData/scanpy and Seurat pipelines
- **Spatial transcriptomics** — guidance for Visium, MERFISH, Xenium, and other spatial CCC analysis
- **Beyond protein LR** — metabolite-mediated CCC (MEBOCOST) and multicellular programs (DIALOGUE)
- **Multi-sample comparison** — tensor decomposition, MOFA+, CellChat merged, MultiNicheNet workflows
- **Pitfall warnings** — common mistakes and parameter tuning advice per tool
- Triggers automatically when user asks about cell-cell communication, ligand-receptor analysis, or signaling inference

## Multi-Platform Support

| Platform | Status | Details |
|----------|--------|---------|
| **Claude Code** | ✅ Full support | Native SKILL.md loading |
| **Codex** | ✅ Full support | `agents/openai.yaml` configuration |

## Tools Covered

| Tool | Language | Core Use Case | Guide |
|------|----------|---------------|-------|
| **LIANA+** | Python | Multi-method aggregation, spatial bivariate, MISTy, tensor/MOFA+ | [guides/liana.md](guides/liana.md) |
| **CellPhoneDB** | Python | Curated LR database, permutation testing, v5 scoring | [guides/cellphonedb.md](guides/cellphonedb.md) |
| **CellChat** | R | Pathway hierarchy, rich visualization, condition comparison | [guides/cellchat.md](guides/cellchat.md) |
| **NicheNet** | R | Ligand→TF→target prediction, MultiNicheNet multi-sample | [guides/nichenet.md](guides/nichenet.md) |
| **COMMOT** | Python | Optimal transport spatial CCC, signaling direction vector fields | [guides/commot.md](guides/commot.md) |
| **Squidpy** | Python | scverse ecosystem LR analysis, OmniPath database, zero-friction | [guides/squidpy.md](guides/squidpy.md) |
| **MEBOCOST** | Python | Metabolite-mediated CCC (non-protein signals) | [guides/mebocost.md](guides/mebocost.md) |
| **DIALOGUE** | R | Multicellular programs (cross-cell-type coordination) | [guides/dialogue.md](guides/dialogue.md) |

## Comparison: With skill vs. Without

| Capability | Native Agent | This Skill |
|-----------|-------------|------------|
| Knows which tool fits which data/goal | Roughly | ✅ Precise decision tree |
| Correct API calls and parameters | Often wrong | ✅ Verified code templates |
| Knows common pitfalls (z-scale, counts_data default, etc.) | ❌ | ✅ 5-10 warnings per tool |
| Multi-tool combo workflows | ❌ | ✅ Complete pipeline templates |
| Spatial CCC parameter tuning (dis_thr, bandwidth, etc.) | ❌ | ✅ Per-technology recommendations |

## Prerequisites

Python tools:
```bash
pip install liana cellphonedb commot scanpy
```

R tools:
```r
devtools::install_github("jinworks/CellChat")
devtools::install_github("saeyslab/nichenetr")
# Multi-sample:
devtools::install_github("saeyslab/multinichenetr")
```

## Skill Installation

### Claude Code

```bash
# Global (available in all projects)
git clone https://github.com/Agents365-ai/ccc-skill.git ~/.claude/skills/ccc-skill

# Project-level
git clone https://github.com/Agents365-ai/ccc-skill.git .claude/skills/ccc-skill
```

### Codex

```bash
git clone https://github.com/Agents365-ai/ccc-skill.git ~/.codex/skills/ccc-skill
```

### Installation paths summary

| Platform | Global path | Project path |
|----------|-------------|--------------|
| Claude Code | `~/.claude/skills/ccc-skill/` | `.claude/skills/ccc-skill/` |
| Codex | `~/.codex/skills/ccc-skill/` | N/A |

## Usage

Just describe your CCC analysis needs in natural language:

```
> I have Visium spatial transcriptomics data and want to analyze cell-cell communication. What tool should I use?

> Write LIANA+ rank_aggregate code for my scRNA-seq h5ad data

> How do I compare communication differences between two conditions in CellChat?

> I want to predict which ligands caused transcriptional changes in receiver cells. What method should I use?

> How should I set COMMOT's dis_thr for MERFISH data?
```

The skill picks the right tool from the decision tree and generates correct code from the corresponding guide.

## File Structure

```
ccc-skill/
├── SKILL.md              # Main skill file — decision tree + tool comparison + universal principles
├── guides/
│   ├── liana.md          # LIANA+ complete workflow guide
│   ├── cellphonedb.md    # CellPhoneDB complete workflow guide
│   ├── cellchat.md       # CellChat complete workflow guide
│   ├── nichenet.md       # NicheNet complete workflow guide
│   ├── commot.md         # COMMOT complete workflow guide
│   ├── squidpy.md        # Squidpy CCC complete workflow guide
│   ├── mebocost.md       # MEBOCOST metabolite CCC guide
│   └── dialogue.md       # DIALOGUE multicellular programs guide
├── agents/
│   └── openai.yaml       # Codex/OpenAI platform configuration
├── README.md             # English documentation (this file)
└── README_CN.md          # Chinese documentation
```

- `SKILL.md` — **core file**, loaded by all platforms as skill instructions
- `guides/*.md` — detailed code templates and parameter reference per tool, referenced by SKILL.md on demand

## FAQ

### LLMs already know about these tools. Why do I need a skill?

LLMs have general awareness of these tools but frequently get details wrong:

1. **Parameter defaults** — e.g., CellPhoneDB's `counts_data` defaults to `'ensembl'` not symbol; LIANA+'s `use_raw=True` reads from `adata.raw`
2. **Tool selection** — which tool for which data type/goal is often imprecise without guidance
3. **Spatial parameters** — `dis_thr`, `bandwidth`, `contact.range` must be tuned per technology platform
4. **API changes** — CellChat v2 vs v1, Seurat v5 vs v4 API differences
5. **Hidden pitfalls** — z-scaling destroys zeros, numeric cell type names, missing `subsetData()`, etc.

The skill injects verified code templates and pitfall warnings into the agent's context, avoiding trial-and-error cycles.

## Known Limitations

- Code templates are based on 2025-2026 tool versions; APIs may change with updates
- R code templates use the Seurat pipeline; SCE users may need adaptation
- Does not include full API documentation for every parameter — focuses on common workflows and critical parameters

## License

MIT

## Support

If this skill helps you, consider supporting the author:

<table>
  <tr>
    <td align="center">
      <img src="https://raw.githubusercontent.com/Agents365-ai/images_payment/main/qrcode/wechat-pay.png" width="180" alt="WeChat Pay">
      <br>
      <b>WeChat Pay</b>
    </td>
    <td align="center">
      <img src="https://raw.githubusercontent.com/Agents365-ai/images_payment/main/qrcode/alipay.png" width="180" alt="Alipay">
      <br>
      <b>Alipay</b>
    </td>
    <td align="center">
      <img src="https://raw.githubusercontent.com/Agents365-ai/images_payment/main/qrcode/buymeacoffee.png" width="180" alt="Buy Me a Coffee">
      <br>
      <b>Buy Me a Coffee</b>
    </td>
  </tr>
</table>

## Author

**Agents365-ai**

- Bilibili: https://space.bilibili.com/441831884
- GitHub: https://github.com/Agents365-ai
