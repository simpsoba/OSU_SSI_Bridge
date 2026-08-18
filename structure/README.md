# structure/

Pier, deck, pile cap, and piles. Each piece is two scripts:

| | When | What |
|---|---|---|
| `Build*Nodes.tcl` | `BuildModel.tcl`, before soil | nodes and nodal mass |
| `Build*Elements.tcl` | after soil gravity, or immediately in `PlotModel.tcl` | beam-columns / ZLS |

`BuildStructElements.tcl` sources the four `*Elements` files (pier, deck, cap, piles). Soil gravity needs the pile/cap **nodes** for p-y ties, but not beam-column stiffness (those wait until `FoldStructNodes.tcl`).

Shared nodes (same tag, stacked mass, no `equalDOF`): cap TC = pier base; deck soffit CL = pier top; pile heads = cap BL / BC / BR.

Sections: `PierSection.tcl`, `BuildPileSection.tcl`. Helper: `IncrMass.tcl`.
