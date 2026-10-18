# Validation

The final internal models (13 learned paths and two baselines) and the four external routes were refitted from the prepared inputs with the recorded configurations. The new predictions matched the archived predictions exactly (maximum absolute difference 0).

| File | Content |
|---|---|
| `internal_refit_comparison.csv` | Refitted versus archived internal predictions |
| `external_expression_refit_comparison.csv`, `external_reduced_refit_comparison.csv` | Refitted versus archived external predictions |
| `manuscript_numeric_checks.csv` | 113 values in Tables 2, 3, 5, 6 and 7 recalculated and compared with the manuscript |
| `session_info.txt` | R session used for the check |

The full candidate grid of the cross-validation was not rerun for this check.
