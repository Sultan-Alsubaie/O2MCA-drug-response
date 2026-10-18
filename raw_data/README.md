# raw_data

The four original CSV files (`Gene_expression.csv`, `Methylation.csv`, `Drug_response_IC50.csv`, `PubchemFP.csv`) are restored here automatically from `data/inputs.rds.part*` the first time any script runs. `01_read_and_align_data.R` checks them against `config/raw_input_hashes.json`.
